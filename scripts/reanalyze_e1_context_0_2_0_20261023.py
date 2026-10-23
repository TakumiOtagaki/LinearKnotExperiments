#!/usr/bin/env python3
"""Audit all profile-only tasks and reuse the adopted E1 analysis definitions.

Only context columns are replaced. Non-context metrics/resources retain their
certified saved values. The existing stratified and local-error scripts read
the original BPP/accessibility and a checked new-profile overlay. No RNA CLI
is executed. Large intermediates live in one temporary node-local SQLite DB.
"""
import argparse
from collections import Counter, defaultdict
import csv
from functools import lru_cache
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import pickle
import sqlite3
import subprocess
import sys
import tarfile
import time
import zlib

import analyze_accuracy as accuracy
import analyze_e1_production as production
import analyze_e1_pk_stratified as stratified
import analyze_e1_local_errors_89baea52 as local
import pseudoloop_context_array_20261023 as runner
from run_cli import digest, now
from task_storage import ArchiveReader, MemoryFile, RunReader, local_workspace

ROOT = Path(__file__).resolve().parents[1]
OLD = {c: ROOT/'results/local'/n for c, n in runner.OLD.items()}
SUFFIX = '0.2.0-20261023'
NAMES = {c: f'e1-bprna-{c}-context-{SUFFIX}' for c in OLD}
OVERLAY = None
CONNECTIONS = {}
ORIGINAL_READ_GROUP = production.read_group


def cache_archive_reads():
    """Read each task tar once; share the frozen input index across tasks.

    Artifact checksums still pass through RunReader.artifact. A task archive is
    at most 156 MB in the short cohort and is freed with its RunReader.
    """
    def initialize(self, path):
        self.path = Path(path)
        self.members = {}
        self._bytes = self.path.read_bytes()
        with tarfile.open(fileobj=io.BytesIO(self._bytes), mode='r:') as tf:
            for member in tf:
                if member.isfile():
                    assert member.name not in self.members
                    self.members[member.name] = (member.offset_data, member.size)
                tf.members.clear()

    def read(self, name):
        offset, size = self.members[str(Path(name))]
        data = self._bytes[offset:offset+size]
        assert len(data) == size
        return data

    ArchiveReader.__init__ = initialize
    ArchiveReader.read = read
    original = RunReader.archive

    @lru_cache(maxsize=2)
    def inputs(path):
        return ArchiveReader(path)

    def archive(self, path):
        path = Path(path).resolve()
        return inputs(path) if path.name == 'inputs.tar' else original(self, path)

    RunReader.archive = archive


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def overlay_read_group(reader, jobs, name, mode, record, provenance):
    parsed, check = ORIGINAL_READ_GROUP(reader, jobs, name, mode, record, provenance)
    pid = os.getpid()
    if pid not in CONNECTIONS:
        CONNECTIONS[pid] = sqlite3.connect(f'file:{OVERLAY}?mode=ro', uri=True)
    row = CONNECTIONS[pid].execute('SELECT value FROM profiles WHERE name=? AND mode=?', (name, mode)).fetchone()
    assert row is not None
    parsed['LinearPKCapR'] = pickle.loads(zlib.decompress(row[0]))
    # Keep the original saved consistency measurements, but independently check
    # the combined old BPP/accessibility and new profile as well.
    assert accuracy.check_consistency(name, mode, 0, parsed)['status'] == 'success'
    return parsed, check


def read_metrics(directory):
    with (directory/'metrics.tsv').open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        rows = {(r['input'], r['mode']): r for r in reader}
    return fields, rows


def audit_and_metrics(plan_path, output, db):
    p = runner.load(plan_path)
    plan_hash = digest(plan_path)
    records = {}
    for c, a in p['input_archives'].items():
        names = {v['member']: k for k, v in p['inputs'].items() if v['cohort'] == c}
        with tarfile.open(ROOT/a['path'], 'r|') as tf:
            for member in tf:
                if member.name in names:
                    data = tf.extractfile(member).read()
                    name = names[member.name]
                    assert hashlib.sha256(data).hexdigest() == p['inputs'][name]['sha256']
                    lines = data.decode().splitlines()
                    records[name] = (name, lines[1], lines[2])
    assert set(records) == set(p['inputs'])
    saved = {c: read_metrics(d) for c, d in OLD.items()}
    baseline_hashes = {c: digest(d/'metrics.tsv') for c, d in OLD.items()}
    for c, d in OLD.items():
        assert baseline_hashes[c] == json.loads((d/'provenance.json').read_text())['metrics_sha256']
    changed = defaultdict(int)
    maxima = Counter()
    audited = []
    @lru_cache(maxsize=3)
    def old_archive(path):
        return ArchiveReader(ROOT/path)
    db.execute('CREATE TABLE profiles (name TEXT, mode TEXT, value BLOB, PRIMARY KEY(name,mode))')
    total = 0
    started = time.monotonic()
    for number, task in enumerate(p['tasks'], 1):
        d = plan_path.parent/'tasks'/task['tier']/f"{task['task_id']:06d}"
        meta = json.loads((d/'task.json').read_text())
        assert meta['status'] == 'success' and not meta['pilot']
        assert meta['plan_sha256'] == plan_hash and meta['runner_sha256'] == p['script_sha256']
        assert meta['inputs'] == task['inputs'] and meta['completed_jobs'] == task['profile_jobs']
        data = (d/'task_results.jsonl').read_bytes()
        assert data.endswith(b'\n')
        rows = [json.loads(line) for line in data.splitlines()]
        expected = {(name, m['name']) for name in task['inputs'] for m in p['modes'][p['inputs'][name]['cohort']]}
        assert len(rows) == len(expected) == task['profile_jobs']
        assert {(r['input'], r['mode']) for r in rows} == expected
        assert meta['status_counts'] == {'success': len(rows)}
        members = {}
        by_name = defaultdict(dict)
        end = 0
        for r in rows:
            name, mode = r['input'], r['mode']
            inp = p['inputs'][name]
            assert r['status'] == 'success' and r['returncode'] == 0 and r['error'] is None
            assert r['input_sha256'] == inp['sha256'] and r['binary_sha256'] == p['binary_sha256']
            assert r['parameter_sha256'] == p['parameter_sha256'] and not r['warmup'] and r['repetition'] == 0
            assert r['program'] == 'LinearPKCapR' and r['subcommand'] == 'profile'
            args = next(m['args'] for m in p['modes'][inp['cohort']] if m['name'] == mode)
            assert r['command'][1:3] == ['profile', '-i'] and r['command'][4:] == args
            assert Path(r['command'][3]).name == name+'.txt'
            assert r['archive_end'] > end and r['archive_end'] % 512 == 0
            end = r['archive_end']
            assert set(r['artifacts']) == {r['id']+'/'+f for f in ('stdout.txt', 'stderr.txt', 'command.json', 'time.tsv')}
            for m, info in r['artifacts'].items():
                assert m not in members
                members[m] = (r, info)
        seen = set()
        with tarfile.open(d/'task_outputs.tar', 'r|', bufsize=1024*1024) as tf:
            for m in tf:
                if not m.isfile():
                    continue
                assert m.name in members and m.name not in seen
                seen.add(m.name)
                r, info = members[m.name]
                raw = tf.extractfile(m).read()
                assert len(raw) == info['bytes'] and hashlib.sha256(raw).hexdigest() == info['sha256']
                if m.name.endswith('/command.json'):
                    assert json.loads(raw) == r['command']
                if not m.name.endswith('/stdout.txt'):
                    continue
                assert info['sha256'] == r['stdout_summary']['sha256']
                name, mode = r['input'], r['mode']
                record = records[name]
                profile = accuracy.read_contexts(MemoryFile(raw), record)
                old = next(o for o in p['inputs'][name]['observations'] if o['mode'] == mode)
                old_raw = old_archive(old['archive']).read(old['member'])
                assert hashlib.sha256(old_raw).hexdigest() == old['stdout_sha256']
                previous = accuracy.read_contexts(MemoryFile(old_raw), record)
                for c in accuracy.CONTEXTS:
                    if c in ('Multiloop', 'PseudoLoop'):
                        continue
                    error = max(abs(a-b) for a, b in zip(profile[c], previous[c]))
                    maxima[c] = max(maxima[c], error)
                    assert error <= 1e-9
                error = max(abs(a+b-c-e) for a,b,c,e in zip(profile['Multiloop'],profile['PseudoLoop'],previous['Multiloop'],previous['PseudoLoop']))
                maxima['Multiloop+PseudoLoop'] = max(maxima['Multiloop+PseudoLoop'], error)
                assert error <= 1e-9
                by_name[name][mode] = profile
                db.execute('INSERT INTO profiles VALUES (?,?,?)', (name, mode, zlib.compress(pickle.dumps(profile, protocol=5))))
        assert seen == set(members) and (d/'task_outputs.tar').stat().st_size >= end
        for name, profiles in by_name.items():
            cohort = p['inputs'][name]['cohort']
            reference = profiles['exact' if cohort == 'short' else 'beam400_outer']
            record = records[name]
            for mode, profile in profiles.items():
                if mode == 'exact':
                    continue
                row = saved[cohort][1][name, mode]
                assert row['input_sha256'] == p['inputs'][name]['sha256']
                values = accuracy.compare_contexts(profile, reference, record)
                update = {'context_tv_mean_Gfree': values['context_tv_mean_Gfree'],
                          **{'context_'+c+'_mae_Gfree': v for c,v in values['context_mae_Gfree'].items()}}
                for key, value in update.items():
                    new = 'NA' if value is None else str(value)
                    changed[key] += row[key] != new
                    row[key] = new
        db.commit()
        total += len(rows)
        audited.append(dict(tier=task['tier'], task_id=task['task_id'], jobs=len(rows),
                            checkpoint_sha256=digest(d/'task_results.jsonl'),
                            task_metadata_sha256=digest(d/'task.json'),
                            archive_sha256=digest(d/'task_outputs.tar')))
        if number % 10 == 0 or number == len(p['tasks']):
            print(f'profile audit {number}/{len(p["tasks"])} tasks; {total} profiles; {time.monotonic()-started:.1f}s', flush=True)
    assert total == p['profile_jobs'] == 131753
    audit = dict(all_passed=True, tasks=len(audited), profiles=total, inputs=len(records),
                 failed=0, missing=0, per_base_invariance_tolerance=1e-9,
                 per_base_invariance_max_absolute_differences=dict(maxima),
                 plan_sha256=plan_hash, source_commit=p['source_commit'], tasks_audited=audited)
    write_json(output/'task_output_validation.json', audit)
    for c, (fields, rows) in saved.items():
        dest = output/NAMES[c]
        dest.mkdir()
        stratified.write_tsv(dest/'metrics.tsv', fields, rows.values())
        _, prior = read_metrics(OLD[c])
        non_context = [k for k in fields if not k.startswith('context_')]
        assert all(all(row[k] == prior[key][k] for k in non_context) for key,row in rows.items())
        prov = json.loads((OLD[c]/'provenance.json').read_text())
        prov.update(generated_utc=now(), metrics_sha256=digest(dest/'metrics.tsv'),
                    context_source_commit=p['source_commit'], context_source_version='0.2.0',
                    context_plan_sha256=plan_hash, context_audit_sha256=digest(output/'task_output_validation.json'),
                    retained_noncontext_source=str(OLD[c].relative_to(ROOT)),
                    retained_noncontext_metrics_sha256=baseline_hashes[c],
                    noncontext_columns_byte_identical=True,
                    script_sha256={**prov['script_sha256'],Path(__file__).name:digest(Path(__file__))})
        write_json(dest/'provenance.json', prov)
        production.summarize(dest, make_plots=False)
        # Use the existing scalar/domain validator, then record the overlay audit.
        subprocess.run([sys.executable, '-B', str(ROOT/'scripts/validate_e1_production.py'), str(dest)], check=True)
    write_json(output/'context_update_validation.json', dict(all_passed=True,
               changed_context_columns=dict(changed), noncontext_columns_byte_identical=True,
               source_metrics_sha256=baseline_hashes))


def main():
    global OVERLAY
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, default=runner.DEFAULT/'array_plan.json')
    ap.add_argument('--out', type=Path, default=ROOT/'results/local'/f'e1-context-{SUFFIX}')
    ap.add_argument('--scratch-dir', type=Path, default=Path('/tmp'))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--profile-db', type=Path,
                    help='Resume the downstream analyses with an already audited node-local profile DB')
    args = ap.parse_args()
    if args.profile_db:
        assert not str(args.profile_db.resolve()).startswith('/gs/bs/')
        audit = json.loads((args.out/'task_output_validation.json').read_text())
        assert audit['all_passed'] and audit['plan_sha256'] == digest(args.plan)
    else:
        args.out.mkdir(parents=True, exist_ok=False)
    with local_workspace(args.scratch_dir) as scratch:
        OVERLAY = args.profile_db or scratch/'profiles.sqlite'
        db = sqlite3.connect(OVERLAY)
        if not args.profile_db:
            audit_and_metrics(args.plan, args.out, db)
        assert db.execute('SELECT count(*) FROM profiles').fetchone()[0] == 131753
        db.close()
        subprocess.run([sys.executable, '-B', str(ROOT/'scripts/export_e1_context_reference_descriptors_0_2_0_20261023.py'),
                        '--plan', str(args.plan), '--db', str(OVERLAY), '--out', str(args.out)], check=True)
        cache_archive_reads()
        production.read_group = overlay_read_group
        short = ROOT/'results/runs'/runner.OLD['short']
        pkout = args.out/f'e1-pk-stratified-context-{SUFFIX}'
        import reanalyze_e1_pk_context_parallel
        reanalyze_e1_pk_context_parallel.run(short, OVERLAY, args.out/NAMES['short'], OLD['short'],
            ROOT/'results/local/e1-pk-stratified-89baea52-20261018', pkout, args.workers, ORIGINAL_READ_GROUP)
        localout = args.out/f'e1-local-errors-context-{SUFFIX}'
        local.analyze(argparse.Namespace(run=short, out=localout, workers=args.workers))
        for out in (pkout, localout):
            prov = json.loads((out/'provenance.json').read_text())
            prov.update(context_source_version='0.2.0',context_source_commit=runner.load(args.plan)['source_commit'],
                        context_plan_sha256=digest(args.plan),context_overlay_script_sha256=digest(Path(__file__)),
                        context_audit_sha256=digest(args.out/'task_output_validation.json'))
            write_json(out/'provenance.json',prov)
        for connection in CONNECTIONS.values():
            connection.close()
    print(f'completed: {args.out}', flush=True)


if __name__ == '__main__':
    main()
