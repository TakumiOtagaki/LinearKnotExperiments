#!/usr/bin/env python3
"""Check archive integrity and recompute the report without running kernels."""
import collections
import hashlib
import json
from pathlib import Path
import tarfile
from analyze import collect, pairs
from run import HERE, dump, fingerprint, sha


def main():
    manifest = json.loads((HERE / 'manifest.json').read_text())
    assert manifest['integrity_ok'] and manifest['before'] == manifest['after'] == fingerprint()
    for name, expected in manifest['artifacts_sha256'].items():
        assert sha(HERE / name) == expected, name
    runs = [json.loads(x) for x in (HERE / 'runs.jsonl').read_text().splitlines()]
    assert len(runs) == manifest['process_count']
    assert all(r['returncode'] == 0 and not r['timed_out'] for r in runs)
    members, draw_count, raw_counts = {}, 0, {}
    with tarfile.open(HERE / 'raw.tar.gz', 'r:gz') as archive:
        for member in archive:
            assert member.isfile() and member.name not in members
            f = archive.extractfile(member)
            h = hashlib.sha256(); size = 0
            if member.name.startswith('sampling/') and member.name.endswith('.stdout'):
                counts = collections.Counter()
                for i, line in enumerate(f):
                    h.update(line); size += len(line)
                    if i >= 2 and line.strip():
                        counts[frozenset(pairs(line.decode().strip()))] += 1
                        draw_count += 1
                raw_counts[Path(member.name).stem] = counts
            else:
                for chunk in iter(lambda: f.read(1048576), b''):
                    h.update(chunk); size += len(chunk)
            members[member.name] = dict(sha256=h.hexdigest(), bytes=size)
    assert members == manifest['members']
    assert set(members) == {r[k] for r in runs for k in ['stdout_member', 'stderr_member']} | {'build/probe', 'cli/input.txt'}
    report = json.loads((HERE / 'report.json').read_text())
    summary = json.loads((HERE / 'summary.json').read_text())
    rescored_report, rescored_summary = collect()
    assert report == rescored_report and summary == rescored_summary
    for row in report['sampling']:
        case = next(r for r in report['decoding'] if r['id'] == row['id'])
        counts = raw_counts[row['id'] + '_' + row['mode']]
        fixed = pairs(case['scaffold'])
        expected = {frozenset(fixed | pairs(gp)): v['observed'] for gp, v in row['frequencies'].items() if v['observed']}
        assert counts == expected and sum(counts.values()) == row['n']
    assert draw_count == summary['sampling']['draws']
    result = dict(ok=True, protected_fingerprints_unchanged=True, archive_members=len(members),
                  successful_processes=len(runs), recounted_draws=draw_count,
                  decoding_cases=len(report['decoding']), report_recomputed=True, summary_recomputed=True,
                  tools_sha256={n: sha(Path(__file__).resolve().parent / n) for n in ['probe.cc', 'run.py', 'analyze.py', 'validate.py']},
                  artifact_sha256={n: sha(HERE / n) for n in ['report.json', 'summary.json', 'raw.tar.gz', 'runs.jsonl', 'config.json', 'manifest.json']})
    dump(HERE / 'validation.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
