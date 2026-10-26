#!/usr/bin/env python3
"""Small, bounded structure-prediction pilot; all per-process files stay in /tmp.

References are complete coordinate sets, including noncanonical pairs. Crossing
subsets are computed independently in reference and prediction. No bracket-type
matching or one-base tolerance is used. This is a development pilot, not E1.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import random
import resource
import signal
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time
import zipfile

from prepare_biological_inputs import read_fasta, threshold_scaffold

ROOT = Path(__file__).resolve().parents[1]
BINS = [(50, 200), (200, 500), (500, 1000), (1000, 3001)]
SEED = 20261025
CANONICAL = {"AU", "UA", "CG", "GC", "GU", "UG"}
ENVIRONMENT = dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                   MKL_NUM_THREADS="1", LC_ALL="C", PYTHONDONTWRITEBYTECODE="1")
METHODS = ["LinearKnot-centroid", "LinearKnot-MFE", "LinearPartition-MEA",
           "ThreshKnot", "IPknot-LPV", "Spark", "Knotty"]


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def bracket_pairs(text):
    opening = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    closing = ")]}>abcdefghijklmnopqrstuvwxyz"
    stacks, result = {c: [] for c in opening}, set()
    for j, c in enumerate(text, 1):
        if c in opening:
            stacks[c].append(j)
        elif c in closing:
            stack = stacks[opening[closing.index(c)]]
            if not stack:
                raise ValueError("unmatched closing bracket")
            result.add((stack.pop(), j))
        elif c != ".":
            raise ValueError("unknown structure symbol: " + c)
    if any(stacks.values()):
        raise ValueError("unmatched opening bracket")
    return result


def validate_pairs(pairs, n):
    endpoints = []
    for i, j in pairs:
        if not 1 <= i < j <= n:
            raise ValueError("pair outside sequence")
        endpoints.extend((i, j))
    if len(endpoints) != len(set(endpoints)):
        raise ValueError("multiple partners for a base")


def crossing_pairs(pairs):
    ordered, result = sorted(pairs), set()
    for k, (i, j) in enumerate(ordered):
        for a, b in ordered[k + 1:]:
            if a >= j:
                break
            if i < a < j < b:
                result.update(((i, j), (a, b)))
    return result


def score(reference, prediction):
    tp = len(reference & prediction)
    nr, np = len(reference), len(prediction)
    return dict(tp=tp, reference=nr, predicted=np,
                sensitivity=tp / nr if nr else None,
                ppv=tp / np if np else None,
                f1=2 * tp / (nr + np) if nr + np else None)


def metrics(reference, prediction):
    cr, cp = crossing_pairs(reference), crossing_pairs(prediction)
    return dict(all=score(reference, prediction), crossing=score(cr, cp),
                reference_crossing_found_anywhere=len(cr & prediction))


def prepare(args):
    args.out.mkdir(parents=True, exist_ok=False)
    fasta = ROOT / "results/local/data-sources/bprna-local-clean-20261008/bpRNA_1m_90_filtered_clean.fa"
    sequences = dict(read_fasta(fasta))
    sources, candidates, excluded = [], {}, Counter()
    if args.bprna_zip:
        sources.append(dict(path=str(args.bprna_zip), sha256=sha(args.bprna_zip)))
        with zipfile.ZipFile(args.bprna_zip) as z:
            members = sorted(n for n in z.namelist() if n.endswith(".dbn"))
            for member in members:
                name = Path(member).stem
                if name not in sequences:
                    excluded["not_in_E1"] += 1
                    continue
                lines = [l.strip() for l in z.read(member).decode().splitlines()
                         if l.strip() and not l.startswith(("#", ">"))]
                if len(lines) != 2:
                    raise ValueError("expected sequence and DBN: " + member)
                seq, ref = lines
                seq = seq.upper().replace("T", "U")
                if seq != sequences[name]:
                    raise ValueError("E1 sequence mismatch: " + name)
                candidates[name] = dict(id=name, seq=seq, ref=ref,
                                        source_member=member, source="official bpRNA zip")
        scope = "E1 collection joined to complete official annotations"
    else:
        for filename in ("pk-recovery-bprna-20261016/pk_sel.json",
                         "pk-recovery-long-rrna-20261016/mid300_700_selection.json"):
            path = ROOT / "results/local" / filename
            sources.append(dict(path=str(path), sha256=sha(path)))
            for r in json.loads(path.read_text()):
                if sequences.get(r["id"]) != r["seq"]:
                    raise ValueError("E1 sequence mismatch")
                row = dict(id=r["id"], seq=r["seq"], ref=r["ref"], source=str(path.relative_to(ROOT)))
                if r["id"] in candidates and candidates[r["id"]]["ref"] != r["ref"]:
                    raise ValueError("inconsistent reference")
                candidates[r["id"]] = row
        scope = "fallback: saved 200+80 reference cohorts; NOT the full eligible E1 population"
    eligible = []
    for name, r in sorted(candidates.items()):
        if len(r["seq"]) != len(r["ref"]):
            raise ValueError("reference length mismatch")
        pairs = bracket_pairs(r["ref"])
        validate_pairs(pairs, len(r["seq"]))
        if not crossing_pairs(pairs):
            excluded["no_crossing_pairs"] += 1
            continue
        if not 50 <= len(r["seq"]) <= 3000:
            excluded["outside_length_range"] += 1
            continue
        r.update(length=len(r["seq"]), reference_pairs=sorted(pairs),
                 reference_crossing_pairs=sorted(crossing_pairs(pairs)),
                 noncanonical_reference_pairs=sum(r["seq"][i-1]+r["seq"][j-1] not in CANONICAL for i,j in pairs),
                 sequence_sha256=hashlib.sha256(r["seq"].encode()).hexdigest())
        eligible.append(r)
    rng, selected, bins = random.Random(SEED), [], []
    for lo, hi in BINS:
        pool = [r for r in eligible if lo <= r["length"] < hi]
        chosen = sorted(rng.sample(pool, min(args.per_bin, len(pool))), key=lambda r: (r["length"], r["id"]))
        label = f"{lo}-{hi if hi < 3001 else 3000}"
        for r in chosen:
            r["bin"] = label
        bins.append(dict(label=label, min_inclusive=lo, max_exclusive=hi,
                         available=len(pool), selected_ids=[r["id"] for r in chosen]))
        selected.extend(chosen)
    with (args.out / "cohort.jsonl").open("w") as f:
        for r in selected:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    dump(args.out / "selection.json", dict(created_utc=now(), seed=SEED,
        selection="one Python random.Random(seed); ID-sorted pools, sample without replacement per bin in ascending bin order",
        scope=scope, sources=sources, E1_fasta_sha256=sha(fasta), E1_records=len(sequences),
        requested_per_bin=args.per_bin, bins=bins, excluded=dict(excluded),
        eligible_ids=[r["id"] for r in eligible], reference_filter="at least one geometric crossing; no canonical-pair filter"))
    print(json.dumps(bins), flush=True)


def parse_prediction(stdout, seq, method, scaffold=None):
    n = len(seq)
    if method == "ThreshKnot":
        table = {}
        for line in stdout.splitlines():
            f = line.split()
            if len(f) == 3 and f[0].isdigit() and f[2].isdigit():
                i, j = int(f[0]), int(f[2])
                if i in table or not 1 <= i <= n or f[1] != seq[i-1]:
                    raise ValueError("invalid BPSEQ position/base")
                table[i] = j
        if set(table) != set(range(1, n + 1)):
            raise ValueError("incomplete BPSEQ")
        if any(j and table.get(j) != i for i,j in table.items()):
            raise ValueError("asymmetric BPSEQ")
        pairs = {(i,j) for i,j in table.items() if j > i}
    else:
        found = []
        for line in stdout.splitlines():
            f = line.split()
            if not f:
                continue
            s = f[0]
            if len(s) == n and s != seq and set(s) <= set(".()[]{}<>ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"):
                try:
                    pairs = bracket_pairs(s)
                except ValueError:
                    continue
                found.append((s, pairs))
        if method.startswith("LinearKnot"):
            # The CLI prints sequence, scaffold, and decoded structure.
            if len(found) != 2 or found[0][0] != scaffold:
                raise ValueError("unexpected LinearKnot output")
            pairs = found[1][1]
        else:
            if len(found) != 1:
                raise ValueError("expected exactly one predicted structure")
            pairs = found[0][1]
    validate_pairs(pairs, n)
    # All comparison programs should emit only canonical pairs.
    if any(seq[i-1]+seq[j-1] not in CANONICAL for i,j in pairs):
        raise ValueError("noncanonical predicted pair")
    if method == "LinearPartition-MEA" and crossing_pairs(pairs):
        raise ValueError("MEA baseline unexpectedly pseudoknotted")
    if method.startswith("LinearKnot") and not bracket_pairs(scaffold) <= pairs:
        raise ValueError("LinearKnot did not preserve scaffold")
    return pairs


def archive_bytes(tar, name, data):
    header = tarfile.TarInfo(name)
    header.size, header.mtime = len(data), int(time.time())
    tar.addfile(header, io.BytesIO(data))
    return hashlib.sha256(data).hexdigest()


def run(args):
    if args.scratch.resolve().is_relative_to(Path("/gs/bs")):
        raise ValueError("scratch must be node-local")
    args.scratch.mkdir(parents=True, exist_ok=True)
    records = [json.loads(l) for l in (args.out/"cohort.jsonl").read_text().splitlines()]
    if (args.out/"results.jsonl").exists():
        raise ValueError("refuse overwriting results; choose a new output")
    env = dict(os.environ, **ENVIRONMENT)
    env["CPARTY_FIXED_ENERGY_PARAM_FILE"] = str(args.parameter.resolve())
    jobs, order = [], []
    for r in records:
        jobs.append((r, "scaffold"))
        methods = METHODS[:]
        random.Random(SEED + len(order)).shuffle(methods)
        order.append(dict(id=r["id"], methods=methods))
        jobs.extend((r, method) for method in methods)
    setup = json.loads(args.setup.read_text()) if args.setup else {}
    binaries = {k: str(getattr(args,k).resolve()) if getattr(args,k) else None
                for k in ("linearknot", "linearpartition", "spark", "ipknot", "loader", "knotty", "parameter")}
    dump(args.out/"run.json", dict(status="running", started_utc=now(), purpose="small development pilot",
        hostname=os.uname().nodename, job_id=os.environ.get("JOB_ID"), nslots=os.environ.get("NSLOTS"),
        workers=1, warmups=0, repetitions=1, env=ENVIRONMENT, seed=SEED, schedule=order,
        invocation=sys.argv,
        timeouts_seconds=dict(default=120, Spark=30, Knotty=30), memory_limit_mib=4096,
        bin_timeout_stop=2, Spark_total_budget_seconds=180, total_budget_seconds=7200,
        settings=dict(LP_beam=200, scaffold_threshold=0.7, min_helix=2, LK_beam=100,
                      max_hairpin=100, max_case3_depth=100, gamma=1, LP_MEA_gamma=1,
                      ThreshKnot_threshold=0.3, IPknot_model="LinearPartition-V", IPknot_refinements=1),
        binaries=binaries, binary_sha256={k:sha(p) for k,p in binaries.items() if p}, setup=setup,
        script_sha256=sha(__file__), cohort_sha256=sha(args.out/"cohort.jsonl")))
    start = time.monotonic()
    scaffolds, scaffold_resources, timeout_counts = {}, {}, Counter()
    spark_seconds = 0.
    with tarfile.open(args.out/"raw_outputs.tar", "w") as tar, (args.out/"results.jsonl").open("w") as f:
        def save(row, folder=None):
            row["raw_sha256"] = {}
            if folder:
                for p in sorted(folder.iterdir()):
                    if p.is_file():
                        member = f"jobs/{row['id']}/{row['method']}/{p.name}"
                        row["raw_sha256"][member] = archive_bytes(tar, member, p.read_bytes())
            tar.fileobj.flush()
            os.fsync(tar.fileobj.fileno())
            row["archive_end"] = tar.offset
            f.write(json.dumps(row, ensure_ascii=False, allow_nan=False)+"\n")
            f.flush()
            os.fsync(f.fileno())
        for idx, (r, method) in enumerate(jobs):
            name, seq, label = r["id"], r["seq"], r["bin"]
            row = dict(id=name, method=method, bin=label, length=len(seq), started_utc=now(),
                       elapsed_seconds=None, time_wall_seconds=None, peak_rss_kib=None,
                       status=None, metrics=None, prediction_pairs=None)
            reason = None
            if time.monotonic()-start > 7200:
                reason = "total pilot wall budget exhausted"
            elif method == "Knotty" and (not args.knotty or len(seq)>200):
                reason = "unavailable binary" if not args.knotty else "length >200 nt"
            elif timeout_counts[method,label] >= 2:
                reason = "two earlier timeouts in this method/bin"
            elif method == "Spark" and spark_seconds >= 180:
                reason = "Spark total 180-second budget exhausted"
            elif method in ("LinearKnot-centroid", "LinearKnot-MFE", "Spark") and name not in scaffolds:
                reason = "scaffold failed"
            elif method == "IPknot-LPV" and not args.ipknot:
                reason = "unavailable binary"
            if reason:
                row.update(status="skipped", error=reason, finished_utc=now())
                save(row)
                continue
            with tempfile.TemporaryDirectory(prefix="job-", dir=args.scratch) as tmp:
                folder = Path(tmp)
                input_file, bpp_file = folder/"input.fa", folder/"bpp.tsv"
                input_file.write_text(f">{name}\n{seq}\n")
                scaffold = scaffolds.get(name)
                if method == "scaffold":
                    argv = [str(args.linearpartition), "-V", "-b", "200", "--dangles", "2", "-c", "0.7", "-r", str(bpp_file)]
                elif method.startswith("LinearKnot"):
                    input_file.write_text(f">{name}\n{seq}\n{scaffold}\n")
                    subcommand = "centroid" if method.endswith("centroid") else "fold"
                    argv = [str(args.linearknot), subcommand, "-i", str(input_file), "-b", "100", "--max-hairpin", "100", "--max-case3-depth", "100", "--ranking", "outer"]
                    if subcommand == "centroid":
                        argv += ["--gamma", "1"]
                elif method in ("ThreshKnot", "LinearPartition-MEA"):
                    argv = [str(args.linearpartition), "-V", "-b", "200", "--dangles", "2"]
                    argv += ["-T", "--threshold", "0.3"] if method == "ThreshKnot" else ["-M", "-g", "1"]
                elif method == "IPknot-LPV":
                    argv = [str(args.ipknot), "-e", "LinearPartition-V", "--beam-size", "200", "-n", "1", str(input_file)]
                    if args.loader:
                        argv = [str(args.loader), "--library-path", args.library_path, *argv]
                elif method == "Spark":
                    argv = [str(args.spark), "-P", str(args.parameter), "-d2", "-r", scaffold, seq]
                else:
                    argv = [str(args.knotty), seq, "-w"]
                timeout = 30 if method in ("Spark", "Knotty") else 120
                row.update(argv=argv, cwd=str(folder), timeout_seconds=timeout, memory_limit_mib=4096)
                dump(folder/"command.json", row)
                (folder/"stdin.txt").write_text(seq+"\n")
                timed = ["/usr/bin/time", "-f", "%e\t%M\t%U\t%S", "-o", str(folder/"time.tsv"), "--", *argv]
                def limits():
                    cap = 4096*1024*1024
                    resource.setrlimit(resource.RLIMIT_AS, (cap,cap))
                    resource.setrlimit(resource.RLIMIT_CORE, (0,0))
                then = time.monotonic()
                proc = subprocess.Popen(timed, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, env=env, cwd=folder,
                                        start_new_session=True, preexec_fn=limits)
                try:
                    stdout, stderr = proc.communicate((seq+"\n").encode(), timeout=timeout)
                    row["status"] = "success" if proc.returncode == 0 else "failed"
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    stdout, stderr = proc.communicate()
                    row.update(status="timeout", error="wall timeout")
                    timeout_counts[method,label] += 1
                row.update(elapsed_seconds=time.monotonic()-then, returncode=proc.returncode, finished_utc=now())
                if method == "Spark":
                    spark_seconds += row["elapsed_seconds"]
                (folder/"stdout.txt").write_bytes(stdout)
                (folder/"stderr.txt").write_bytes(stderr)
                tpath = folder/"time.tsv"
                if tpath.exists():
                    fields = tpath.read_text().strip().splitlines()[-1].split("\t")
                    if len(fields)==4:
                        row.update(time_wall_seconds=float(fields[0]), peak_rss_kib=int(fields[1]),
                                   user_seconds=float(fields[2]), system_seconds=float(fields[3]))
                if row["status"] == "success":
                    try:
                        if method == "scaffold":
                            scaffold = threshold_scaffold(seq, bpp_file, 0.7, 2)
                            scaffolds[name] = scaffold
                            scaffold_resources[name] = {k: row[k] for k in ("elapsed_seconds", "peak_rss_kib")}
                            row["scaffold"] = scaffold
                            pairs = bracket_pairs(scaffold)
                            row["scaffold_pairs"] = sorted(pairs)
                            row["metrics"] = metrics({tuple(p) for p in r["reference_pairs"]}, pairs)
                            (folder/"scaffold.txt").write_text(scaffold+"\n")
                        else:
                            pairs = parse_prediction(stdout.decode(), seq, method, scaffold)
                            row["prediction_pairs"] = sorted(pairs)
                            row["metrics"] = metrics({tuple(p) for p in r["reference_pairs"]}, pairs)
                            if method in ("LinearKnot-centroid", "LinearKnot-MFE", "Spark"):
                                row["scaffold_resources"] = scaffold_resources[name]
                                row["end_to_end_seconds"] = row["elapsed_seconds"] + scaffold_resources[name]["elapsed_seconds"]
                                row["end_to_end_peak_rss_kib"] = max(row["peak_rss_kib"],scaffold_resources[name]["peak_rss_kib"])
                            else:
                                row["end_to_end_seconds"] = row["elapsed_seconds"]
                                row["end_to_end_peak_rss_kib"] = row["peak_rss_kib"]
                            if method == "Spark":
                                row["missing_scaffold_pairs"] = sorted(bracket_pairs(scaffold)-pairs)
                                row["accuracy_reference_only"] = True
                    except (ValueError, IndexError) as e:
                        row.update(status="invalid_output", error=str(e))
                dump(folder/"result.json", row)
                save(row, folder)
                print(f"{idx+1}/{len(jobs)} {name} n={len(seq)} {method} {row['status']} {row['elapsed_seconds']:.3f}s", flush=True)
        for path in [Path(__file__), args.out/"cohort.jsonl", args.out/"selection.json", args.out/"run.json"]:
            archive_bytes(tar,"metadata/"+path.name,path.read_bytes())
    manifest = json.loads((args.out/"run.json").read_text())
    manifest.update(status="complete", finished_utc=now(), total_elapsed_seconds=time.monotonic()-start,
                    raw_outputs_sha256=sha(args.out/"raw_outputs.tar"),
                    results_sha256=sha(args.out/"results.jsonl"),
                    final_binary_sha256={k:sha(p) for k,p in binaries.items() if p})
    if manifest["binary_sha256"] != manifest["final_binary_sha256"]:
        raise ValueError("binary changed during run")
    dump(args.out/"run.json",manifest)


def summarize(args):
    cohort = [json.loads(l) for l in (args.out/"cohort.jsonl").read_text().splitlines()]
    rows = [json.loads(l) for l in (args.out/"results.jsonl").read_text().splitlines()]
    summary = []
    for lo, hi in BINS:
        label = f"{lo}-{hi if hi<3001 else 3000}"
        for method in ["scaffold", *METHODS]:
            subset = [r for r in rows if r["bin"]==label and r["method"]==method]
            ok = [r for r in subset if r["status"]=="success"]
            result = dict(bin=label,method=method,selected=sum(r["bin"]==label for r in cohort),
                          status_counts=dict(Counter(r["status"] for r in subset)), success=len(ok), metrics={})
            for family in ("all","crossing"):
                tp = sum(r["metrics"][family]["tp"] for r in ok)
                nr = sum(r["metrics"][family]["reference"] for r in ok)
                np = sum(r["metrics"][family]["predicted"] for r in ok)
                result["metrics"][family] = dict(tp=tp,reference=nr,predicted=np,
                    sensitivity=tp/nr if nr else None, ppv=tp/np if np else None,
                    f1=2*tp/(nr+np) if nr+np else None,
                    macro_f1=statistics.mean(r["metrics"][family]["f1"] for r in ok) if ok else None)
            for field in ("elapsed_seconds","peak_rss_kib","end_to_end_seconds","end_to_end_peak_rss_kib"):
                values=[r[field] for r in ok if r.get(field) is not None]
                result["median_"+field]=statistics.median(values) if values else None
            result["missing_scaffold_RNAs"] = sum(bool(r.get("missing_scaffold_pairs")) for r in ok)
            summary.append(result)
    # Also retain a common-success comparison rather than silently pairing
    # Spark's easier survivors with the full cohort's results.
    paired = []
    for label in {r["bin"] for r in cohort}:
        methods = [m for m in METHODS if m != "Knotty"]
        common = set.intersection(*[{r["id"] for r in rows if r["bin"]==label and r["method"]==m and r["status"]=="success"} for m in methods])
        paired.append(dict(bin=label, methods=methods, ids=sorted(common), count=len(common)))
    dump(args.out/"summary.json", dict(created_utc=now(),groups=summary,common_success=paired,
        aggregation="micro sums across successful RNAs; failed and skipped denominators shown separately; macro F1 also retained",
        empty_prediction="PPV null, sensitivity=0, F1=0 for nonempty references",
        crossing_definition="pairs crossing at least one other pair in their own structure; TP is intersection of the two crossing subsets"))
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def validate(args):
    """Recompute every metric from archived outputs and check durable hashes."""
    cohort = {r['id']:r for r in map(json.loads,(args.out/'cohort.jsonl').read_text().splitlines())}
    rows = list(map(json.loads,(args.out/'results.jsonl').read_text().splitlines()))
    manifest = json.loads((args.out/'run.json').read_text())
    assert manifest['status']=='complete'
    assert manifest['raw_outputs_sha256']==sha(args.out/'raw_outputs.tar')
    assert manifest['results_sha256']==sha(args.out/'results.jsonl')
    assert manifest['cohort_sha256']==sha(args.out/'cohort.jsonl')
    assert manifest['binary_sha256']==manifest['final_binary_sha256']
    assert len({(r['id'],r['method']) for r in rows})==len(rows)==len(cohort)*8
    def slow_cross(pairs):
        return {p for p in pairs if any(p[0]<q[0]<p[1]<q[1] or q[0]<p[0]<q[1]<p[1] for q in pairs)}
    checked_members, parsed, errors = 0, 0, []
    with tarfile.open(args.out/'raw_outputs.tar') as tar:
        members={m.name:m for m in tar}
        assert len(members)==len(tar.getmembers()), 'duplicate tar member'
        archived_script=tar.extractfile('metadata/'+Path(__file__).name).read()
        assert hashlib.sha256(archived_script).hexdigest()==manifest['script_sha256']
        for row in rows:
            for name, digest in row['raw_sha256'].items():
                assert hashlib.sha256(tar.extractfile(members[name]).read()).hexdigest()==digest
                checked_members += 1
            if row['status']!='success':
                continue
            r=cohort[row['id']]
            reference={tuple(p) for p in r['reference_pairs']}
            assert slow_cross(reference)=={tuple(p) for p in r['reference_crossing_pairs']}
            prefix=f"jobs/{row['id']}/{row['method']}/"
            if row['method']=='scaffold':
                prediction=bracket_pairs(row['scaffold'])
                assert not slow_cross(prediction)
                posterior={}
                for line in tar.extractfile(prefix+'bpp.tsv').read().decode().splitlines():
                    f=line.split()
                    if len(f)==3 and f[0].isdigit():
                        posterior[int(f[0]),int(f[1])]=float(f[2])
                above={p for p,v in posterior.items() if v>=0.7}
                independent={p for p in above if (p[0]-1,p[1]+1) in above or (p[0]+1,p[1]-1) in above}
                assert prediction==independent
            else:
                stdout=tar.extractfile(prefix+'stdout.txt').read().decode()
                scaffold=None
                if row['method'].startswith('LinearKnot'):
                    scaffold=tar.extractfile(f"jobs/{row['id']}/scaffold/scaffold.txt").read().decode().strip()
                prediction=parse_prediction(stdout,r['seq'],row['method'],scaffold)
                assert prediction=={tuple(p) for p in row['prediction_pairs']}
                parsed += 1
            for family,ref,pred in [('all',reference,prediction),('crossing',slow_cross(reference),slow_cross(prediction))]:
                got=row['metrics'][family]
                assert got['tp']==sum(p in ref for p in pred)
                assert got['reference']==len(ref) and got['predicted']==len(pred)
                expected=2*len(ref&pred)/(len(ref)+len(pred)) if ref or pred else None
                assert got['f1']==expected
            assert row['peak_rss_kib']>0 and row['elapsed_seconds']>=0
    result=dict(status='passed',created_utc=now(),RNAs=len(cohort),rows=len(rows),
                raw_members_verified=checked_members,predictions_reparsed=parsed,
                scaffold_rule_independently_checked=True,metrics_independently_recounted=True,
                binary_hashes_unchanged=True,errors=errors)
    dump(args.out/'validation.json',result)
    print(json.dumps(result))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest="command",required=True)
    p=subs.add_parser("prepare");p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bprna-zip",type=Path);p.add_argument("--per-bin",type=int,default=20)
    p=subs.add_parser("run");p.add_argument("--out",type=Path,required=True)
    p.add_argument("--scratch",type=Path,required=True);p.add_argument("--setup",type=Path)
    for name in ("linearknot","linearpartition","spark","parameter"):
        p.add_argument("--"+name,type=Path,required=True)
    for name in ("ipknot","loader","knotty"):
        p.add_argument("--"+name,type=Path)
    p.add_argument("--library-path")
    for name in ('summarize','validate'):
        p=subs.add_parser(name);p.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    globals()[args.command](args)


if __name__=="__main__":
    main()
