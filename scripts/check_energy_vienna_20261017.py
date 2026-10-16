#!/usr/bin/env python3
"""Independent, small empty-scaffold ViennaRNA comparison; no source edits.

Builds a read-only observer in node-local scratch, plus a separately labelled
diagnostic observer changing ONLY the two multiloop-closing neighbor guards.
The latter is not an adopted implementation or a manuscript measurement.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import resource
import subprocess
import tarfile
import tempfile
import time

import RNA


ROOT = Path(__file__).resolve().parents[1]
HEADER = Path("scfg2/replay/w_final_exact_inside.hh")
OLD = ("si = j < owner_.length() ? owner_.encoded_base(j-1) : -1;\n"
       "        sj = i > 1 ? owner_.encoded_base(i+1) : -1;")
NEW = ("si = owner_.encoded_base(j-1);\n"
       "        sj = owner_.encoded_base(i+1);")


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def pairs(db):
    stack, result = [], []
    for i, x in enumerate(db, 1):
        if x == "(":
            stack.append(i)
        elif x == ")":
            result.append((stack.pop(), i))
        elif x != ".":
            raise ValueError("expected pseudoknot-free dot-bracket")
    if stack:
        raise ValueError("unbalanced structure")
    return sorted(result)


def db_from_pairs(n, values):
    db = ["."] * n
    for i, j in values:
        db[i-1], db[j-1] = "(", ")"
    return "".join(db)


def loop_tree(db):
    p = pairs(db)
    children = {i: [] for i, _ in p}
    children[0] = []
    stack = [0]
    ends = dict(p)
    for i, j in p:
        while stack[-1] and i > ends[stack[-1]]:
            stack.pop()
        children[stack[-1]].append(i)
        stack.append(i)
    return p, children


def boundary_multiloops(db):
    p, children = loop_tree(db)
    return [(i, j) for i, j in p if len(children[i]) >= 2 and (i == 1 or j == len(db))]


def make_md(dangles=2, smooth=1):
    md = RNA.md()
    md.temperature, md.betaScale, md.dangles = 37., 1., dangles
    md.pf_smooth, md.min_loop_size = smooth, 3
    md.noLP = md.noGU = md.noGUclosure = md.logML = md.circ = md.gquad = 0
    md.special_hp, md.compute_bpp, md.uniq_ML = 1, 1, 1
    # ViennaRNA's reference salt concentration gives no salt correction.
    md.salt = 1.021
    return md


def pf_logz(fc, n):
    rounded_energy = fc.pf()[1]
    q = fc.exp_matrices.q1k[n]
    if q <= 0 or not math.isfinite(q):
        raise ValueError("ViennaRNA PF is unreachable or nonfinite")
    # pf() returns a C float in this binding. Read the double precision root.
    logz = math.log(q) + n * math.log(fc.exp_params.pf_scale)
    return logz, -fc.exp_params.kT/1000 * logz, rounded_energy


def vienna_ensemble(sequence, dangles=2):
    fc = RNA.fold_compound(sequence, make_md(dangles))
    mfe_db, mfe_rounded = fc.mfe()
    logz, energy, rounded = pf_logz(fc, len(sequence))
    return fc, dict(vienna_logz=logz, vienna_pf_energy=energy,
                    vienna_pf_energy_float=rounded, vienna_rt=fc.exp_params.kT/1000,
                    vienna_mfe_structure=mfe_db,
                    vienna_mfe_energy=round(float(mfe_rounded), 2))


def structure_reference(sequence, db, dangles):
    """Independent exact single-structure PF score and integer eval_structure."""
    fc = RNA.fold_compound(sequence, make_md(dangles))
    integer_energy = fc.eval_structure(db)
    # Every dot is forced unpaired and every pair is enforced, leaving one R.
    constraint = db.replace(".", "x")
    if not fc.hc_add_from_db(constraint,
                            RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP):
        raise ValueError("ViennaRNA rejected complete structure constraint")
    logz, pf_energy, rounded = pf_logz(fc, len(sequence))
    # Recheck that each forced pair has probability one. This also guards the
    # binding-specific hard-constraint options rather than assuming their effect.
    bpp = fc.bpp()
    forced_error = max([abs(bpp[i][j]-1.) for i, j in pairs(db)] + [0.])
    if forced_error > 1e-9:
        raise ValueError("hard constraint did not enforce every target pair")
    pt = RNA.ptable(db)
    _, children = loop_tree(db)
    breakdown = []
    for i in [0] + sorted(children.keys()-{0}):
        if i == 0:
            kind = "exterior"
        elif not children[i]:
            kind = "hairpin"
        elif len(children[i]) == 1:
            kind = "stack_or_internal"
        else:
            kind = "multiloop"
        breakdown.append(dict(i=i, kind=kind, energy=fc.eval_loop_pt(i, pt)/100.))
    if abs(sum(x["energy"] for x in breakdown)-integer_energy) > 2e-5:
        raise ValueError("ViennaRNA loop energies do not sum to eval_structure")
    return dict(vienna_eval_energy=float(integer_energy),
                vienna_structure_pf_energy=pf_energy, vienna_structure_logz=logz,
                vienna_structure_pf_float=rounded,
                vienna_constraint_max_pair_error=forced_error,
                vienna_loop_energies=breakdown)


def archive_bytes(archive, name, value):
    data = value.encode() if isinstance(value, str) else value
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mtime = 0
    archive.addfile(info, io.BytesIO(data))


def child_limits():
    # Keep this explicitly small: one serial observer, 2 GiB, 40 CPU seconds.
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (40, 41))


def invoke(binary, param, record, mode, dangles, targets, archive, name, trace=False):
    argv = [str(binary), str(param), record["sequence"], record.get("scaffold", "."*record["length"]),
            str(dangles), mode, str(int(trace))]
    stdin = "".join(x+"\n" for x in targets)
    start = time.monotonic()
    try:
        run = subprocess.run(argv, input=stdin, text=True, capture_output=True,
                             timeout=45, preexec_fn=child_limits)
        stdout, stderr, code = run.stdout, run.stderr, run.returncode
        status = "success" if code == 0 else "failed"
    except subprocess.TimeoutExpired as e:
        stdout = e.stdout or b""
        stderr = e.stderr or b""
        stdout = stdout.decode() if isinstance(stdout, bytes) else stdout
        stderr = stderr.decode() if isinstance(stderr, bytes) else stderr
        code, status = None, "timeout"
    meta = dict(command=argv, status=status, returncode=code, elapsed_seconds=time.monotonic()-start,
                stdin=stdin, stdout_sha256=hashlib.sha256(stdout.encode()).hexdigest(),
                stderr_sha256=hashlib.sha256(stderr.encode()).hexdigest())
    archive_bytes(archive, name+"/command.json", json.dumps(meta, indent=2)+"\n")
    archive_bytes(archive, name+"/stdout.tsv", stdout)
    archive_bytes(archive, name+"/stderr.txt", stderr)
    if code != 0:
        return dict(**meta, error=stderr[-1500:])
    out = dict(**meta, bpp={}, targets=[], rules=[])
    for line in stdout.splitlines():
        fields = line.split("\t")
        tag, f = fields[0], fields[1:]
        if tag == "META":
            out["rt"], out["temperature"] = map(float, f[:2])
            out["dangles"], out["pf_smooth"], out["min_hairpin"], out["max_internal"] = map(int, f[2:])
        elif tag == "MLPARAM":
            out["ml_parameters"] = list(map(float, f))
        elif tag == "Z":
            out["logz"], out["pf_energy"] = map(float, f)
        elif tag == "MFE":
            out["mfe_energy"], out["mfe_structure"] = float(f[0]), f[1]
        elif tag == "BPP":
            out["bpp"][int(f[0]), int(f[1])] = float(f[2])
        elif tag == "EVAL":
            out["targets"].append(dict(index=int(f[0]), same_pairs=bool(int(f[1])),
                                       reachable=bool(int(f[2])), energy=float(f[3])))
        elif tag == "RULE":
            out["rules"].append(dict(index=int(f[0]), rule=f[1], i=int(f[2]), j=int(f[3]),
                                     ip=int(f[4]), jp=int(f[5]), k=int(f[6]), l=int(f[7]), energy=float(f[8])))
    return out


def constructed_inputs():
    result = []
    ml = "(.(((....))).(((....))).)"
    adjacent = "((((....)))(((....))))"
    for closing in ["GC", "AU", "GU"]:
        for name, core in [("gapped", ml), ("adjacent", adjacent)]:
            for left, right in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                db = "."*left+core+"."*right
                sequence = ["A"]*len(db)
                for i, j in pairs(db):
                    sequence[i-1], sequence[j-1] = "G", "C"
                sequence[left], sequence[len(db)-right-1] = closing
                result.append(dict(id=f"ml_{name}_{closing}_pad{left}{right}", group="constructed",
                                   sequence="".join(sequence), length=len(db), target=db))
    for stem in ["GC", "AU", "GU"]:
        for left, right in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            db = "."*left+"(...)"+"."*right
            sequence = "A"*left+stem[0]+"AAA"+stem[1]+"A"*right
            result.append(dict(id=f"hairpin_{stem}_pad{left}{right}", group="constructed",
                               sequence=sequence, length=len(db), target=db))
    # Exterior branches directly abut each other. dangles=2 uses paired neighbors.
    result.append(dict(id="exterior_adjacent_stems", group="constructed", sequence="GGGAAACCCGGGAAACCC",
                       length=18, target="(((...)))(((...)))"))
    return result


def select_inputs(fasta, seed, random_count, bprna_count):
    rng = random.Random(seed)
    result = []
    for index in range(random_count):
        n = [20, 40, 80, 120, 200][index % 5]
        s = "".join(rng.choice("ACGU") for _ in range(n))
        result.append(dict(id=f"random_{index:03d}_n{n}", group="random", sequence=s, length=n))
    records = []
    identifier, chunks = None, []
    for line in fasta.read_text().splitlines()+[">"]:
        if line.startswith(">"):
            if identifier is not None:
                s = "".join(chunks).upper().replace("T", "U")
                if 50 <= len(s) <= 200 and set(s) <= set("ACGU"):
                    records.append(dict(id=identifier, group="bpRNA", sequence=s, length=len(s)))
            identifier, chunks = line[1:].split()[0] if line[1:] else None, []
        else:
            chunks.append(line.strip())
    for k, (lo, hi) in enumerate([(50, 99), (100, 149), (150, 200)]):
        cohort = sorted((r for r in records if lo <= r["length"] <= hi), key=lambda r: r["id"])
        count = bprna_count // 3 + (k < bprna_count % 3)
        result.extend(rng.sample(cohort, count))
    for r in result:
        r["sequence_sha256"] = hashlib.sha256(r["sequence"].encode()).hexdigest()
    return result


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w") as f:
        w = csv.DictWriter(f, columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def summarize(rows, structures):
    summary = []
    for group in sorted({r["group"] for r in rows}):
        for variant in ["pinned", "inner_neighbors_diagnostic"]:
            selected = [r for r in rows if r["group"] == group and r["variant"] == variant]
            success = [r for r in selected if r["status"] == "success"]
            if not success:
                continue
            summary.append(dict(group=group, variant=variant, inputs=len(selected), successes=len(success),
                                min_length=min(r["length"] for r in success),
                                max_length=max(r["length"] for r in success),
                                max_abs_logz_difference=max(abs(r["delta_logz"]) for r in success),
                                max_abs_pf_energy_difference=max(abs(r["delta_pf_energy"]) for r in success),
                                max_abs_bpp_difference=max(r["max_abs_bpp_difference"] for r in success),
                                mfe_same_pair_set=sum(r["mfe_same_pair_set"] for r in success),
                                max_abs_mfe_vs_integer_difference=max(abs(r["mfe_delta_vs_vienna_integer"]) for r in success),
                                logz_difference_above_1e_9=sum(abs(r["delta_logz"]) > 1e-9 for r in success),
                                bpp_difference_above_1e_9=sum(r["max_abs_bpp_difference"] > 1e-9 for r in success)))
    return dict(ensemble=summary, structure_evaluations=len(structures),
                structure_unreachable=sum(not r["reachable"] for r in structures),
                boundary_multiloop_targets=sum(bool(json.loads(r["boundary_multiloops"])) for r in structures),
                max_abs_pinned_structure_pf_difference=max(abs(r["delta_pinned_vs_vienna_pf"]) for r in structures if r["reachable"]),
                max_abs_diagnostic_structure_pf_difference=max(abs(r["delta_diagnostic_vs_vienna_pf"]) for r in structures if r["reachable"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT/"deps/LinearCParty-d02d261c")
    parser.add_argument("--build", type=Path, default=ROOT/"build/linear-cparty-d02d261c")
    parser.add_argument("--fasta", type=Path, default=ROOT/"results/local/data-sources/bprna-local-clean-20261008/bpRNA_1m_90_filtered_clean.fa")
    parser.add_argument("--output", type=Path, default=ROOT/"results/local/energy-vienna-check-20261017")
    parser.add_argument("--random-count", type=int, default=150)
    parser.add_argument("--bprna-count", type=int, default=150)
    parser.add_argument("--seed", type=int, default=20261017)
    parser.add_argument("--sample-inputs", type=int, default=12)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("output must be new (preserve previous observations)")
    if args.random_count > 200 or args.bprna_count > 200:
        raise SystemExit("this observer is intentionally limited to a small local comparison")
    args.output.mkdir(parents=True)
    source, build = args.source.resolve(), args.build.resolve()
    param = source/"params/rna_Turner04.par"
    probe_source = Path(__file__).with_name("energy_vienna_probe_20261017.cc")
    dependencies = [Path(__file__), probe_source, source/"src"/HEADER,
                    source/"src/scfg2/weights/log_loops.hh", source/"src/ViennaRNA/params/params.c",
                    build/"libcparty_scfg2_exact.a", build/"libRNA.a", param, args.fasta]
    hashes = {str(p): digest(p) for p in dependencies}
    records = select_inputs(args.fasta, args.seed, args.random_count, args.bprna_count)+constructed_inputs()
    for r in records:
        r.setdefault("sequence_sha256", hashlib.sha256(r["sequence"].encode()).hexdigest())
    (args.output/"inputs.jsonl").write_text("".join(json.dumps(r)+"\n" for r in records))
    provenance = dict(started_utc=now(), purpose="development independent energy audit, not performance measurement",
                      seed=args.seed, source=str(source), build=str(build),
                      source_commit=subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
                      vienna_version=RNA.__version__, vienna_module=RNA.__file__,
                      command=os.sys.argv, dependencies_sha256_before=hashes,
                      model_details=str(make_md()), limits=dict(threads=1, child_address_space_gib=2,
                      child_cpu_seconds=40, child_wall_seconds=45),
                      diagnostic_header_replacement=dict(path=str(HEADER), old=OLD, new=NEW),
                      warnings=["MFE in LinearCParty optimizes PF-smoothed weights; Vienna mfe/eval_structure use integer parameters.",
                                "Nonempty-scaffold partition functions are not compared: no public no-crossing-extension switch was identified."])
    (args.output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    RNA.params_load(str(param))
    RNA.init_rand(args.seed)
    ensemble_rows, structure_rows = [], []
    with tempfile.TemporaryDirectory(prefix="energy-vienna-", dir=os.environ.get("TMPDIR", "/tmp")) as scratch_name, \
            tarfile.open(args.output/"raw_outputs.tar.gz", "w:gz") as archive:
        scratch = Path(scratch_name)
        if str(scratch.resolve()).startswith("/gs/bs/"):
            raise SystemExit("scratch must be node-local")
        header_text = (source/"src"/HEADER).read_text()
        if header_text.count(OLD) != 1:
            raise ValueError("pinned closing-guard witness no longer matches")
        shadow = scratch/"shadow"/HEADER
        shadow.parent.mkdir(parents=True)
        shadow.write_text(header_text.replace(OLD, NEW))
        provenance["diagnostic_header_sha256"] = digest(shadow)
        binaries = {}
        for variant in ["pinned", "inner_neighbors_diagnostic"]:
            binary = scratch/variant
            cmd = ["c++", "-std=c++17", "-O2"]
            if variant != "pinned":
                cmd += ["-I"+str(scratch/"shadow")]
            cmd += ["-I"+str(source/"src"), "-I"+str(build/"generated"), str(probe_source),
                    str(build/"libcparty_scfg2_exact.a"), str(build/"libRNA.a"), "-lm", "-o", str(binary)]
            run = subprocess.run(cmd, text=True, capture_output=True, timeout=120)
            archive_bytes(archive, variant+"/build.json", json.dumps(dict(command=cmd, returncode=run.returncode,
                          stdout=run.stdout, stderr=run.stderr), indent=2))
            run.check_returncode()
            binaries[variant] = binary
        provenance["observer_binary_sha256"] = {k: digest(p) for k, p in binaries.items()}
        for index, record in enumerate(records):
            rid, sequence, n = record["id"], record["sequence"], record["length"]
            fc, reference = vienna_ensemble(sequence)
            vbpp = fc.bpp()
            archive_bytes(archive, rid+"/vienna_reference.json", json.dumps(dict(
                sequence=sequence, model_details=str(make_md()), **reference)))
            archive_bytes(archive, rid+"/vienna_bpp.tsv", "i\tj\tprobability\n"+"".join(
                f"{i}\t{j}\t{vbpp[i][j]:.17g}\n" for i in range(1, n+1)
                for j in range(i+1, n+1) if vbpp[i][j] > 0.))
            observed = {}
            for variant, binary in binaries.items():
                run = invoke(binary, param, record, "ensemble", 2, [], archive, rid+"/"+variant)
                observed[variant] = run
                row = dict(id=rid, group=record["group"], length=n, variant=variant, status=run["status"],
                           elapsed_seconds=run["elapsed_seconds"], **reference)
                if run["status"] == "success":
                    differences = [(abs(run["bpp"].get((i, j), 0.)-vbpp[i][j]), i, j)
                                   for i in range(1, n+1) for j in range(i+1, n+1)]
                    err, i, j = max(differences)
                    row.update(logz=run["logz"], pf_energy=run["pf_energy"], delta_logz=run["logz"]-reference["vienna_logz"],
                               delta_pf_energy=run["pf_energy"]-reference["vienna_pf_energy"],
                               max_abs_bpp_difference=err, max_bpp_pair=f"{i},{j}",
                               max_bpp_lcp=run["bpp"].get((i, j), 0.), max_bpp_vienna=vbpp[i][j],
                               mfe_energy=run["mfe_energy"], mfe_structure=run["mfe_structure"],
                               mfe_same_pair_set=pairs(run["mfe_structure"]) == pairs(reference["vienna_mfe_structure"]),
                               mfe_delta_vs_vienna_integer=run["mfe_energy"]-reference["vienna_mfe_energy"])
                    if abs(run["rt"]-reference["vienna_rt"]) > 1e-12:
                        raise ValueError("thermal energy differs")
                else:
                    row["error"] = run.get("error")
                ensemble_rows.append(row)
            # Every sequence gets its two independently predicted MFE targets.
            # Selected inputs also get ViennaRNA Boltzmann samples. Constructed
            # cases get detailed traces, dangles=0 controls and G partitions.
            targets = {reference["vienna_mfe_structure"]: {"vienna_mfe"}}
            for variant, run in observed.items():
                if run["status"] == "success":
                    targets.setdefault(run["mfe_structure"], set()).add(variant+"_mfe")
            if index < args.sample_inputs:
                for db in fc.pbacktrack(3):
                    targets.setdefault(db, set()).add("vienna_boltzmann_sample")
            if "target" in record:
                targets.setdefault(record["target"], set()).add("constructed_target")
            target_list = list(targets)
            for dangles in ([2, 0] if "target" in record else [2]):
                evals = {variant: invoke(binary, param, record, "eval", dangles, target_list, archive,
                                         rid+f"/eval_d{dangles}_"+variant, trace="target" in record)
                         for variant, binary in binaries.items()}
                if any(r["status"] != "success" for r in evals.values()):
                    raise RuntimeError("structure evaluator failed; raw failure is preserved")
                for ti, db in enumerate(target_list):
                    ref = structure_reference(sequence, db, dangles)
                    pinned = evals["pinned"]["targets"][ti]
                    diagnostic = evals["inner_neighbors_diagnostic"]["targets"][ti]
                    sr = dict(id=rid, group=record["group"], length=n, structure=db, scaffold="."*n,
                              roles=",".join(sorted(targets[db])), dangles=dangles,
                              boundary_multiloops=json.dumps(boundary_multiloops(db)),
                              reachable=pinned["same_pairs"] and pinned["reachable"],
                              pinned_energy=pinned["energy"], diagnostic_energy=diagnostic["energy"],
                              delta_pinned_vs_vienna_pf=pinned["energy"]-ref["vienna_structure_pf_energy"],
                              delta_diagnostic_vs_vienna_pf=diagnostic["energy"]-ref["vienna_structure_pf_energy"],
                              delta_pinned_vs_vienna_eval=pinned["energy"]-ref["vienna_eval_energy"],
                              **{k: v for k, v in ref.items() if k != "vienna_loop_energies"})
                    structure_rows.append(sr)
                    if "target" in record:
                        archive_bytes(archive, rid+f"/target{ti}_d{dangles}_vienna.json", json.dumps(ref, indent=2))
            # Check identical pseudoknot-free completed structures with the target
            # pairs all assigned to G, and with only the outer pair assigned to G.
            # This checks fixed structures; it is NOT a nonempty-G PF comparison.
            if "target" in record and record["id"].startswith("ml_"):
                db = record["target"]
                pp = pairs(db)
                for label, fixed in [("all_pairs_in_G", pp), ("outer_pair_in_G", pp[:1])]:
                    g = db_from_pairs(n, fixed)
                    extension = db_from_pairs(n, set(pp)-set(fixed))
                    rec = dict(record, scaffold=g)
                    evals = {variant: invoke(binary, param, rec, "eval", 2, [extension], archive,
                                             rid+"/"+label+"_"+variant, trace=True)
                             for variant, binary in binaries.items()}
                    if any(r["status"] != "success" for r in evals.values()):
                        raise RuntimeError("nonempty G fixed-structure evaluation failed")
                    ref = structure_reference(sequence, db, 2)
                    a, b = (evals[k]["targets"][0] for k in binaries)
                    structure_rows.append(dict(id=rid, group=record["group"], length=n, structure=db, scaffold=g,
                        roles=label, dangles=2, boundary_multiloops=json.dumps(boundary_multiloops(db)),
                        reachable=a["same_pairs"] and a["reachable"], pinned_energy=a["energy"], diagnostic_energy=b["energy"],
                        delta_pinned_vs_vienna_pf=a["energy"]-ref["vienna_structure_pf_energy"],
                        delta_diagnostic_vs_vienna_pf=b["energy"]-ref["vienna_structure_pf_energy"],
                        delta_pinned_vs_vienna_eval=a["energy"]-ref["vienna_eval_energy"],
                        **{k: v for k, v in ref.items() if k != "vienna_loop_energies"}))
            if (index+1) % 10 == 0:
                print(json.dumps(dict(completed=index+1, total=len(records), id=rid,
                                      max_pinned_logz=max(abs(r.get("delta_logz", 0.)) for r in ensemble_rows if r["variant"] == "pinned"))), flush=True)
                write_tsv(args.output/"ensemble_metrics.tsv", ensemble_rows)
                write_tsv(args.output/"structure_metrics.tsv", structure_rows)
    write_tsv(args.output/"ensemble_metrics.tsv", ensemble_rows)
    write_tsv(args.output/"structure_metrics.tsv", structure_rows)
    summary = summarize(ensemble_rows, structure_rows)
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    write_tsv(args.output/"summary.tsv", summary["ensemble"])
    provenance.update(finished_utc=now(), dependencies_unchanged=all(digest(p) == hashes[str(p)] for p in dependencies),
                      dependencies_sha256_after={str(p): digest(p) for p in dependencies},
                      result_sha256={p.name: digest(p) for p in args.output.iterdir() if p.name != "provenance.json"})
    (args.output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    if not provenance["dependencies_unchanged"]:
        raise RuntimeError("dependencies changed during comparison")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
