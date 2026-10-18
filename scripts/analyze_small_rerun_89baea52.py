#!/usr/bin/env python3
"""Stream the new packed recovery/membership runs and export example data tars."""
import argparse
import io
import json
import os
from pathlib import Path
import statistics
import tarfile
import tempfile

import analyze_pk_membership as membership
from analyze_e3_neat1 import pairs_of
from prepare_paper_rerun_89baea52 import COMMIT, ROOT
from run_cli import digest, now, write_json
from task_storage import RunReader


def groups(name):
    run = ROOT / "results/runs" / (name+"-89baea52-20261018")
    plan = json.loads((run / "array_plan.json").read_text())
    assert plan["software_before"]["commit"] == COMMIT
    for task in plan["tasks"]:
        reader = RunReader(run / task["directory"] / "run")
        assert reader.manifest["status"] == "success"
        assert reader.manifest["software_before"]["commit"] == COMMIT
        assert reader.manifest["software_lock_valid_after_run"] is True
        for job in reader.manifest["jobs"]:
            if job["warmup"]:
                continue
            if job["status"] != "success":
                raise ValueError("failed/missing observation: "+job["id"])
            yield reader, job


def mfe_square(text, n):
    for line in text.splitlines():
        fields = line.split()
        token = fields[0] if fields else ""
        if len(fields) >= 2 and len(token) == n and set(token) <= set(".()[]{}<>"):
            return {(i,j) for i,j,o in membership.parse_structure(token) if o == "["}
    raise ValueError("MFE structure missing")


def recovery(out):
    cohorts = {"200": json.loads((ROOT / "results/local/pk-recovery-bprna-20261016/pk_sel.json").read_text()),
               "80": json.loads((ROOT / "results/local/pk-recovery-long-rrna-20261016/mid300_700_selection.json").read_text())}
    results = []
    for count,param,name in (("200","turner04","pk-recovery-200-turner04"),
                             ("200","dp09","pk-recovery-200-dp09"), ("80","turner04","pk-recovery-80")):
        selected = {r["id"]:r for r in cohorts[count]}
        seen = set()
        for reader,job in groups(name):
            assert job["program"] == "LinearCParty"
            r = selected[job["input"]]; assert r["id"] not in seen; seen.add(r["id"])
            p = {}
            for line in reader.artifact(job,"bpp.tsv").decode().splitlines():
                cols = line.split()
                if len(cols) == 3 and cols[0].isdigit():
                    p[int(cols[0]),int(cols[1])] = float(cols[2])
            pk = set(map(tuple,r["pk"]))
            square = mfe_square(reader.artifact(job,"stdout.txt").decode(),len(r["seq"]))
            results.append(dict(id=r["id"], cohort=count, par=param, n=len(r["seq"]),
                pk_probs=[p.get(tuple(pair),0.) for pair in r["pk"]], n_pk=len(pk),
                mfe_crossing=len(square), mfe_crossing_correct=len(square & pk),
                median_span=statistics.median(j-i for i,j in pk),
                source_revision=COMMIT, bpp_sha256=digest_bytes(reader.artifact(job,"bpp.tsv"))))
        assert seen == selected.keys()
    write_json(out / "pk_results.json", [r for r in results if r["cohort"] == "200"])
    write_json(out / "mid300_700_results.json", [r for r in results if r["cohort"] == "80"])
    summary = []
    for count,param in (("200","turner04"),("200","dp09"),("80","turner04")):
        selected = [r for r in results if r["cohort"] == count and r["par"] == param]
        reference_pairs = sum(r["n_pk"] for r in selected)
        predicted = sum(r["mfe_crossing"] for r in selected)
        matched = sum(r["mfe_crossing_correct"] for r in selected)
        row = dict(cohort=count, parameter=param, rnas=len(selected), reference_pairs=reference_pairs,
            reference_pairs_gt_half=sum(sum(p > .5 for p in r["pk_probs"]) for r in selected),
            rnas_mean_probability_gt_half=sum(statistics.mean(r["pk_probs"]) > .5 for r in selected),
            mfe_square_pairs=predicted, matched_square_pairs=matched,
            sensitivity=matched/reference_pairs, ppv=matched/predicted if predicted else None,
            distance_groups={})
        for label,condition in (("median_span_le_100",lambda r:r["median_span"] <= 100),
                                ("median_span_gt_100",lambda r:r["median_span"] > 100)):
            subset = [r for r in selected if condition(r)]
            row["distance_groups"][label] = dict(rnas=len(subset),
                rna_mean_probability=statistics.mean(statistics.mean(r["pk_probs"]) for r in subset) if subset else None)
        summary.append(row)
    write_json(out / "recovery_summary.json", summary)
    rrna = []
    old = ROOT / "results/local/pk-recovery-long-rrna-20261016/outputs.tar.gz"
    with tarfile.open(old) as tar:
        info = json.load(tar.extractfile("info.json"))
    write_json(out / "rrna_reference_metadata.json", info)
    for reader,job in groups("pk-recovery-rrna"):
        row = dict(input=job["input"], program=job["program"],
                   artifacts=job.get("artifacts"), status=job["status"])
        if job["program"] == "LinearCParty":
            p = {}
            for line in reader.artifact(job,"bpp.tsv").decode().splitlines():
                c = line.split()
                if len(c) == 3 and c[0].isdigit():
                    p[int(c[0]),int(c[1])] = float(c[2])
            row["reference_pk_pairs"] = info[job["input"]]["pk"]
            row["reference_probabilities"] = [p.get(tuple(pair),0.) for pair in row["reference_pk_pairs"]]
        rrna.append(row)
    write_json(out / "rrna_artifact_index.json", rrna)


def digest_bytes(data):
    import hashlib
    return hashlib.sha256(data).hexdigest()


def membership_summary(out):
    records = []
    for reader,job in groups("pk-membership"):
        r = json.loads(reader.artifact(job,"stdout.txt"))
        assert r["returncode"] == 0 and r["probability_source_revision"] == COMMIT
        records.append(r)
    cohort = json.loads((ROOT / "results/local/pk-recovery-bprna-20261016/pk_sel.json").read_text())
    assert len(records) == len(cohort) == 200
    previous = json.loads((ROOT / "results/local/small-rerun-89baea52-20261018/pk_results.json").read_text())
    assert all(r["source_revision"] == COMMIT for r in previous)
    membership.aggregate([membership.prepare(r) for r in cohort], records, previous, out,
                         revision="89baea52", previous_probability_source="same_release_recovery_exact")
    (out / "raw.jsonl").write_text("".join(json.dumps(r)+"\n" for r in records))


def examples(out):
    for group in ("terc", "figure1"):
        with tarfile.open(out / (group+"_data.tar"), "w") as tar:
            names = set()
            def add(name,data):
                assert name not in names; names.add(name)
                info = tarfile.TarInfo(name); info.size = len(data); tar.addfile(info,io.BytesIO(data))
            old = ROOT / ("results/local/figure4-terc-20261017" if group == "terc" else "results/figures/figure1_20261016/data")
            for n in ("rec.json","in.txt"):
                add(n,(old/n).read_bytes())
            mapping = {"LinearCParty":"lc.out","LinearPKCapR":"capr","LinearPKRaccess":"racc",
                       "LinearPKSampling":"samp","LinearPKMEA":"mea","LinearPKSubopt":"subopt"}
            for reader,job in groups(group):
                prefix = ("b200" if job["mode"] == "beam200_outer" else "exact")+"."
                ext = mapping[job["program"]]
                if group == "figure1":
                    prefix = ""; ext = {"capr":"capr.out","racc":"racc.out","samp":"samp.out","mea":"mea.out","subopt":"subopt.out"}.get(ext,ext)
                if group != "terc" or job["program"] != "LinearPKRaccess":
                    add(prefix+ext,reader.artifact(job,"stdout.txt"))
                if job["program"] == "LinearCParty":
                    add(prefix+"bpp.tsv",reader.artifact(job,"bpp.tsv"))
            if group == "terc":
                for reader,job in groups("terc-byloop"):
                    prefix = "b200" if job["mode"] == "beam200_outer" else "exact"
                    add(prefix+".racc",reader.artifact(job,"stdout.txt"))
            add("README.md",f"Release {COMMIT}. Raw outputs streamed from frozen packed plans; commands and hashes in the task JSONL. Turner04; Terc seed 20261017, Figure 1 seed 20261016.\n".encode())
    # Use the existing numerical parsers in local scratch. Rendering and
    # manuscript replacement remain separate from remeasuring these values.
    os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ["TMPDIR"])/"rerun-89baea52-mpl"))
    import make_figure4_terc as terc
    import make_figure1 as figure1
    with tempfile.TemporaryDirectory(prefix="example-summary-89baea52-", dir=os.environ["TMPDIR"]) as td:
        scratch = Path(td)
        for group in ("terc", "figure1"):
            local = scratch / group; local.mkdir()
            with tarfile.open(out / (group+"_data.tar")) as archive:
                for member in archive:
                    assert member.isfile() and Path(member.name).name == member.name
                    (local/member.name).write_bytes(archive.extractfile(member).read())
            if group == "terc":
                rec = json.loads((local/"rec.json").read_text())
                _, seq, scaffold = (local/"in.txt").read_text().splitlines()
                pk = set(map(tuple,rec["pk"]))
                stems = dict(outer={p for p in pk if p[0] <= 68}, inner={p for p in pk if p[0] >= 69}, both=pk)
                settings = {s:terc.read_setting(local,s,len(seq),scaffold) for s in ("exact","b200")}
                write_json(out/"terc_summary.json", dict(release=COMMIT,
                    settings={s:terc.summarize(v,stems,len(seq)) for s,v in settings.items()},
                    beam_exact_comparison=terc.compare(settings["b200"],settings["exact"])))
            else:
                sample_rows = [line.split() for line in (local/"samp.out").read_text().splitlines()]
                samples = sorted([(int(c[0]),c[1]) for c in sample_rows if len(c)==2 and c[0].isdigit()], reverse=True)
                summary = dict(release=COMMIT, total_draws=sum(c for c,_ in samples),
                               unique_structures=len(samples))
                assert summary["total_draws"] == 1000
                if len(samples) >= 2:
                    parsed = figure1.read_outputs(local)
                    summary.update(shown_sample_counts=[c for c,_ in parsed["shown"]],
                        other_draws=parsed["total"]-sum(c for c,_ in parsed["shown"]),
                        visible_pairs_strict_p_gt_0_01=sum(p>.01 for p in parsed["bpp"].values()),
                        pseudoloop_dominant_region=parsed["pseudo_region"], context_sum_error=parsed["context_error"])
                else:
                    summary["figure_layout_review"] = "Only one sampled structure; update the two-sample panel before rendering."
                write_json(out/"figure1_summary.json",summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind",choices=("recovery","membership","examples"))
    args = parser.parse_args()
    out = ROOT / ("results/local/pk-membership-89baea52-20261018" if args.kind == "membership" else "results/local/small-rerun-89baea52-20261018")
    if args.kind == "membership":
        out.mkdir(exist_ok=False)
    else:
        out.mkdir(exist_ok=True)
        target = "pk_results.json" if args.kind == "recovery" else "terc_data.tar"
        if (out/target).exists():
            raise ValueError("preserve the existing analysis")
    {"recovery": recovery,"membership": membership_summary,"examples": examples}[args.kind](out)
    write_json(out / (args.kind+"_provenance.json"),dict(release=COMMIT,created_utc=now(),script_sha256=digest(Path(__file__))))


if __name__ == "__main__":
    main()
