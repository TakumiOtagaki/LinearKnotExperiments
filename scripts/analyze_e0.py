#!/usr/bin/env python3
"""Keep numerical parity and closure model differences as separate E0 results."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path

from rna_records import read_record, scaffold_pairs
from run_cli import digest, experiment_path, now, write_json


def read_observation(path, record):
    value = json.loads(path.read_text())
    for key in ("log_z", "free_energy", "rt", "mfe"):
        if not math.isfinite(value[key]):
            raise ValueError("nonfinite E0 observation")
    if value["rt"] <= 0 or abs(value["log_z"] + value["free_energy"]/value["rt"]) > 1e-7:
        raise ValueError("thermodynamic metadata inconsistent")
    pairs = {}
    for i,j,p in value["pairs"]:
        if not (type(i) is int and type(j) is int and 1 <= i < j <= len(record[1])) or (i,j) in pairs:
            raise ValueError("invalid BPP coordinates")
        if not math.isfinite(p) or not -1e-9 <= p <= 1+1e-9:
            raise ValueError("invalid BPP probability")
        pairs[i,j]=p
    for pair in scaffold_pairs(record[2]):
        if abs(pairs.get(pair,0)-1)>1e-8:
            raise ValueError("G does not have probability one")
    incident=[0.0]*len(record[1])
    for (i,j),p in pairs.items():
        incident[i-1]+=p;incident[j-1]+=p
    if any(p>1+1e-8 for p in incident):raise ValueError("incident BPP exceeds one")
    if any(incident[i]>1e-8 for i,c in enumerate(record[2]) if c in "xX"):
        raise ValueError("forced-unpaired base is paired")
    return {**value,"pair_map":pairs}


def compare(candidate, reference):
    if abs(candidate["rt"]-reference["rt"])>1e-12:raise ValueError("RT differs")
    keys=set(candidate["pair_map"])|set(reference["pair_map"])
    errors=[abs(candidate["pair_map"].get(k,0)-reference["pair_map"].get(k,0)) for k in keys]
    return {"delta_log_z":candidate["log_z"]-reference["log_z"],
            "delta_free_energy":candidate["free_energy"]-reference["free_energy"],
            "delta_mfe":candidate["mfe"]-reference["mfe"],
            "bpp_max_absolute_difference":max(errors,default=0),"bpp_l1_sum":sum(errors),
            "coordinate_union_count":len(keys)}


def analyze(run):
    manifest=json.loads((run/"run.json").read_text())
    if manifest["status"]=="running" or manifest.get("summary_only"):
        raise ValueError("finish and aggregate observations before E0 analysis")
    records={}
    for entry in manifest["inputs"]:
        path=Path(entry["copied_path"])
        if digest(path)!=entry["sha256"]:raise ValueError("input changed")
        records[entry["name"]]=read_record(path)
    config=manifest["config"]["analysis"]
    parsed={};errors={}
    jobs=[j for j in manifest["jobs"] if not j["warmup"]]
    for j in jobs:
        if j["status"]!="success":continue
        try:parsed[j["input"],j["program"]]=read_observation(run/j["directory"]/"stdout.txt",records[j["input"]])
        except (ValueError,KeyError,OSError) as exc:errors[j["input"],j["program"]]=str(exc)
    rows=[]
    for name in records:
        for program,kind in (("E0LegacyCompatible","numerical_parity"),("E0Default","closure_model_difference")):
            keys=[(name,program),(name,"E0CPartyOracle")]
            row={"input":name,"comparison_kind":kind,"program":program,"reference_program":"E0CPartyOracle"}
            if any(k in errors for k in keys):row.update(status="artifact_error",errors=[errors[k] for k in keys if k in errors])
            elif any(k not in parsed for k in keys):row["status"]="observation_unavailable"
            else:
                m=compare(parsed[keys[0]],parsed[keys[1]])
                parity=(abs(m["delta_log_z"])<=config["log_z_absolute_tolerance"] and
                        abs(m["delta_free_energy"])<=config["free_energy_absolute_tolerance"] and
                        abs(m["delta_mfe"])<=config["mfe_absolute_tolerance"] and
                        m["bpp_max_absolute_difference"]<=config["bpp_absolute_tolerance"])
                row.update(metrics=m,status="success" if kind=="closure_model_difference" or parity else "parity_failed")
                if kind=="closure_model_difference":
                    row["restriction_direction_consistent"]=m["delta_log_z"]<=config["log_z_absolute_tolerance"] and m["delta_mfe"]>=-config["mfe_absolute_tolerance"]
                    compatible=parsed.get((name,"E0LegacyCompatible"))
                    if compatible is not None:
                        row["default_minus_compatible"]=compare(parsed[keys[0]],compatible)
                        baseline=compare(compatible,parsed[keys[1]])
                        row["compatible_minus_oracle"]=baseline
                        row["closure_attribution_valid"]=(abs(baseline["delta_log_z"])<=config["log_z_absolute_tolerance"] and
                            abs(baseline["delta_mfe"])<=config["mfe_absolute_tolerance"] and
                            baseline["bpp_max_absolute_difference"]<=config["bpp_absolute_tolerance"])
                profile=parsed[keys[0]].get("contexts")
                if profile:
                    row["profile_sum_max_deviation"]=max(abs(sum(x)-1) for x in profile)
                    row["profile_independent_reference_available"]=False
            rows.append(row)
    return {"schema_version":1,"created_utc":now(),"run_manifest_sha256":digest(run/"run.json"),
            "analysis_script_sha256":digest(Path(__file__)),"planned_inputs":len(records),
            "job_status_counts":dict(Counter(j["status"] for j in jobs)),"observations":rows}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("run");parser.add_argument("--out",required=True)
    args=parser.parse_args();out=experiment_path(args.out)
    if out.exists():parser.error("use new analysis output")
    result=analyze(experiment_path(args.run));out.parent.mkdir(parents=True,exist_ok=True);write_json(out,result)
    print(json.dumps(dict(Counter(r["status"] for r in result["observations"]))))
    return int(any(r["status"]!="success" for r in result["observations"]))


if __name__=="__main__":raise SystemExit(main())
