#!/usr/bin/env python3
"""Freeze cohorts and calibrate scheduling configs from saved pilot observations."""

import copy
import hashlib
import json
from pathlib import Path
import random

from biological_cohort import acquired_records
from run_cli import ROOT, digest, git_state, now, write_json

DATE = "20261011"
SOURCE = "manifests/data-sources/bprna-local-20261008.json"
DATASET = "results/local/bprna-full-scaffolds-20261011/dataset.json"
INTERVALS = [(50,99),(100,199),(200,299),(300,399),(400,799),(800,1599),(1600,3199),(3200,4381)]


def save(relative, value):
    path = ROOT / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ValueError("do not overwrite frozen preparation: " + relative)
    write_json(path, value)


def cohort(name, records, selection):
    relative = "manifests/cohorts/" + name + ".json"
    save(relative, {"schema_version": 1, "created_utc": now(), "source_manifest": SOURCE,
                   "source_manifest_sha256": digest(ROOT / SOURCE), "scaffold_dataset": DATASET,
                   "selection": selection, "records": [{"id": n, "length": len(s),
                       "sequence_sha256": hashlib.sha256(s.encode()).hexdigest()} for n,s in records]})
    return relative


def fit(rows, mode, low, high):
    selected = [r for r in rows if r["mode"] == mode and low <= r["length"] <= high]
    if not selected:
        selected = [r for r in rows if r["mode"] == mode]
    return {"constant_seconds": 0.1, "seconds_coefficient": max(1e-10, max(
                max(0,r["time_wall_seconds"]-0.1)/r["length"]**3 for r in selected)),
            "length_exponent": 3, "memory_base_mib": 128,
            "memory_coefficient": 1.5*max(1e-8,max(max(0,r["peak_rss_kib"]/1024-128)/r["length"]**3 for r in selected)),
            "memory_exponent": 3}


def tier(name, low, high, model, resource, cap, workers, target, wall, maximum):
    return {"name": name, "min_length": low, "max_length": high, "resource": resource,
            "h_rt": wall, "target_seconds": target, "memory_budget_mib":
                {"cpu_4":8192,"cpu_8":16384,"cpu_16":32768,"cpu_40":81920,
                 "node_q":172032,"node_h":327680}[resource],
            "task_overhead_mib":512, "workers":workers, "max_inputs_per_task":maximum,
            "job_limits":{"memory_limit_mib":cap,"timeout_seconds":600 if high<400 else 7200},
            "order_kind":"blocked_random" if high<400 else "beam_ascending", "cost_model":model}


def main():
    info, _, records = acquired_records(ROOT / SOURCE)
    records.sort(key=lambda r:(len(r[1]),r[0]))
    if len(records)!=26063 or sum(len(s)<400 for _,s in records)!=24625 or (len(records[0][1]),len(records[-1][1]))!=(50,4381):
        raise ValueError("author-approved clean cohort changed")
    full = cohort("bprna-full-20261011",records,{"kind":"all_clean_records","seed":None})
    rng=random.Random(20261011)
    subset=[]; pilot=[]
    for low,high in INTERVALS:
        pool=sorted([r for r in records if low<=len(r[1])<=high])
        chosen=rng.sample(pool, min(75 if high<400 else 3,len(pool)))
        (subset if high<400 else pilot).extend(chosen)
    subset.sort(key=lambda r:(len(r[1]),r[0]));pilot.sort(key=lambda r:(len(r[1]),r[0]))
    e0cohort=cohort("bprna-e0-subset-20261011",subset,{"kind":"75_per_short_tier","seed":20261011,
                     "pool_order":"ID sorted; Random.sample; final length/ID sorted"})
    pilotcohort=cohort("bprna-long-beam-pilot-20261011",pilot,{"kind":"3_per_long_tier","seed":20261011,
                        "pool_order":"same seeded generator after short E0 selection; ID sorted pools"})
    base=json.loads((ROOT/"configs/e1_bprna_15aa7742_20261011.json").read_text())
    run=ROOT/"results/runs/e1-bprna-15aa7742-20261011/run.json"
    measured=json.loads(run.read_text())
    lengths={i["name"]:len(Path(i["copied_path"]).read_text().splitlines()[1]) for i in measured["inputs"]}
    rows=[{k:j[k] for k in ("id","input","mode","program","status","time_wall_seconds","peak_rss_kib")}
          | {"length":lengths[j["input"]]} for j in measured["jobs"] if not j["warmup"]]
    if len(rows)!=324 or any(r["status"]!="success" for r in rows):raise ValueError("pilot incomplete")
    long_run=ROOT/"results/runs/e2-rnacentral-58a2b796-20261010/run.json"
    long_manifest=json.loads(long_run.read_text())
    ll={i["name"]:len(Path(i["copied_path"]).read_text().splitlines()[1]) for i in long_manifest["inputs"]}
    long_rows=[{k:j.get(k) for k in ("id","input","program","status","time_wall_seconds","peak_rss_kib")}
               | {"length":ll[j["input"]]} for j in long_manifest["jobs"] if not j["warmup"]]
    good=[r for r in long_rows if r["status"]=="success"]
    long_time=max(r["time_wall_seconds"]/r["length"] for r in good)
    long_memory=1.5*max(r["peak_rss_kib"]/1024/r["length"] for r in good)
    calibration={"schema_version":1,"created_utc":now(),"short_run":str(run.relative_to(ROOT)),
      "short_run_sha256":digest(run),"short_observations":rows,"long_run":str(long_run.relative_to(ROOT)),
      "long_run_sha256":digest(long_run),"long_observations_including_censored":long_rows,
      "method":"Per short tier and mode: upper envelope max((wall-0.1)/n^3), max((RSS_MiB-128)/n^3)*1.5. Pool all three CLIs; fallback to all 12 RNAs for unobserved tiers. Summed sequential work; safety factor 1.5. Long time: max successful wall/n at b100; scale by (b/100)^2. Long RSS uses the same b^2 planning envelope and 1.5 headroom. No fit to timeout durations; retain them as feasibility warnings. No speedup assumed.",
      "long_seconds_per_nt_b100":long_time,"long_mib_per_nt_b100_with_headroom":long_memory,
      "long_nocap_time_multiplier":2,"limitations":"12 RNAs, different scaffold geometries; concurrent accuracy times. RNAcentral 58a2b796 is a cross-version prior, not a 15aa7742 measurement. O(nb^2) with fixed caps does not bound uncapped search; x2 is provisional. Heavy-tail exact and AD oracle memory need task failures retained and pilot refit."}
    save("manifests/calibration/bprna-array-20261011.json",calibration)
    mainconfig=copy.deepcopy(base)
    mainconfig.update(name="e1-bprna-full-15aa7742-20261011",purpose="formal",inputs=[],input_cohort=full,
                      memory_limit_mib=65536,manifest_checkpoint_jobs=1000)
    mainconfig["modes"][0]["length_range"]={"min_length":50,"max_length":399}
    mainconfig["analysis"].update(exact_max_length_exclusive=400,
          long_reference_modes={"default":"beam200_default","nocap":"beam200_nocap"},
          long_comparison="Changes relative to largest measured width within each cap regime; no exact error, no monotonicity assumed.")
    save("configs/e1_bprna_full_15aa7742_20261011.json",mainconfig)
    tiers=[]
    resources=[("cpu_4",1024,4,3600,"02:00:00",500), ("cpu_8",2048,4,7200,"03:00:00",300),
               ("cpu_16",8192,2,7200,"03:00:00",150), ("cpu_40",16384,2,10800,"05:00:00",30),
               ("cpu_16",16384,1,7200,"03:00:00",8), ("cpu_40",32768,1,14400,"06:00:00",4),
               ("cpu_40",65536,1,18000,"08:00:00",2), ("cpu_40",65536,1,24000,"10:00:00",1)]
    for (low,high),(res,cap,workers,target,wall,maximum) in zip(INTERVALS,resources):
        per={}
        for mode in mainconfig["modes"]:
            if high<400:per[mode["name"]]=fit(rows,mode["name"],low,high)
            elif mode["name"]!="exact":
                b=int(mode["args"][1]);f=(b/100)**2
                per[mode["name"]]={"constant_seconds":0.1,"seconds_coefficient":long_time*f*(2 if mode["name"].endswith("nocap") else 1),
                    "length_exponent":1,"memory_base_mib":128,"memory_coefficient":long_memory*f,"memory_exponent":1}
        model={"basis":"manifests/calibration/bprna-array-20261011.json: observed envelopes plus safety, not a complexity fit",
               "seconds_per_nt_condition":1,"mib_per_nt":0,"safety_factor":1.5,"per_mode":per}
        tiers.append(tier(f"nt{low}_{high}",low,high,model,res,cap,workers,target,wall,maximum))
    save("configs/arrays/e1_bprna_full_20261011.json",{"schema_version":1,"measurement":"accuracy","seed":20261011,"compact_plan":True,"tiers":tiers})
    for width in (300,400):
        config=copy.deepcopy(mainconfig)
        config.update(name=f"e1-bprna-long-beam{width}-pilot-20261011",input_cohort=pilotcohort)
        config["modes"]=[]
        for group,limit in (("default",30),("nocap",100000)):
            config["modes"].append({"name":f"beam{width}_{group}","args":["-b",str(width),"--max-hairpin",str(limit),"--max-case3-depth",str(limit)]})
        config["analysis"].update(long_reference_modes={g:f"beam{width}_{g}" for g in ("default","nocap")},
            phase="feasibility_pilot",external_comparison_runs=["e1-bprna-full-15aa7742-20261011"]+(["e1-bprna-long-beam300-pilot-20261011"] if width==400 else []),
            expansion_gate="Run b300 after reviewing b20/50/100/200. Run b400 only after reviewing b300 per-bin success, wall and RSS. Expand full long cohort only after both pilots and revised cost budgets; largest completed width is the sensitivity reference, never exact.")
        pts=[]
        for low,high in INTERVALS[4:]:
            per={m["name"]:{"constant_seconds":0.1,"seconds_coefficient":long_time*(width/100)**2*(2 if m["name"].endswith("nocap") else 1),
                  "length_exponent":1,"memory_base_mib":128,"memory_coefficient":long_memory*(width/100)**2,"memory_exponent":1} for m in config["modes"]}
            model={"basis":"cross-version RNAcentral b100 prior scaled by b^2; provisional uncapped multiplier2",
                   "seconds_per_nt_condition":1,"mib_per_nt":0,"safety_factor":1.5,"per_mode":per}
            t=tier(f"nt{low}_{high}",low,high,model,"node_q" if width==300 else "node_h",
                   131072 if width==300 else 196608,1,40000 if width==300 else 60000,
                   "16:00:00" if width==300 else "24:00:00",1)
            t["job_limits"]["timeout_seconds"]=10800 if width==300 else 18000
            pts.append(t)
        save(f"configs/e1_bprna_long_beam{width}_pilot_20261011.json",config)
        save(f"configs/arrays/e1_bprna_long_beam{width}_pilot_20261011.json",{"schema_version":1,"measurement":"accuracy","seed":20261011,"tiers":pts})
    # Witness copies make E0 independent of the manuscript repository after preparation.
    witnessdir=ROOT/"data/e0-witnesses-20261011";witnessdir.mkdir(exist_ok=False)
    origin=ROOT.parent/"textbook/examples/verify_15aa7742"
    witness_sources=[];inputs=[]
    for name in ("closure14","nested45"):
        path=origin/"inputs"/(name+".txt")
        lines=path.read_text().splitlines()
        destination=witnessdir/(name+".txt")
        destination.write_text(">"+name+"\n"+lines[0]+"\n"+lines[1]+"\n")
        inputs.append({"name":name,"path":str(destination.relative_to(ROOT))})
        witness_sources.append({"id":name,"source":"textbook/examples/verify_15aa7742/inputs/"+path.name,"source_sha256":digest(path),"input_sha256":digest(destination)})
    longpath=origin/"inputs/long_cases.json"
    for entry in json.loads(longpath.read_text()):
        if entry["id"] not in ("BE41","BE65"):continue
        destination=witnessdir/(entry["id"]+".txt")
        destination.write_text(">"+entry["id"]+"\n"+entry["sequence"]+"\n"+entry["scaffold"]+"\n")
        inputs.append({"name":entry["id"],"path":str(destination.relative_to(ROOT))})
        witness_sources.append({"id":entry["id"],"source":"textbook/examples/verify_15aa7742/inputs/long_cases.json","source_sha256":digest(longpath),"input_sha256":digest(destination),"transformation":"none; preserve X forced-unpaired positions in both methods"})
    save("data/e0-witnesses-20261011/dataset.json",{"schema_version":1,"records":witness_sources})
    e0=copy.deepcopy(mainconfig)
    e0.update(name="e0-bprna-legacy-comparison-20261011",inputs=inputs,input_cohort=e0cohort,
              software_lock="manifests/software-e0-15aa7742-20261011.json",modes=[{"name":"exact"}],
              programs=[{"name":n} for n in ("E0LegacyCompatible","E0Default","E0CPartyOracle")])
    e0["analysis"]={"comparisons":{"compatible_vs_oracle":"Same seven-patch derivation set: logZ, BPP, MFE must agree within predeclared tolerance",
       "default_vs_oracle":"Model difference due to completed WMB weak closure; report differences, not numerical failures"},
       "log_z_absolute_tolerance":1e-7,"free_energy_absolute_tolerance":1e-7,"mfe_absolute_tolerance":1e-7,
       "bpp_absolute_tolerance":1e-8,"profile_scope":"Linear profiles for both policies only; upstream has no independent analytic profile API",
       "bpp_oracle":"Reverse derivative of upstream PF arithmetic, external observation copy, V/VP pair-emission tags; validate energy against uninstrumented oracle before arrays"}
    save("configs/e0_bprna_legacy_20261011.json",e0)
    et=[]
    intervals=[(1,49)]+INTERVALS[:4]
    for low,high in intervals:
        model={"basis":"E1 exact upper envelope x10 time, x6 memory for reverse-AD oracle; 14/45/41/65 nt smoke validates output only, not bpRNA resource scaling",
          "seconds_per_nt_condition":1,"mib_per_nt":0,"safety_factor":1.5,"per_mode":{"exact":fit(rows,"exact",max(50,low),high)}}
        f=model["per_mode"]["exact"];f["seconds_coefficient"]*=10;f["constant_seconds"]*=10;f["memory_coefficient"]*=6
        res,cap=("cpu_4",4096) if high<100 else ("cpu_8",8192) if high<200 else ("cpu_40",32768) if high<300 else ("node_q",131072)
        et.append(tier(f"nt{low}_{high}",low,high,model,res,cap,1,18000,"08:00:00",50 if high<200 else 10))
    save("configs/arrays/e0_bprna_legacy_20261011.json",{"schema_version":1,"measurement":"accuracy","seed":20261011,"tiers":et})
    dataset=json.loads((ROOT/"data/bprna-pilot-20261008/dataset.json").read_text())
    lp=ROOT/"deps/LinearPartition-20261008/linearpartition"
    lock={"schema_version":1,"source_commit":dataset["linearpartition_source"]["commit"],
          "source_git_state":git_state(lp.parent),"wrapper_sha256":digest(lp),
          "vienna_binary_sha256":digest(lp.parent/"bin/linearpartition_v"),
          "extra_files":[{"path":"scripts/"+name,"sha256":digest(ROOT/"scripts"/name)} for name in
                        ("prepare_biological_inputs.py","rna_records.py","scaffold_array.py","biological_cohort.py")]}
    save("manifests/software-scaffold-b450fb3e-20261011.json",lock)
    gt=[]
    for low,high in INTERVALS:
        model={"basis":"LP bprna elapsed envelope 0.05+0.0012n below400; RNAcentral successful elapsed envelope 0.05+0.007n above400, x1.5. See saved dataset hashes in calibration.",
               "seconds_constant_per_condition":0.05,"seconds_per_nt_condition":0.0012 if high<400 else 0.007,
               "memory_base_mib":128,"mib_per_nt":1,"safety_factor":1.5}
        t=tier(f"nt{low}_{high}",low,high,model,"cpu_4",8192,1,900,"01:00:00",500 if high<400 else 100)
        t["memory_budget_mib"]=8704
        t["job_limits"]["timeout_seconds"]=120
        gt.append(t)
    save("configs/scaffold_bprna_full_20261011.json",{"schema_version":1,"purpose":"formal","source_manifest":SOURCE,
         "linearpartition":"deps/LinearPartition-20261008/linearpartition","software_lock":"manifests/software-scaffold-b450fb3e-20261011.json",
         "scaffold_dataset":DATASET,"timeout_seconds":120,"array":{"schema_version":1,"measurement":"accuracy","seed":20261011,"tiers":gt}})
    print(json.dumps({"full":len(records),"exact":sum(len(s)<400 for _,s in records),"e0_subset":len(subset),"long_pilot":len(pilot),
          "long_seconds_per_nt_b100":long_time,"long_memory_per_nt_b100":long_memory},indent=2))


if __name__=="__main__":main()
