#!/usr/bin/env python3
"""Save small review artifacts for prepared bpRNA arrays, without execution."""

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import shutil
import subprocess

from plan_array import wall_seconds
from run_cli import ROOT, digest, experiment_path, mode_applies, write_json

RATES={"cpu_4":.015,"cpu_8":.03,"cpu_16":.06,"cpu_40":.15,"node_q":.25,"node_h":.5}


def report(plans, output):
    output.mkdir(parents=True,exist_ok=False)
    result={"schema_version":1,"preparation_only":True,"priority":0,
      "billing_source":"https://www.somuka.titech.ac.jp/reiki_int/reiki_honbun/x385RG00001693.html",
      "billing_checked_date":"2026-10-11",
      "point_formula":"sum(resource_coefficient * (0.7*max(estimated_task_wall_seconds,300)+0.1*h_rt_seconds)/3600), priority0. No parallel speedup or queue wait included. h_rt envelope substitutes h_rt for actual wall. Estimates omit per-charge 0.0001-point truncation, disk charges, and retries.",
      "plans":[]}
    for label,path in plans:
        plan=json.loads(path.read_text());groups=defaultdict(list)
        for task in plan["tasks"]:groups[task["tier"]].append(task)
        rows=[];assignment=["tier\ttask_id\tinputs\tjobs\testimated_seconds\testimated_memory_mib\tresource\th_rt\tinput_ids"]
        def count(task):
            if plan["run_kind"]=="scaffold_array":return len(task["inputs"])
            config=plan["config"]
            return sum(sum(mode_applies(m,i) for m in config["modes"]) for i in task["inputs"])*len(config["programs"])*(config.get("warmups",0)+config.get("repetitions",1))
        for name,tasks in groups.items():
            coefficient=RATES[tasks[0]["requested_resource"].split("=")[0]]
            seconds=sum(t["estimated_seconds"] for t in tasks)
            points=sum(coefficient*(.7*max(t["estimated_seconds"],300)+.1*wall_seconds(t["h_rt"]))/3600 for t in tasks)
            cap_points=sum(coefficient*.8*wall_seconds(t["h_rt"])/3600 for t in tasks)
            rows.append({"tier":name,"inputs":sum(len(t["inputs"]) for t in tasks),"tasks":len(tasks),
              "jobs":sum(count(t) for t in tasks),"estimated_hours":seconds/3600,"estimated_points":points,
              "h_rt_envelope_points":cap_points,"resource":tasks[0]["requested_resource"],"h_rt":tasks[0]["h_rt"],
              "workers":tasks[0]["workers"],"inputs_per_task_min":min(len(t["inputs"]) for t in tasks),
              "inputs_per_task_max":max(len(t["inputs"]) for t in tasks),
              "task_max_estimated_memory_mib":max(t["estimated_memory_mib"] for t in tasks)})
            for t in tasks:
                assignment.append("\t".join(map(str,[name,t["task_id"],len(t["inputs"]),count(t),t["estimated_seconds"],
                    t["estimated_memory_mib"],t["requested_resource"],t["h_rt"],",".join(i["name"] for i in t["inputs"])])))
        directory=output/label;directory.mkdir()
        (directory/"tasks.tsv").write_text("\n".join(assignment)+"\n")
        (directory/"array_plan.sha256").write_text(digest(path)+"\n")
        scripts=[]
        for script in plan["job_scripts"]:
            source=Path(script["path"])
            subprocess.run(["bash","-n",str(source)],check=True)
            shutil.copyfile(source,directory/source.name)
            scripts.append({"tier":script["tier"],"path":str(source),"sha256":digest(source)})
        totals={key:sum(r[key] for r in rows) for key in ("inputs","tasks","jobs","estimated_hours","estimated_points","h_rt_envelope_points")}
        entry={"label":label,"path":str(path),"array_plan_sha256":digest(path),
               "pending_scaffolds":plan.get("pending_scaffolds",False),"tiers":rows,"totals":totals,"scripts":scripts}
        write_json(directory/"summary.json",entry);result["plans"].append(entry)
        saved=json.loads((path.parent/("summary.json" if plan["run_kind"]=="scaffold_array" else "run.json")).read_text())
        expected=saved["planned_jobs"] if plan["run_kind"]=="scaffold_array" else saved["planned_measured_jobs"]
        if totals["jobs"]!=expected:raise ValueError("report differs from aggregate")
    write_json(output/"summary.json",result)
    for entry in result["plans"]:print(entry["label"],json.dumps(entry["totals"]))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",required=True)
    parser.add_argument("plans",nargs="+",help="label=plan-path")
    args=parser.parse_args()
    plans=[(item.split("=",1)[0],experiment_path(item.split("=",1)[1])) for item in args.plans]
    report(plans,experiment_path(args.out))


if __name__=="__main__":main()
