"""Boundary, stage recovery, deferred G, and independent oracle observation checks."""

from contextlib import redirect_stdout
import copy
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from analyze_e0 import compare, read_observation
from biological_cohort import LP_COMMIT, SETTINGS, expand_cohort
from plan_array import estimate, load_plan, prepare, task_jobs
from run_array_task import run_task
from run_cli import digest, plan_jobs, write_json
import scaffold_array
from analyze_accuracy import analyze as analyze_accuracy
import test_array


class FullPlanTests(unittest.TestCase):
    def setUp(self):
        self.fixture=test_array.ArrayTests(methodName="test_dry_run_freezes_inputs_lock_and_plan_without_execution")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.f=self.fixture

    def test_exact_boundary_and_calibrated_mode_cost(self):
        config=copy.deepcopy(self.f.config)
        config.update(warmups=0,repetitions=1,inputs=[{"name":"a","length":399},{"name":"b","length":400}])
        config["modes"][0]["length_range"]={"min_length":1,"max_length":399}
        jobs=plan_jobs(config)
        self.assertEqual(len(jobs),6)
        self.assertFalse(any(j["input"]=="b" and j["mode"]=="exact" for j in jobs))
        tier=copy.deepcopy(self.f.settings["tiers"][0])
        tier["cost_model"].update(safety_factor=2,per_mode={
            "exact":{"seconds_coefficient":.001,"length_exponent":2,"memory_base_mib":2,"memory_coefficient":.01},
            "beam":{"seconds_coefficient":.01,"memory_coefficient":.01}})
        config["memory_limit_mib"]=64
        self.assertAlmostEqual(estimate(399,config,tier)["estimated_seconds"],2*2*(399**2*.001+399*.01))
        self.assertAlmostEqual(estimate(400,config,tier)["estimated_seconds"],2*2*400*.01)

    def test_compact_schedule_and_tier_limits(self):
        self.f.settings["compact_plan"]=True
        self.f.settings["tiers"][0]["job_limits"]={"timeout_seconds":4,"memory_limit_mib":128}
        self.f.prepare(dry_run=True)
        task=self.f.plan["tasks"][0]
        self.assertNotIn("planned_jobs",task)
        self.assertEqual(len(task_jobs(self.f.plan,task)),24)
        c=json.loads(Path(task["config_path"]).read_text())
        self.assertEqual((c["timeout_seconds"],c["memory_limit_mib"]),(4,128))
        self.assertEqual(load_plan(self.f.plan_path),self.f.plan)
        Path(task["planned_jobs_path"]).write_text("[]\n")
        with self.assertRaisesRegex(ValueError,"frozen dependency"):load_plan(self.f.plan_path)

    def test_increasing_beam_checkpointed_resume_keeps_completed_failures(self):
        self.f.config.update(warmups=0,repetitions=1,manifest_checkpoint_jobs=1000,
            modes=[{"name":"beam50","args":["-b","50"]},{"name":"beam20","args":["-b","20"]}])
        self.f.settings["tiers"][0]["order_kind"]="beam_ascending"
        self.f.prepare()
        with patch.dict(os.environ,self.f.env),redirect_stdout(io.StringIO()):
            self.assertEqual(run_task(self.f.plan_path,"short",through_beam=20),130)
        directory=self.f.output/self.f.plan["tasks"][0]["directory"]/"run"
        manifest=json.loads((directory/"run.json").read_text())
        self.assertEqual({j["mode"] for j in manifest["jobs"]},{"beam20"})
        kept={j["id"]:(directory/j["directory"]/"result.json").read_bytes() for j in manifest["jobs"]}
        # Simulate a controller killed before its final checkpoint; result markers
        # and the append journal are enough to recover starts and failures.
        manifest["jobs"]=[];manifest["execution_order"]=[]
        write_json(directory/"run.json",manifest)
        with patch.dict(os.environ,self.f.env),redirect_stdout(io.StringIO()):
            self.assertEqual(run_task(self.f.plan_path,"short",resume=True,through_beam=50),1)
        resumed=json.loads((directory/"run.json").read_text())
        self.assertEqual(len(resumed["execution_order"]),8)
        self.assertEqual([j["mode"] for j in resumed["jobs"]].count("beam20"),4)
        for name,payload in kept.items():self.assertEqual((directory/"jobs"/name/"result.json").read_bytes(),payload)

    def test_deferred_scaffold_plan_never_runs_dummy_dots(self):
        self.f.config["inputs"][0]["scaffold_pending"]=True
        self.f.config["inputs"][0]["planning_sequence"]="A"*5
        self.f.prepare(dry_run=True)
        # Explicit pending state is normally set by expand_cohort; exercise the
        # runner's independent guard even if a manually constructed config bypasses it.
        with self.assertRaisesRegex(ValueError,"pending scaffolds"):
            self.f.run_task()
        self.assertFalse(self.f.events.exists())


class CohortAndScaffoldTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.fasta=self.root/"clean.fa"
        self.records=[("a","GGAAAACC"),("b","GGAAAAACC"),("c","GGAAAAAACC")]
        self.fasta.write_text("".join(">"+n+"\n"+s+"\n" for n,s in reversed(self.records)))
        self.acquisition=self.root/"source.json"
        write_json(self.acquisition,{"copied_path":str(self.fasta),"sha256":digest(self.fasta),"records":3})
        self.lp=self.root/"linearpartition"
        self.lp.write_text(f"#!{sys.executable}\nimport pathlib,sys\ns=sys.stdin.read().strip()\n"
                           "p=pathlib.Path(sys.argv[sys.argv.index('-r')+1])\n"
                           "p.write_text(f'1 {len(s)} 0.9\\n2 {len(s)-1} 0.8\\n')\n")
        self.lp.chmod(0o755);(self.root/"bin").mkdir();(self.root/"bin/linearpartition_v").write_text("fixture")
        self.lock=self.root/"lp.json"
        write_json(self.lock,{"source_commit":LP_COMMIT,"wrapper_sha256":digest(self.lp),
                              "vienna_binary_sha256":digest(self.root/"bin/linearpartition_v")})
        self.state={"commit":LP_COMMIT,"dirty":False}
        self.config={"schema_version":1,"purpose":"development","source_manifest":str(self.acquisition),
           "linearpartition":str(self.lp),"software_lock":str(self.lock),"timeout_seconds":2,
           "array":{"schema_version":1,"measurement":"accuracy","seed":1,"tiers":[{
             "name":"short","min_length":1,"max_length":100,"resource":"cpu_16","h_rt":"00:01:00",
             "target_seconds":10,"memory_budget_mib":24576,"workers":1,"max_inputs_per_task":2,
             "cost_model":{"basis":"tiny fixture","seconds_per_nt_condition":.1,"mib_per_nt":0}}]}}
        self.config_path=self.root/"config.json";write_json(self.config_path,self.config)
        self.output=self.root/"plan";self.jobs=self.root/"jobs"
        self.patches=[patch("scaffold_array.git_state",return_value=self.state),
                      patch("prepare_biological_inputs.git_state",return_value=self.state)]
        for p in self.patches:p.start();self.addCleanup(p.stop)
        env=patch.dict(os.environ,{k:v for k,v in os.environ.items() if k not in ("JOB_ID","SGE_TASK_ID")},clear=True)
        env.start();self.addCleanup(env.stop)

    def test_full_g_array_resume_and_verified_cohort_resolution(self):
        plan=scaffold_array.plan_array(self.config_path,self.output,self.jobs)
        self.assertEqual([len(t["inputs"]) for t in plan["tasks"]],[2,1])
        path=self.output/"array_plan.json"
        cohort=self.root/"cohort.json";dataset=self.root/"generated/dataset.json"
        write_json(cohort,{"source_manifest":str(self.acquisition),"source_manifest_sha256":digest(self.acquisition),
             "scaffold_dataset":str(dataset),"records":[{"id":n,"length":len(s),
                "sequence_sha256":hashlib.sha256(s.encode()).hexdigest()} for n,s in self.records]})
        config={"input_cohort":str(cohort),"inputs":[]}
        expanded,_,pending=expand_cohort(config,True)
        self.assertTrue(pending);self.assertTrue(all(i["scaffold_pending"] for i in expanded["inputs"]))
        with self.assertRaisesRegex(ValueError,"pending"):expand_cohort(config,False)
        with redirect_stdout(io.StringIO()):
            self.assertEqual(scaffold_array.run(path,"short",1),0)
            before=(self.output/plan["tasks"][0]["directory"]/"run/scaffold_logs/a/result.json").read_bytes()
            saved=self.output/plan["tasks"][0]["directory"]/"run/dataset.json"
            interrupted=json.loads(saved.read_text());interrupted["records"]=[];write_json(saved,interrupted)
            self.assertEqual(scaffold_array.run(path,"short",1,resume=True),0)
            self.assertEqual(scaffold_array.run(path,"short",2),0)
        self.assertEqual((self.output/plan["tasks"][0]["directory"]/"run/scaffold_logs/a/result.json").read_bytes(),before)
        summary=scaffold_array.aggregate(path,dataset.parent)
        self.assertEqual((summary["status"],summary["completed_jobs"]),("success",3))
        expanded,provenance,pending=expand_cohort(config,False)
        self.assertFalse(pending);self.assertIsNotNone(provenance["scaffold_dataset_sha256"])
        self.assertEqual(Path(expanded["inputs"][0]["path"]).read_text().splitlines()[2],"((....))")
        Path(expanded["inputs"][0]["path"]).write_text("changed\n")
        with self.assertRaises(ValueError):expand_cohort(config,False)

    def test_incomplete_g_never_publishes_a_reduced_dataset(self):
        scaffold_array.plan_array(self.config_path,self.output,self.jobs)
        with self.assertRaisesRegex(ValueError,"incomplete"):
            scaffold_array.aggregate(self.output/"array_plan.json",self.root/"generated")
        self.assertFalse((self.root/"generated").exists())


class E0ObservationTests(unittest.TestCase):
    def test_pair_union_and_forced_unpaired_validation(self):
        reference={"rt":.6,"log_z":1,"free_energy":-.6,"mfe":-.5,"pair_map":{(1,7):.3}}
        candidate={**reference,"log_z":.9,"free_energy":-.54,"pair_map":{(2,8):.2}}
        result=compare(candidate,reference)
        self.assertAlmostEqual(result["bpp_l1_sum"],.5);self.assertAlmostEqual(result["bpp_max_absolute_difference"],.3)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"out.json"
            write_json(p,{"log_z":1,"free_energy":-.6,"rt":.6,"mfe":-.5,"pairs":[[1,8,.3]]})
            with self.assertRaisesRegex(ValueError,"forced-unpaired"):
                read_observation(p,("a","GAAAAAAC","X......."))

    @unittest.skipUnless(shutil.which("c++"),"C++ compiler unavailable")
    def test_oracle_reverse_derivative_against_independent_two_structure_ensemble(self):
        source='''#include "legacy_ad.hh"
#include <cmath>
int main() {
  using namespace oracle_ad;
  const double rt=.616;
  auto a=tag(Real(2),1,6), b=tag(Real(3),2,7);
  auto z=a*b + Real(4); // alternatives with weights 6 and 4
  auto energy=-log(z)*rt;
  reverse(energy);
  for(auto t:tags) if(std::abs(-tape[t.node].adj*t.value/rt-.6)>1e-12) return 1;
  return 0;
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp);(path/"probe.cc").write_text(source)
            subprocess.run(["c++","-std=c++17","-I"+str(ROOT/"tools"),str(path/"probe.cc"),"-o",str(path/"probe")],check=True,capture_output=True)
            subprocess.run([str(path/"probe")],check=True,capture_output=True)


class LongSensitivityTests(unittest.TestCase):
    def test_long_rnas_use_each_cap_regime_and_largest_planned_width(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);record=root/"rna.txt";record.write_text(">rna\n"+"A"*400+"\n"+"."*400+"\n")
            config={"model_settings":{"temperature_celsius":37},"analysis":{
              "exact_max_length_exclusive":400,"long_reference_modes":{"default":"beam200_default","nocap":"beam200_nocap"}}}
            def save_run(run,modes):
                run.mkdir();jobs=[]
                for mode,status in modes:
                    directory=run/"jobs"/mode;directory.mkdir(parents=True)
                    bpp=directory/"bpp.tsv";bpp.write_text(">rna\n")
                    jobs.append({"id":mode,"input":"rna","mode":mode,"program":"LinearCParty","warmup":False,"repetition":0,
                         "status":status,"directory":"jobs/"+mode,"command":["LinearCParty","--bpp",str(bpp),"--bpp-cutoff","0"]})
                write_json(run/"run.json",{"schema_version":1,"status":"success","purpose":"development","config":config,
                    "inputs":[{"name":"rna","copied_path":str(record),"sha256":digest(record)}],"jobs":jobs,
                    "software_before":{"commit":"1"*40},"experiment_before":{"commit":"2"*40},"config_sha256":"3"*64,
                    "software_lock":{"energy_parameter_sha256":"4"*64},"software_changed_during_run":False,"binaries_changed_during_run":False})
            save_run(root/"main",[(m,"success") for m in ("beam20_default","beam200_default","beam20_nocap","beam200_nocap")])
            baseline=analyze_accuracy(root/"main","exact")
            self.assertTrue(all(o["comparison_kind"]=="beam_width_sensitivity" for o in baseline["observations"]))
            for row in baseline["observations"]:
                self.assertEqual(row["reference_mode"],"beam200_nocap" if row["mode"].endswith("nocap") else "beam200_default")
            save_run(root/"wider",[("beam400_default","failed"),("beam400_nocap","success")])
            wider=analyze_accuracy(root/"main","exact",[root/"wider"])
            default=next(o for o in wider["observations"] if o["mode"]=="beam20_default")
            self.assertEqual((default["reference_mode"],default["analysis_status"]),("beam400_default","reference_unavailable"))
            nocap=next(o for o in wider["observations"] if o["mode"]=="beam20_nocap")
            self.assertEqual((nocap["reference_mode"],nocap["analysis_status"]),("beam400_nocap","success"))


if __name__=="__main__":unittest.main()
