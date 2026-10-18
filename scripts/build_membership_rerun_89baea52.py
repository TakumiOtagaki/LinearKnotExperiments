#!/usr/bin/env python3
"""Build the small external membership observer; perform no RNA computation."""
import copy
import json
import shlex
import subprocess

from prepare_paper_rerun_89baea52 import BUILD, COMMIT, LOCK, ROOT, SOURCE, release_ready
from run_cli import build_settings, digest, now, write_json


def main():
    release_ready()
    output = ROOT / "build/membership-89baea52"
    lock_path = ROOT / "manifests/software-membership-89baea52-20261018.json"
    if output.exists() or lock_path.exists():
        raise ValueError("preserve the existing observer and lock")
    output.mkdir()
    cpp = ROOT / "scripts/pk_membership_driver.cc"
    libs = [BUILD / n for n in ("libcparty_scfg2_exact.a", "libRNA.a")]
    hashes = {str(p): digest(p) for p in libs}
    command = ["c++", "-O2", "-std=c++17", "-I"+str(SOURCE / "src"),
               "-I"+str(BUILD / "generated"), str(cpp), *map(str, libs),
               "-o", str(output / "membership_driver")]
    with (output / "build.log").open("w") as log:
        log.write(shlex.join(command)+"\n"); log.flush()
        subprocess.run(command, stdout=log, stderr=log, check=True)
    cohort = ROOT / "results/local/pk-recovery-bprna-20261016/pk_sel.json"
    wrapper = output / "E0Membership"
    wrapper.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec "${LCP_PYTHON:-python3}" ' +
        shlex.quote(str(ROOT / "scripts/rerun_membership_cli.py")) + " --driver " +
        shlex.quote(str(output / "membership_driver")) + " --cohort " + shlex.quote(str(cohort)) + ' "$@"\n')
    wrapper.chmod(0o755)
    (output / "CMakeCache.txt").write_text("CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_CXX_COMPILER:FILEPATH=/bin/c++\n")
    assert hashes == {str(p): digest(p) for p in libs}
    release_ready()
    write_json(output / "build_provenance.json", {"created_utc": now(), "command": command,
        "source_commit": COMMIT, "libraries_sha256": hashes, "base_lock_sha256": digest(ROOT / LOCK)})
    lock = copy.deepcopy(json.loads((ROOT / LOCK).read_text()))
    lock.update(created_utc=now(), purpose="formal", binaries={"E0Membership": digest(wrapper)},
                build_settings=build_settings(output), cmake_cache_sha256=digest(output / "CMakeCache.txt"),
                verification_method="External observation-only driver linked to certified release archives")
    extras = [*libs, cpp, cohort, output / "membership_driver", output / "build_provenance.json",
              ROOT / LOCK, ROOT / "scripts/rerun_membership_cli.py", ROOT / "scripts/analyze_pk_membership.py",
              ROOT / "scripts/rna_records.py", ROOT / "scripts/build_membership_rerun_89baea52.py"]
    lock["extra_files"] = [{"path": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in extras]
    write_json(lock_path, lock)


if __name__ == "__main__":
    main()
