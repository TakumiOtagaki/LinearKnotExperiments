#!/usr/bin/env python3
"""Render the adopted non-E1 figures from the fresh 89baea52 analyses."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
REVISION = "89baea52"
DATE = "20261018"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(script, *args):
    subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / script),
                    *map(str, args)], cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", type=Path,
        default=ROOT / f"results/local/small-rerun-{REVISION}-{DATE}")
    parser.add_argument("--examples-only", action="store_true")
    parser.add_argument("--skip-neat1", action="store_true", help="Render ready examples and E2 while NEAT1 is pending")
    args = parser.parse_args()
    scratch_root = Path(os.environ["TMPDIR"]).resolve()
    assert not str(scratch_root).startswith("/gs/bs/")
    os.environ.update(PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1",
                      OPENBLAS_NUM_THREADS="1", MPLCONFIGDIR=str(scratch_root / "non-e1-89baea52-mpl"))
    with tempfile.TemporaryDirectory(prefix="figure-data-89baea52-", dir=scratch_root) as td:
        for group in ("figure1", "terc"):
            archive = args.examples / (group + "_data.tar")
            data = Path(td) / group
            data.mkdir()
            with tarfile.open(archive) as stream:
                for member in stream:
                    assert member.isfile() and Path(member.name).name == member.name
                    (data / member.name).write_bytes(stream.extractfile(member).read())
            name = "figure1" if group == "figure1" else "figure4_terc"
            out = ROOT / f"results/figures/{name}_{REVISION}_{DATE}"
            if group == "figure1":
                run("make_figure1.py", "--data", data, "--out", out)
            else:
                run("make_figure4_terc.py", "--data", data, "--out", out,
                    "--software-lock", ROOT / "manifests/software-pilot-89baea52.json",
                    "--identity-source", ROOT / "results/local/figure4-terc-20261017/README.md")
            # Keep the regenerated caption and validation, with durable raw-data
            # references rather than a temporary extraction path or old release.
            (out / "README.md").write_text(
                f"# {name}: {REVISION} ({DATE})\n\n"
                f"All probabilities, structures and sample counts were recomputed with {REVISION}.\n"
                f"Raw archive: `{archive.relative_to(ROOT) if archive.is_relative_to(ROOT) else archive}`.\n"
                "The archive contains the frozen input and new packed-task outputs.\n"
                "Only node-local TMPDIR is used for extraction.\n"
                "The caption and JSON report record the derived values and figure checks.\n\n"
                "Reproduce from the experiment repository:\n\n"
                "```bash\npython3 -B scripts/render_non_e1_89baea52.py\n```\n")
            metadata = dict(revision=REVISION, raw_archive=str(archive),
                raw_archive_sha256=sha(archive), renderer_sha256=sha(Path(__file__)),
                data_sha256={p.name: sha(p) for p in data.iterdir()},
                tar_extraction="node-local TMPDIR only")
            (out / "rerun_provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
            if group == "terc":
                provenance_path = out / "provenance.json"
                provenance = json.loads(provenance_path.read_text())
                for source in provenance["sources"]:
                    basename = Path(source["path"]).name
                    if basename in metadata["data_sha256"] and source["sha256"] == metadata["data_sha256"][basename]:
                        source["path"] = str(archive) + "::" + basename
                provenance["software"]["association_source"] = str(archive) + "::README.md"
                provenance["outputs"] = [dict(path=p.name, sha256=sha(p), bytes=p.stat().st_size)
                    for p in sorted(out.iterdir()) if p.is_file() and p.name != "provenance.json"]
                provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    if not args.examples_only:
        run("make_figures.py", "--review-layout", "--figure3-only",
            "--e2", ROOT / f"results/local/e2-{REVISION}-{DATE}",
            "--e2-revision", REVISION,
            "--out", ROOT / f"results/figures/main_{REVISION}_{DATE}")
        resource_out = ROOT / f"results/figures/main_{REVISION}_{DATE}"
        (resource_out / "README.md").write_text(
            f"# Time and memory: {REVISION}\n\n"
            f"All current measurements use `results/local/e2-{REVISION}-{DATE}/`.\n"
            "The renderer reads no E1 measurements in this invocation.\n"
            "Four figures give the main resource plot, independent mode memory, "
            "all synthetic curves, and the concurrency comparison.\n"
            "`plot_data.tsv` gives plotted values and ratio numerators/denominators.\n"
            "`manifest.json` gives source and script hashes; `validation.json` gives geometry and font checks.\n"
            "Main width is 178 mm; SI width is 165 mm; PNGs are 300 dpi.\n"
            "Each program requests the outputs specified in `caption.md`.\n"
            "Measured growth over these lengths does not establish an asymptotic order.\n\n"
            "Reproduce with `python3 -B scripts/render_non_e1_89baea52.py`.\n")
        manifest_path = resource_out / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["derived_files"]["README.md"] = sha(resource_out / "README.md")
        manifest["render_wrapper_sha256"] = sha(Path(__file__))
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        if not args.skip_neat1:
            run("make_figure4.py", "--data", ROOT / f"results/local/e3-neat1-{REVISION}-{DATE}",
                "--revision", REVISION,
                "--out", ROOT / f"results/figures/figure_s3_neat1_{REVISION}_{DATE}")


if __name__ == "__main__":
    main()
