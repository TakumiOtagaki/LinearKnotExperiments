#!/usr/bin/env python3
"""Render the adopted E2 layout with new 4f97fd3 measurements and SI checks."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/main_4f97fd3_20261021")
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        raise ValueError("use a fresh figure directory")
    subprocess.run([sys.executable, "-B", str(ROOT / "scripts/make_figures.py"),
                    "--review-layout", "--figure3-only", "--e2",
                    str(ROOT / "results/local/e2-4f97fd3-20261021"),
                    "--e2-revision", "4f97fd3", "--out", str(out)], check=True)
    # The author-approved output roles and seed explanations apply unchanged.
    caption = (ROOT / "results/figures/main_89baea52_20261020_fb/caption.md").read_text()
    caption += """

**Figure S3. Peak memory of independent calculation modes.**
各点は LinearKnot 4f97fd3 の独立したプロセス一回の peak RSS である。
各モードは準備と順方向計算を含むため、最大値の差はメモリ成分の内訳ではない。
All points are peak RSS from separate one-shot LinearKnot 4f97fd3 processes.
Each mode includes preparation and forward calculation; differences between peaks do not measure memory components.

**Figure S4. Effect of concurrent execution on measured resources.**
4,793、5,482、6,990 nt の同じ入力と出力で、1 process と最大9 processを比較する。
塗りつぶしは同時実行の対照、白抜きは本文の測定を、単一 process の中央値で割った比である。
These are matched 4,793-, 5,482-, and 6,990-nt inputs and requested outputs, with one or at most nine active processes.
Filled symbols divide the concurrent-control median by the single-process median; open symbols use the main resource measurements as numerator.
Hosts differ; these ratios describe the observations and do not correct the main measurements.
"""
    (out / "caption.md").write_text(caption)
    (out / "README.md").write_text("""# LinearKnot 4f97fd3 time and memory

The four PDF/PNG pairs replace the resource figure and its SI companions.
`figure3_main` retains the author-approved 178 x 55 mm layout.
`figure3_synthetic_all_cli_si` shows all five lengths for the seven main shapes;
`figure3_mode_memory_si` and `figure3_parallel_impact_si` show the stage and concurrency checks.
Only time/RSS measurements change. Accuracy and structure figures retain their adopted sources.
All measurement points use `results/local/e2-4f97fd3-20261021/`.
Paired Japanese/English captions are in `caption.md`; plot scalars and validation are retained.

Reproduce: `python3 -B scripts/plot_timing_4f97fd3.py` from experiments/.
""")
    from run_cli import digest
    manifest = json.loads((out / "manifest.json").read_text())
    manifest["timing_wrapper_sha256"] = digest(Path(__file__))
    for name in ("README.md", "caption.md"):
        manifest["derived_files"][name] = digest(out / name)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
