#!/usr/bin/env python3
"""Freeze the NEAT1 RefSeq transcripts and build E3 scaffolds.

The scaffold uses the same rule as E1/E2: LinearPartition-V (b=200, dangles 2)
posterior pairs with probability >= 0.7, keeping helices of at least two pairs
(prepare_biological_inputs.threshold_scaffold). Sequences are read from FASTA
files that the caller downloaded from NCBI; their hashes and the download URLs
are recorded. Writes one input file per transcript and one manifest.
"""
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_biological_inputs import threshold_scaffold  # noqa: E402

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text"
TRANSCRIPTS = {
    "NEAT1_1": "NR_028272.1",  # MENepsilon, short isoform
    "NEAT1_2": "NR_131012.1",  # MENbeta, long isoform
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def read_fasta(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0]
    sequence = "".join(line.strip() for line in lines[1:]).upper().replace("T", "U")
    if set(sequence) - set("ACGU"):
        raise ValueError(f"{path}: non-ACGU symbols")
    return header, sequence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-dir", type=Path, required=True)
    parser.add_argument("--linearpartition", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--retrieved-utc", required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("output directory exists; use a new one")
    lp_root = args.linearpartition.resolve().parent
    lp_commit = subprocess.run(["git", "-C", str(lp_root), "rev-parse", "HEAD"],
                               capture_output=True, text=True, check=True).stdout.strip()
    lp_dirty = subprocess.run(["git", "-C", str(lp_root), "status", "--short"],
                              capture_output=True, text=True, check=True).stdout.strip()
    if lp_dirty:
        raise SystemExit("LinearPartition checkout is not clean")
    args.out.mkdir(parents=True)
    manifest = {"schema_version": 1, "dataset": "NEAT1 RefSeq transcripts for E3",
                "retrieved_utc": args.retrieved_utc,
                "scaffold_method": "LinearPartition-V posterior threshold, no MEA",
                "scaffold_settings": {"beam": 200, "dangles": 2, "pair_probability_threshold": 0.7, "min_helix": 2},
                "linearpartition": {"path": str(args.linearpartition), "commit": lp_commit,
                                    "binary_sha256": sha256(args.linearpartition.read_bytes())},
                "records": []}
    for name, accession in TRANSCRIPTS.items():
        fasta = args.fasta_dir / f"{accession}.fa"
        header, sequence = read_fasta(fasta)
        with tempfile.TemporaryDirectory() as tmp:
            bpp = Path(tmp) / "pairs.tsv"
            command = [str(args.linearpartition.resolve()), "-V", "-b", "200", "--dangles", "2",
                       "-c", "0.7", "-r", str(bpp)]
            started = time.perf_counter()
            subprocess.run(command, input=sequence + "\n", capture_output=True, text=True, check=True)
            seconds = time.perf_counter() - started
            scaffold = threshold_scaffold(sequence, bpp, 0.7, 2)
        record_path = args.out / f"{name}.txt"
        record_path.write_text(f">{name}\n{sequence}\n{scaffold}\n")
        manifest["records"].append({
            "name": name, "accession": accession, "fasta_header": header,
            "source_url": EFETCH.format(accession), "fasta_sha256": sha256(fasta.read_bytes()),
            "length": len(sequence), "sequence_sha256": sha256(sequence.encode()),
            "scaffold_sha256": sha256(scaffold.encode()), "scaffold_pairs": scaffold.count("("),
            "path": record_path.name, "linearpartition_seconds": round(seconds, 3),
            "linearpartition_command": command[:-1] + ["<tmp>/pairs.tsv"]})
    a = read_fasta(args.fasta_dir / f"{TRANSCRIPTS['NEAT1_1']}.fa")[1]
    b = read_fasta(args.fasta_dir / f"{TRANSCRIPTS['NEAT1_2']}.fa")[1]
    manifest["isoform_shared_prefix_nt"] = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({r["name"]: [r["length"], r["scaffold_pairs"]] for r in manifest["records"]}))


if __name__ == "__main__":
    main()
