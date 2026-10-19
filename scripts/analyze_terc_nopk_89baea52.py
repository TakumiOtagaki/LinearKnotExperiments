#!/usr/bin/env python3
"""Compare the frozen Terc exact ensemble with a G-constrained ViennaRNA ensemble.

Run from the experiments repository. All disposable files stay in TMPDIR.
No implementation checkout or manuscript file is changed.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import tempfile
import time

import RNA

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "deps/LinearCParty-89baea52"
BUILD = ROOT / "build/linear-cparty-89baea52"
OLD = ROOT / "results/local/small-rerun-89baea52-20261018/terc_data.tar"
FIG = ROOT / "results/figures/figure4_terc_89baea52_20261018"
CONTEXTS = ["Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop"]
LOOPS = dict(E="Exterior", H="Hairpin", B="Bulge", I="Internal", M="Multiloop", P="PseudoLoop")
HC = dict(Exterior=RNA.CONSTRAINT_CONTEXT_EXT_LOOP,
          Hairpin=RNA.CONSTRAINT_CONTEXT_HP_LOOP,
          Bulge_or_Internal=RNA.CONSTRAINT_CONTEXT_INT_LOOP,
          Multiloop=RNA.CONSTRAINT_CONTEXT_MB_LOOP)
REGIONS = dict(outer_left=range(63, 69), outer_right=range(134, 140),
               inner_left=range(69, 72), inner_right=range(130, 133))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def tsv(path, rows):
    with Path(path).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pairset(db):
    pt = RNA.ptable(db)
    return {(i, j) for i, j in enumerate(pt) if i and j > i}


def classify(db):
    """Assign every unpaired base to its immediate enclosing loop.

    A loop with one child is a bulge if exactly one side is empty, and an
    internal loop if both sides are nonempty. Pair endpoints are Stem.
    """
    pt = RNA.ptable(db)
    children = {0: []}
    stack = [0]
    for i, ch in enumerate(db, 1):
        if ch == "(":
            children[stack[-1]].append((i, pt[i]))
            children[i] = []
            stack.append(i)
        elif ch == ")":
            stack.pop()
    assert stack == [0]
    labels = [None] * (len(db) + 1)
    for i, branches in children.items():
        j = pt[i] if i else len(db) + 1
        if i == 0:
            context = "Exterior"
        elif not branches:
            context = "Hairpin"
        elif len(branches) > 1:
            context = "Multiloop"
        else:
            k, l = branches[0]
            context = "Bulge" if k == i + 1 or l == j - 1 else "Internal"
        k = i + 1
        for a, b in branches:
            labels[k:a] = [context] * (a - k)
            labels[a] = labels[b] = "Stem"
            k = b + 1
        labels[k:j] = [context] * (j - k)
    assert all(x in CONTEXTS for x in labels[1:])
    return labels


def make_md():
    md = RNA.md()
    md.temperature, md.betaScale, md.dangles = 37.0, 1.0, 2
    md.pf_smooth, md.min_loop_size = 1, 3
    md.noLP = md.noGU = md.noGUclosure = md.logML = md.circ = md.gquad = 0
    md.special_hp, md.compute_bpp, md.uniq_ML = 1, 1, 1
    md.salt = 1.021
    return md


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "results/local/terc-nopk-89baea52-20261019")
    parser.add_argument("--samples", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20261019)
    args = parser.parse_args()
    assert args.samples > 0
    assert RNA.__version__ == "2.7.2"
    scratch_root = Path(os.environ["TMPDIR"]).resolve()
    assert not str(scratch_root).startswith("/gs/bs/")
    args.out.mkdir(parents=True, exist_ok=False)
    for key in list(os.environ):
        if key.startswith("CPARTY_"):
            del os.environ[key]
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C")
    param = SOURCE / "params/rna_Turner04.par"
    lock = json.loads((ROOT / "manifests/software-pilot-89baea52.json").read_text())
    deps = [Path(__file__), param, OLD, FIG / "provenance.json", FIG / "rerun_provenance.json",
            ROOT / "manifests/software-pilot-89baea52.json", BUILD / "libcparty_scfg2_exact.a",
            BUILD / "libRNA.a", BUILD / "CMakeCache.txt"]
    deps += [BUILD / name for name in ("LinearCParty", "LinearPKCapR", "LinearPKRaccess")]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in deps}
    assert digest(param) == lock["energy_parameter_sha256"]
    for name in ("LinearCParty", "LinearPKCapR", "LinearPKRaccess"):
        assert digest(BUILD / name) == lock["binaries"][name]
    revision = subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True).strip()
    source_status = subprocess.check_output(["git", "-C", str(SOURCE), "status", "--porcelain"], text=True)
    assert revision == lock["source_commit"] and not source_status
    provenance = dict(started_utc=datetime.now(timezone.utc).isoformat(), argv=sys.argv,
                      source_revision=revision, source_status_before=source_status,
                      dependencies_sha256_before=hashes, python=sys.version,
                      vienna_version=RNA.__version__, vienna_module=RNA.__file__,
                      hostname=platform.node(), seed=args.seed, draws_per_ensemble=args.samples,
                      analysis_type="single-input conditional-ensemble comparison; no timing claims")
    write_json(args.out / "provenance.json", provenance)
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="terc-nopk-", dir=scratch_root) as td:
        scratch = Path(td)
        try:
            with tarfile.open(OLD) as archive:
                raw = {name: archive.extractfile(name).read() for name in
                       ("in.txt", "rec.json", "exact.capr", "exact.racc", "exact.bpp.tsv")}
            expected = json.loads((FIG / "rerun_provenance.json").read_text())["data_sha256"]
            for name, data in raw.items():
                assert hashlib.sha256(data).hexdigest() == expected[name]
            for name in ("in.txt", "rec.json"):
                (scratch / name).write_bytes(raw[name])
            seq, g = raw["in.txt"].decode().splitlines()[1:3]
            record = json.loads(raw["rec.json"])
            gpairs = pairset(g)
            n = len(seq)
            assert n == 388 and len(gpairs) == 61 and gpairs == set(map(tuple, record["nested"]))
            RNA.params_load(str(param))
            md = make_md()
            write_json(scratch / "model.json", dict(model_details=str(md), parameter_file=str(param),
                       parameter_sha256=digest(param), hard_constraint_flags=
                       RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP))
            calls = []

            def fold(positions=(), context=None):
                fc = RNA.fold_compound(seq, md)
                assert fc.hc_add_from_db(g, RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
                for i in positions:
                    fc.hc_add_up(i, (context or RNA.CONSTRAINT_CONTEXT_ALL_LOOPS) |
                                 RNA.CONSTRAINT_CONTEXT_ENFORCE)
                rounded = fc.pf()[1]
                q = fc.exp_matrices.q1k[n]
                logz = math.log(q) + n * math.log(fc.exp_params.pf_scale) if q > 0 else -math.inf
                calls.append(dict(positions=list(positions), context=context, q=q,
                                  pf_scale=fc.exp_params.pf_scale, logz=logz,
                                  free_energy=-fc.exp_params.kT / 1000 * logz,
                                  float_binding_free_energy=rounded))
                return fc, logz

            fc, logz = fold()
            rt = fc.exp_params.kT / 1000
            nofe = -rt * logz
            bpp = fc.bpp()
            forced_error = max(abs(bpp[i][j] - 1) for i, j in gpairs)
            assert forced_error < 1e-10
            # Two independently specified hard-constraint interfaces must agree.
            alt = RNA.fold_compound(seq, md)
            for i, j in gpairs:
                alt.hc_add_bp(i, j, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE)
            alt.pf()
            altlogz = math.log(alt.exp_matrices.q1k[n]) + n * math.log(alt.exp_params.pf_scale)
            constraint_logz_error = abs(logz - altlogz)
            assert constraint_logz_error < 1e-10
            nopair = {i: max(0., 1 - sum(bpp[min(i, j)][max(i, j)]
                                       for j in range(1, n + 1) if i != j)) for i in range(1, n + 1)}
            tsv(scratch / "nopk_bpp.tsv", [dict(i=i, j=j, probability=bpp[i][j])
                for i in range(1, n + 1) for j in range(i + 1, n + 1) if bpp[i][j] > 0])
            pkcontexts = {}
            for line in raw["exact.capr"].decode().splitlines():
                tokens = line.split()
                if tokens and tokens[0] in CONTEXTS:
                    assert len(tokens[1:]) == n
                    pkcontexts[tokens[0]] = dict(enumerate(map(float, tokens[1:]), 1))
            assert set(pkcontexts) == set(CONTEXTS)
            pkbpp = {}
            for line in raw["exact.bpp.tsv"].decode().splitlines():
                if line.startswith(">") or not line.strip():
                    continue
                a, b, value = line.split()
                pkbpp[int(a), int(b)] = float(value)
            meta = json.loads((FIG / "provenance.json").read_text())["raw_metadata"]["exact"]
            pkfe = float(meta["ensemble_free_energy"])
            assert abs(rt - float(meta["folding_rt"])) < 1e-14
            pkwindow = None
            for line in raw["exact.racc"].decode().splitlines():
                fields = line.split("\t")
                if fields[0] == "122":  # Raw Raccess starts at zero.
                    total = {int(x.split(",")[0]): float(x.split(",")[1])
                             for x in fields[1].split(";") if x}[4]
                    item = next(x for x in fields[2].removeprefix("byloop:").split(";") if x.startswith("4,"))
                    pkwindow = {LOOPS[k]: float(v) for k, v in
                                (x.split("=") for x in item.split(",")[1:])}
                    pkwindow["total"] = total
            assert pkwindow is not None

            # Check the loop classifier on unambiguous examples of each loop type.
            assert classify("..(((...)))..")[1] == "Exterior"
            assert classify("..(((...)))..")[6] == "Hairpin"
            assert classify("(.((...)))")[2] == "Bulge"
            assert classify("(.((...)).)")[2] == "Internal"
            assert classify("(.((...)).((...)).)")[2] == "Multiloop"
            counts = {i: Counter() for i in range(40, 161)}
            example_structures = []
            sample_stem_errors = []

            def sample(f, seed, label, positions=None):
                RNA.init_rand(seed)
                result = Counter()
                with (scratch / (label + ".samples.txt")).open("w") as stream:
                    for first in range(0, args.samples, 1000):
                        structures = f.pbacktrack(min(1000, args.samples - first))
                        assert len(structures) == min(1000, args.samples - first)
                        for db in structures:
                            assert len(db) == n and gpairs <= pairset(db)
                            stream.write(db + "\n")
                            labels = classify(db)
                            if positions is None:
                                for i in counts:
                                    counts[i][labels[i]] += 1
                                if len(example_structures) < 8:
                                    example_structures.append(db)
                            else:
                                assert len({labels[i] for i in positions}) == 1
                                result[labels[positions[0]]] += 1
                print(label + " sampling complete", flush=True)
                return result

            sample(fc, args.seed, "unconditional")
            # Exact E/H/(B+I)/M masses from separate constrained partition functions.
            # The inexpensive rare B/I split is sampled in the B+I conditional ensemble.
            noctx = {}
            for i in range(40, 161):
                if any(i in pair for pair in gpairs):
                    masses = {c: 0. for c in HC}
                else:
                    masses = {c: math.exp(fold([i], flag)[1] - logz) for c, flag in HC.items()}
                noctx[i] = masses
            _, windowlogz = fold(range(123, 127))
            access = math.exp(windowlogz - logz)
            window = {c: math.exp(fold(range(123, 127), flag)[1] - logz) for c, flag in HC.items()}
            conditional_counts = {}
            for label, positions, seed in (("base123_BI", [123], args.seed + 1),
                                           ("window123_126_BI", list(range(123, 127)), args.seed + 2)):
                conditional_fc, _ = fold(positions, RNA.CONSTRAINT_CONTEXT_INT_LOOP)
                conditional_counts[label] = sample(conditional_fc, seed, label, positions)
                assert set(conditional_counts[label]) <= {"Bulge", "Internal"}

            def split(mass, counter):
                draws = counter["Bulge"] + counter["Internal"]
                if mass == 0:
                    return 0., 0., 0., draws
                if not draws:
                    return None, None, None, draws
                fraction = counter["Bulge"] / draws
                return mass * fraction, mass * (1 - fraction), mass * math.sqrt(fraction * (1 - fraction) / draws), draws

            context_rows = []
            for ensemble in ("crossings_allowed_exact", "no_crossings_exact"):
                for i in range(40, 161):
                    if ensemble == "crossings_allowed_exact":
                        values = {c: pkcontexts[c][i] for c in CONTEXTS}
                        bi, se, draws = values["Bulge"] + values["Internal"], 0., 0
                    else:
                        masses = noctx[i]
                        counter = conditional_counts["base123_BI"] if i == 123 else counts[i]
                        bulge, internal, se, draws = split(masses["Bulge_or_Internal"], counter)
                        values = dict(Stem=1 - nopair[i], Hairpin=masses["Hairpin"],
                                      Bulge=bulge, Internal=internal, Multiloop=masses["Multiloop"],
                                      Exterior=masses["Exterior"], PseudoLoop=0.)
                        bi = masses["Bulge_or_Internal"]
                        sample_stem_errors.append(abs(counts[i]["Stem"] / args.samples - values["Stem"]))
                    context_rows.append(dict(ensemble=ensemble, position=i, nucleotide=seq[i - 1],
                        paired_exact=values["Stem"], unpaired_exact=1 - values["Stem"], **values,
                        Bulge_or_Internal_exact=bi, BI_split_standard_error=se, BI_split_draws=draws))
            tsv(args.out / "context_region.tsv", context_rows)
            bulge, internal, window_se, draws = split(window["Bulge_or_Internal"], conditional_counts["window123_126_BI"])
            window_rows = []
            for label, values, total, se, count in (
                ("crossings_allowed_exact", pkwindow, pkwindow["total"], 0., 0),
                ("no_crossings_exact", dict(Exterior=window["Exterior"], Hairpin=window["Hairpin"],
                    Bulge=bulge, Internal=internal, Multiloop=window["Multiloop"], PseudoLoop=0.),
                    access, window_se, draws)):
                window_rows.append(dict(ensemble=label, start=123, end=126, accessibility_exact=total,
                    **{c: values[c] for c in CONTEXTS if c != "Stem"},
                    Bulge_or_Internal_exact=values["Bulge"] + values["Internal"],
                    BI_split_standard_error=se, BI_split_draws=count))
            tsv(args.out / "window_accessibility.tsv", window_rows)
            energy_rows = [dict(ensemble=label, log_partition_function=z, free_energy_kcal_mol=fe, RT_kcal_mol=rt)
                          for label, z, fe in (("crossings_allowed_exact", -pkfe / rt, pkfe),
                                              ("no_crossings_exact", logz, nofe))]
            tsv(args.out / "ensemble.tsv", energy_rows)
            pair_rows = []
            for i, j in sorted(map(tuple, record["pk"])):
                cross = sorted((a, b) for a, b in gpairs if a < i < b < j or i < a < j < b)
                assert cross and bpp[i][j] == 0
                for label in ("crossings_allowed_exact", "no_crossings_exact"):
                    d = dict(ensemble=label, stem="outer" if i <= 68 else "inner", i=i, j=j,
                             pair_probability=pkbpp.get((i, j), 0.) if label.startswith("crossings") else bpp[i][j],
                             crossing_G_pairs=len(cross), crossing_G_example=f"{cross[0][0]}-{cross[0][1]}")
                    for side, pos in (("left", i), ("right", j)):
                        source = pkbpp if label.startswith("crossings") else {
                            (min(pos, k), max(pos, k)): bpp[min(pos, k)][max(pos, k)] for k in range(1, n + 1) if k != pos}
                        partners = [(b if a == pos else a, value) for (a, b), value in source.items() if pos in (a, b)]
                        partner, value = max(partners, key=lambda item: item[1])
                        d.update({side + "_paired_exact": pkcontexts["Stem"][pos] if label.startswith("crossings") else 1 - nopair[pos],
                                  side + "_top_partner": partner, side + "_top_partner_probability": value})
                    pair_rows.append(d)
            tsv(args.out / "reference_pairs.tsv", pair_rows)
            region_rows = []
            for ensemble in ("crossings_allowed_exact", "no_crossings_exact"):
                for region, positions in REGIONS.items():
                    items = [r for r in context_rows if r["ensemble"] == ensemble and r["position"] in positions]
                    fields = ["paired_exact", "unpaired_exact", *CONTEXTS, "Bulge_or_Internal_exact"]
                    region_rows.append(dict(ensemble=ensemble, region=region, start=min(positions), end=max(positions),
                        **{key: sum(r[key] for r in items) / len(items) if all(r[key] is not None for r in items) else None
                           for key in fields}))
            tsv(args.out / "region_summary.tsv", region_rows)
            write_json(scratch / "pf_calls.json", calls)
            write_json(scratch / "conditional_counts.json", conditional_counts)
            tsv(scratch / "unconditional_context_counts.tsv", [dict(position=i, **{c: counts[i][c] for c in CONTEXTS}) for i in counts])
            write_json(scratch / "exploratory_failure.json", dict(
                description="An initial interactive probe tried exp_matrices.get_Z(), which this Python binding does not expose.",
                error="AttributeError: 'mx_pf' object has no attribute 'get_Z'",
                resolution="Read exp_matrices.q1k[n] (double precision), as in the existing energy audit. No failed value is used."))
            help_runs = []
            for binary in ("LinearCParty", "LinearPKCapR", "LinearPKRaccess"):
                command = [str(BUILD / binary), "--help"]
                result = subprocess.run(command, capture_output=True, text=True, timeout=15)
                (scratch / (binary + ".help.txt")).write_text(result.stdout + result.stderr)
                help_runs.append(dict(argv=command, returncode=result.returncode))
            write_json(scratch / "help_runs.json", help_runs)
            print("Checking representative structure energies with the frozen public API", flush=True)
            binary = scratch / "public_probe"
            command = ["c++", "-O2", "-std=c++17", "-I" + str(SOURCE / "src"),
                       "-I" + str(ROOT / "scripts"), str(ROOT / "scripts/energy_vienna_public_89baea52.cc"),
                       str(BUILD / "libcparty_scfg2_exact.a"), str(BUILD / "libRNA.a"), "-o", str(binary)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=180)
            write_json(scratch / "build_probe.json", dict(argv=command, returncode=result.returncode,
                       stdout=result.stdout, stderr=result.stderr))
            result.check_returncode()
            mfe = RNA.fold_compound(seq, md)
            mfe.hc_add_from_db(g, RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
            mfe_db, mfe_energy = mfe.mfe()
            targets = [mfe_db] + example_structures
            generated = ["".join("." if i in {x for p in gpairs for x in p} else ch
                                  for i, ch in enumerate(db, 1)) for db in targets]
            command = [str(binary), str(param), seq, g, "2", "eval", "0"]
            result = subprocess.run(command, input="\n".join(generated) + "\n", capture_output=True, text=True, timeout=120)
            write_json(scratch / "energy_probe.json", dict(argv=command, stdin="\n".join(generated) + "\n",
                       returncode=result.returncode, stdout=result.stdout, stderr=result.stderr))
            result.check_returncode()
            evaluated = [line.split("\t") for line in result.stdout.splitlines() if line.startswith("EVAL\t")]
            assert len(evaluated) == len(targets)
            structure_rows = []
            for db, line in zip(targets, evaluated):
                check = RNA.fold_compound(seq, md)
                assert check.hc_add_from_db(db.replace(".", "x"), RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
                check.pf()
                z = math.log(check.exp_matrices.q1k[n]) + n * math.log(check.exp_params.pf_scale)
                ref = -rt * z
                assert line[2:4] == ["1", "1"]
                error = abs(float(line[4]) - ref)
                structure_rows.append(dict(structure=db, linearcparty_pf_energy=float(line[4]),
                                           vienna_fixed_structure_pf_energy=ref, abs_error=error))
            tsv(scratch / "energy_spot_checks.tsv", structure_rows)
            # The maximum single-child gap in G bounds internal loops closed by G
            # in any pseudoknot-free superset. Adding pairs can only shorten it.
            pt = RNA.ptable(g)
            gaps = []
            for i, j in sorted(gpairs):
                branches, k = [], i + 1
                while k < j:
                    if pt[k] > k:
                        branches.append((k, pt[k]))
                        k = pt[k] + 1
                    else:
                        k += 1
                if len(branches) == 1:
                    k, l = branches[0]
                    gaps.append(k - i - 1 + j - l - 1)
            context_error = max(abs((1 - nopair[i]) + sum(noctx[i].values()) - 1) for i in noctx)
            validation = dict(G_pairs=len(gpairs), length=n, forced_pair_max_abs_error=forced_error,
                constraint_interfaces_logz_abs_difference=constraint_logz_error,
                context_mass_max_abs_error=context_error,
                window_mass_abs_error=abs(sum(window.values()) - access),
                pk_window_mass_abs_error=abs(sum(pkwindow[c] for c in LOOPS.values()) - pkwindow["total"]),
                sampling_max_stem_abs_error=max(sample_stem_errors),
                classifier_fixtures_passed=5, samples_per_condition=args.samples,
                sampled_structure_G_violations=0, reference_pairs_forbidden_by_crossing_G=9,
                fixed_structure_energy_checks=len(structure_rows),
                fixed_structure_energy_max_abs_error=max(r["abs_error"] for r in structure_rows),
                max_single_child_unpaired_gap_in_G=max(gaps),
                temporary_files_location=str(scratch_root),
                delta_F_allowed_minus_no_crossings=pkfe - nofe,
                accessibility_allowed_minus_no_crossings=pkwindow["total"] - access)
            assert validation["fixed_structure_energy_max_abs_error"] < 1e-8
            assert max(gaps) <= 30
            assert max(context_error, validation["window_mass_abs_error"]) < 1e-9
            write_json(args.out / "validation.json", validation)
            # Keep a compact machine-readable summary for the accompanying README.
            write_json(scratch / "summary.json", dict(energy=energy_rows, window=window_rows,
                base123=[r for r in context_rows if r["position"] == 123], regions=region_rows,
                mfe_structure=mfe_db, mfe_integer_energy=mfe_energy, validation=validation))
            binary.unlink()  # Regenerable executable; preserve build command and sources.
        except BaseException as error:
            write_json(scratch / "failure.json", dict(type=type(error).__name__, message=str(error)))
            raise
        finally:
            with tarfile.open(args.out / "raw_outputs.tar.gz", "w:gz") as archive:
                for path in sorted(scratch.iterdir()):
                    if path.is_file():
                        archive.add(path, arcname=path.name)
                for name in ("energy_vienna_public_89baea52.cc", "energy_vienna_internal_89baea52.cc"):
                    archive.add(ROOT / "scripts" / name, arcname="sources/" + name)
                archive.add(Path(__file__), arcname="sources/" + Path(__file__).name)
            provenance.update(finished_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.monotonic() - start,
                dependencies_sha256_after={str(p.relative_to(ROOT)): digest(p) for p in deps},
                source_status_after=subprocess.check_output(["git", "-C", str(SOURCE), "status", "--porcelain"], text=True),
                raw_archive_sha256=digest(args.out / "raw_outputs.tar.gz"))
            write_json(args.out / "provenance.json", provenance)
    assert provenance["dependencies_sha256_after"] == hashes
    assert not provenance["source_status_after"]
    print(json.dumps(validation, indent=2), flush=True)


if __name__ == "__main__":
    main()
