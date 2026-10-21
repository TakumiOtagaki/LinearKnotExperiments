# LinearKnot 4f97fd3: resource audit and version comparison

12 matrices, 82 tasks, 2,768 processes: all successful, no missing observations.
E2 uses medians of three measurements after one warm-up; NEAT1 uses one measurement per condition.
Time below is GNU time wall seconds; peak RSS is GiB.
Old/new measurements share inputs and requested outputs but differ in source version, date, and host.
These are observed changes, not an isolated implementation speedup estimate.

| Input | Output | b | Old time s | New time s | Old RSS GiB | New RSS GiB |
|---|---|---:|---:|---:|---:|---:|
| URS00023A4DB9 | fold | 100 | 695.12 | 413.95 | 2.75 | 3.22 |
| URS00023A4DB9 | profile | 100 | 464.93 | 275.78 | 2.61 | 3.04 |
| URS00023A4DB9 | access | 100 | 469.74 | 284.40 | 2.61 | 3.04 |
| NEAT1_2 | bpp | 400 | 4,022.34 | 2,794.44 | 8.10 | 9.37 |
| NEAT1_2 | raccess-byloop | 100 | 398.48 | 246.88 | 2.38 | 2.80 |
| NEAT1_2 | sampling | 100 | 488.90 | 356.14 | 15.39 | 15.47 |
| NEAT1_2 | centroid | 400 | 4,760.99 | 3,900.44 | 9.20 | 9.44 |

`version_comparison.tsv` retains all 575 E2/NEAT1/stage/concurrency comparisons.
`pf_version_comparison.tsv` retains all 168 PF-only exact/b100 comparisons.
`neat1_resources.tsv` retains all 33 requested NEAT1 conditions, including profile, total/by-loop access, sampling, and centroid.
`pf_by_length.tsv` distinguishes paired-input medians from all-successful-input medians.
`run_audit.tsv` and `validation.json` retain counts and hashes.
All raw archives were read without filesystem extraction.

CParty has 0/21, 0/21, 12/21, 15/21, 21/21 nonfinite inputs at 250, 500, 1,000, 2,000, 4,000 nt.
On common successful inputs, CParty is faster through 1,000 nt; b100 is faster at 2,000 nt (23.40 vs 34.97 s).
All 51 finite CParty/exact comparisons pass the PF tolerances. Exact above 1,000 nt was excluded by plan.
At 4,000 nt no CParty input is finite; the b100-only full-cohort median is 47.58 s and 931.1 MiB.
The SI table explicitly marks this unpaired value.
