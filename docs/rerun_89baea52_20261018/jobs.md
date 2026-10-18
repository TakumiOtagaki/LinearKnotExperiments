# 89baea52 の job script 一覧

計画時間は既存係数・保守的なallowanceで、新版の計測値ではない。逐次仕事量はtask内の並列speedupを仮定しない。
旧task wallは同じ前回tierの全attemptの実時間を足した参考値で、stage間の待ち・queue待ちは含めない。失敗・再開attemptも含み、終点不明のattemptは0件。[元データとhash](historical_wall_estimates.json)。新しい小規模matrixには旧task metadataがないため — とする。

| Script (experiments/ から) | task数 | 資源 | workers | cap MiB/process | h_rt | 計画逐次仕事量 h | 計画最長task h | 旧task wall 合計/最長 h | 実行 |
|---|---:|---|---:|---:|---|---:|---:|---|---|
| `jobs/e0-bprna-legacy-89baea52-20261018/nt1_49.sh` | 1 | cpu_16=1 | 1 | 4096 | 08:00:00 | 0.006 | 0.006 | 0.000 / 0.000 | array |
| `jobs/e0-bprna-legacy-89baea52-20261018/nt50_99.sh` | 2 | cpu_16=1 | 1 | 4096 | 08:00:00 | 0.124 | 0.077 | 0.007 / 0.004 | array |
| `jobs/e0-bprna-legacy-89baea52-20261018/nt100_199.sh` | 2 | cpu_16=1 | 1 | 8192 | 08:00:00 | 0.828 | 0.438 | 0.019 / 0.010 | array |
| `jobs/e0-bprna-legacy-89baea52-20261018/nt200_299.sh` | 8 | cpu_40=1 | 1 | 32768 | 08:00:00 | 3.163 | 0.626 | 0.080 / 0.016 | array |
| `jobs/e0-bprna-legacy-89baea52-20261018/nt300_399.sh` | 8 | cpu_40=1 | 1 | 80000 | 08:00:00 | 19.994 | 3.271 | 0.286 / 0.046 | array |
| `jobs/e1-bprna-short-89baea52-20261018/nt50_99.sh` | 25 | cpu_16=1 | 2 | 8192 | 02:00:00 | 12.593 | 0.522 | 3.453 / 0.191 | array |
| `jobs/e1-bprna-short-89baea52-20261018/nt100_199.sh` | 34 | cpu_16=1 | 2 | 8192 | 03:00:00 | 41.487 | 1.999 | 7.639 / 0.436 | array |
| `jobs/e1-bprna-short-89baea52-20261018/nt200_299.sh` | 10 | cpu_16=1 | 2 | 16384 | 03:00:00 | 10.358 | 1.636 | 4.555 / 0.684 | array |
| `jobs/e1-bprna-short-89baea52-20261018/nt300_399.sh` | 37 | cpu_16=1 | 2 | 16384 | 08:00:00 | 81.658 | 3.183 | 10.659 / 0.449 | array |
| `jobs/e1-bprna-long-89baea52-20261018/nt400_799.sh` | 138 | cpu_40=1 | 1 | 65536 | 06:00:00 | 505.762 | 4.163 | 28.919 / 1.034 | array |
| `jobs/e1-bprna-long-89baea52-20261018/nt800_1599.sh` | 193 | cpu_40=1 | 1 | 65536 | 08:00:00 | 958.304 | 5.446 | 139.278 / 1.852 | array |
| `jobs/e1-bprna-long-89baea52-20261018/nt1600_3199.sh` | 273 | cpu_40=1 | 1 | 65536 | 12:00:00 | 1467.778 | 8.315 | 252.613 / 3.612 | array |
| `jobs/e1-bprna-long-89baea52-20261018/nt3200_4381.sh` | 90 | cpu_40=1 | 1 | 65536 | 16:00:00 | 562.054 | 10.983 | 115.991 / 4.231 | array |
| `jobs/e2-rnacentral-89baea52-20261018/rna_4_8kb.sh` | 1 | cpu_40=1 | 10 | 8192 | 12:00:00 | 3.171 | 3.171 | 0.498 / 0.498 | array |
| `jobs/e2-rnacentral-89baea52-20261018/rna_8_20kb.sh` | 1 | cpu_40=1 | 10 | 8192 | 12:00:00 | 4.560 | 4.560 | 1.006 / 1.006 | array |
| `jobs/e2-rnacentral-89baea52-20261018/rna_20_30kb.sh` | 2 | cpu_40=1 | 10 | 8192 | 12:00:00 | 9.139 | 5.144 | 1.995 / 1.015 | array |
| `jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n250.sh` | 1 | cpu_40=1 | 10 | 8192 | 12:00:00 | 3.375 | 3.375 | 0.062 / 0.062 | array |
| `jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n500.sh` | 1 | cpu_40=1 | 10 | 8192 | 12:00:00 | 6.188 | 6.188 | 0.302 / 0.302 | array |
| `jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n1000.sh` | 5 | cpu_40=1 | 10 | 8192 | 12:00:00 | 33.750 | 7.875 | 1.026 / 0.255 | array |
| `jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n2000.sh` | 12 | cpu_40=1 | 10 | 8192 | 12:00:00 | 78.300 | 6.525 | 3.850 / 0.376 | array |
| `jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n4000.sh` | 36 | cpu_40=1 | 10 | 8192 | 12:00:00 | 153.900 | 4.275 | 12.097 / 0.499 | array |
| `jobs/e2-memory-breakdown-89baea52-20261018/rna_4_8kb.sh` | 1 | cpu_40=1 | 10 | 8192 | 12:00:00 | 4.572 | 4.572 | 0.058 / 0.058 | array |
| `jobs/e2-memory-breakdown-89baea52-20261018/rna_8_20kb.sh` | 1 | cpu_40=1 | 10 | 8192 | 12:00:00 | 8.023 | 8.023 | 0.117 / 0.117 | array |
| `jobs/e2-memory-breakdown-89baea52-20261018/rna_20_30kb.sh` | 2 | cpu_40=1 | 10 | 8192 | 12:00:00 | 11.049 | 5.538 | 0.275 / 0.140 | array |
| `jobs/e2-parallel-impact-w1-89baea52-20261018/parallel_impact.sh` | 1 | cpu_40=1 | 1 | 8192 | 12:00:00 | 3.171 | 3.171 | 1.375 / 1.375 | array |
| `jobs/e2-parallel-impact-w10-89baea52-20261018/parallel_impact.sh` | 1 | cpu_40=1 | 10 | 8192 | 12:00:00 | 3.171 | 3.171 | 0.519 / 0.519 | array |
| `jobs/e3-neat1-bpp-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 1.590 | 1.590 | 0.149 / 0.149 | array |
| `jobs/e3-neat1-bpp-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 9.501 | 9.501 | 1.138 / 1.138 | array |
| `jobs/e3-neat1-capr-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 1.590 | 1.590 | 0.093 / 0.093 | array |
| `jobs/e3-neat1-capr-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 9.501 | 9.501 | 0.767 / 0.767 | array |
| `jobs/e3-neat1-mea-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 1.590 | 1.590 | 0.175 / 0.175 | array |
| `jobs/e3-neat1-mea-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 9.501 | 9.501 | 1.379 / 1.379 | array |
| `jobs/e3-neat1-raccess-byloop-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 1.590 | 1.590 | 0.096 / 0.096 | array |
| `jobs/e3-neat1-raccess-byloop-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 9.501 | 9.501 | 0.740 / 0.740 | array |
| `jobs/e3-neat1-raccess-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 1.590 | 1.590 | 0.094 / 0.094 | array |
| `jobs/e3-neat1-raccess-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 9.501 | 9.501 | 0.735 / 0.735 | array |
| `jobs/e3-neat1-sampling-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 2 | 40960 | 24:00:00 | 0.434 | 0.434 | 0.049 / 0.049 | array |
| `jobs/e3-neat1-sampling-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 2 | 40960 | 24:00:00 | 1.272 | 1.272 | 0.139 / 0.139 | array |
| `jobs/e3-neat1-variants-bpp-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 3.180 | 3.180 | 0.307 / 0.307 | array |
| `jobs/e3-neat1-variants-bpp-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 6.334 | 6.334 | 2.310 / 2.310 | array |
| `jobs/e3-neat1-variants-capr-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 3.180 | 3.180 | 0.195 / 0.195 | array |
| `jobs/e3-neat1-variants-capr-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 6.334 | 6.334 | 1.564 / 1.564 | array |
| `jobs/e3-neat1-variants-raccess-89baea52-20261018/neat1_short.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 3.180 | 3.180 | 0.197 / 0.197 | array |
| `jobs/e3-neat1-variants-raccess-89baea52-20261018/neat1_long.sh` | 1 | cpu_40=1 | 3 | 20480 | 24:00:00 | 6.334 | 6.334 | 1.494 / 1.494 | array |
| `jobs/pk-recovery-200-turner04-89baea52-20261018/pk200.sh` | 20 | cpu_16=1 | 1 | 8192 | 04:00:00 | 33.333 | 1.667 | — | array |
| `jobs/pk-recovery-200-dp09-89baea52-20261018/pk200.sh` | 20 | cpu_16=1 | 1 | 8192 | 04:00:00 | 33.333 | 1.667 | — | array |
| `jobs/pk-recovery-80-89baea52-20261018/pk80.sh` | 8 | cpu_16=1 | 1 | 8192 | 04:00:00 | 13.333 | 1.667 | — | array |
| `jobs/pk-recovery-rrna-89baea52-20261018/rrna.sh` | 3 | cpu_16=1 | 1 | 16384 | 04:00:00 | 6.000 | 2.000 | — | array |
| `jobs/terc-89baea52-20261018/terc.sh` | 1 | cpu_16=1 | 1 | 8192 | 04:00:00 | 0.400 | 0.400 | — | local |
| `jobs/terc-byloop-89baea52-20261018/terc.sh` | 1 | cpu_16=1 | 1 | 8192 | 04:00:00 | 0.067 | 0.067 | — | local |
| `jobs/figure1-89baea52-20261018/figure1.sh` | 1 | cpu_16=1 | 1 | 8192 | 04:00:00 | 0.200 | 0.200 | — | local |
| `jobs/pk-membership-89baea52-20261018/membership.sh` | 20 | cpu_16=1 | 1 | 8192 | 04:00:00 | 50.000 | 2.500 | — | array |
| `jobs/e0-validation-89baea52-20261018/enumeration.sh` | 1 | cpu_16=1 | 1 | 16384 | 02:00:00 | 0.500 | 0.500 | — | array |
| `jobs/e0-validation-89baea52-20261018/decoding.sh` | 1 | cpu_16=1 | 1 | 16384 | 02:00:00 | 0.500 | 0.500 | — | array |
| `jobs/e0-validation-89baea52-20261018/byloop.sh` | 1 | cpu_16=1 | 1 | 16384 | 02:00:00 | 0.250 | 0.250 | — | array |
| `jobs/e0-validation-89baea52-20261018/vienna.sh` | 1 | cpu_16=1 | 1 | 16384 | 02:00:00 | 0.750 | 0.750 | — | array |
| `jobs/e0-validation-89baea52-20261018/linearspace.sh` | 1 | cpu_16=1 | 1 | 16384 | 02:00:00 | 0.100 | 0.100 | — | array |


## E1 long cpu_16 再計画（2026-10-18）

旧表は保存済みcpu_40計画の記録として残す。比較元の実scriptのh_rtは6/8/12/16時間である。
新runは `e1-bprna-long-89baea52-cpu16-20261018`。全区分8 workers、4,096 MiB RLIMIT_AS/process、33,792 MiBのcap合計。
新taskのh_rtは最大4時間、推定wallは各h_rtの70%以内。以下の時間はallocation枠であり、異なる資源間のポイント換算ではない。

| 区分 | 旧cpu_40 tasks | 旧h_rt | 旧予約h | 新cpu_16 tasks | 新h_rt | 新予約h | 新推定最長task h |
|---|---:|---|---:|---:|---|---:|---:|
| nt400_799 | 138 | 06:00:00 | 828 | 7 | 01:00:00 | 7 | 0.6999 |
| nt800_1599 | 193 | 08:00:00 | 1,544 | 19 | 02:00:00 | 38 | 1.3971 |
| nt1600_3199 | 273 | 12:00:00 | 3,276 | 24 | 03:00:00 | 72 | 2.0997 |
| nt3200_4381 | 90 | 16:00:00 | 1,440 | 10 | 04:00:00 | 40 | 2.7975 |
| 合計 | 694 | — | 7,088 | 60 | 最大4時間 | 157 | — |

試行は開発用で、b=400、元configと同じ引数・出力要求、4 GiB address-space cap。最大2 CLIの試行が重なった。
GNU timeのpeakと、timeout前の `/proc/PID/status` のVmHWMサンプルを区別する。
sampledは打ち切りまでの観測最大であり、完走時のpeakではない。

| 入力 | nt | program | 状態 | CLI経過秒 | RSS GiB |
|---|---:|---|---|---:|---|
| bpRNA_CRW_55322 | 4381 | LinearCParty | timeout | 180.00 | — |
| bpRNA_CRW_55322 | 4381 | LinearCParty | timeout | 600.19 | 1.626 (sampled) |
| bpRNA_CRW_55322 | 4381 | LinearPKRaccess | success | 595.17 | 1.598 (GNU) |
| bpRNA_RFAM_43122 | 3995 | LinearPKRaccess | timeout | 900.10 | 1.721 (sampled) |

メモリ失敗は観測していない。4,381 ntのRaccessは4 GiB以内で完走したが、打ち切ったBPPと3,995 nt入力の完走時のメモリは未確認。
較正は旧runの694 JSONL / 25,884成功CLIの `elapsed_seconds` を同じinput/mode/programへ対応させた。
完走試行の新旧比を転用し、少なくとも旧時間の0.5倍を残す。
各jobに0.5秒を加え、全仕事量を1.5倍する。
各幅の仕事量W、最長job Lに対し `W/8 + 7L/8`（8 job以下ならL）を足し、仮定した並列効率0.8で割り、120秒/taskを加える。
8 workerのthroughputを測ったものではなく、新しいrunの実時間や全入力の4 GiB適合を保証する模型ではない。

[較正と全旧JSONLのhash](cpu16/calibration.json)、[4 scriptsと旧run不変性](cpu16/summary.json)、[投入helper](cpu16/submit_e1_long_cpu16.sh)。
