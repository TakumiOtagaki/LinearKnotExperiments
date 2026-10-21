# Config・plan・投入 script とポイント

すべて experiments/ からの相対パス。未投入。各表の point は優先度0、初回実行のみ。

| Matrix | config / array config / plan | tasks | CLI jobs |
|---|---|---:|---:|
| e2-rnacentral-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e2-rnacentral-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e2-rnacentral-4f97fd3-20261021.json) / [plan](../../results/runs/e2-rnacentral-4f97fd3-20261021/array_plan.json) | 4 | 84 |
| e2-synthetic-shapes-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e2-synthetic-shapes-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e2-synthetic-shapes-4f97fd3-20261021.json) / [plan](../../results/runs/e2-synthetic-shapes-4f97fd3-20261021/array_plan.json) | 55 | 1872 |
| e2-memory-breakdown-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e2-memory-breakdown-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e2-memory-breakdown-4f97fd3-20261021.json) / [plan](../../results/runs/e2-memory-breakdown-4f97fd3-20261021/array_plan.json) | 4 | 35 |
| e2-parallel-impact-w1-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e2-parallel-impact-w1-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e2-parallel-impact-w1-4f97fd3-20261021.json) / [plan](../../results/runs/e2-parallel-impact-w1-4f97fd3-20261021/array_plan.json) | 1 | 36 |
| e2-parallel-impact-w10-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e2-parallel-impact-w10-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e2-parallel-impact-w10-4f97fd3-20261021.json) / [plan](../../results/runs/e2-parallel-impact-w10-4f97fd3-20261021/array_plan.json) | 1 | 36 |
| e3-neat1-bpp-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e3-neat1-bpp-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e3-neat1-bpp-4f97fd3-20261021.json) / [plan](../../results/runs/e3-neat1-bpp-4f97fd3-20261021/array_plan.json) | 2 | 6 |
| e3-neat1-capr-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e3-neat1-capr-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e3-neat1-capr-4f97fd3-20261021.json) / [plan](../../results/runs/e3-neat1-capr-4f97fd3-20261021/array_plan.json) | 2 | 6 |
| e3-neat1-raccess-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e3-neat1-raccess-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e3-neat1-raccess-4f97fd3-20261021.json) / [plan](../../results/runs/e3-neat1-raccess-4f97fd3-20261021/array_plan.json) | 2 | 6 |
| e3-neat1-raccess-byloop-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e3-neat1-raccess-byloop-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e3-neat1-raccess-byloop-4f97fd3-20261021.json) / [plan](../../results/runs/e3-neat1-raccess-byloop-4f97fd3-20261021/array_plan.json) | 2 | 6 |
| e3-neat1-sampling-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e3-neat1-sampling-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e3-neat1-sampling-4f97fd3-20261021.json) / [plan](../../results/runs/e3-neat1-sampling-4f97fd3-20261021/array_plan.json) | 2 | 3 |
| e3-neat1-centroid-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/e3-neat1-centroid-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/e3-neat1-centroid-4f97fd3-20261021.json) / [plan](../../results/runs/e3-neat1-centroid-4f97fd3-20261021/array_plan.json) | 2 | 6 |
| linearknot-pf-4f97fd3-20261021 | [config](../../configs/linearknot_4f97fd3_20261021/linearknot-pf-4f97fd3-20261021.json) / [array](../../configs/linearknot_4f97fd3_20261021/arrays/linearknot-pf-4f97fd3-20261021.json) / [plan](../../results/runs/linearknot-pf-4f97fd3-20261021/array_plan.json) | 5 | 672 |

CLI jobs は warm-up を含む。NEAT1 は0+1、memory breakdown は0+1、他のE2とPFは1+3。

| Script | tasks | 資源 | workers | AS cap MiB/process | 前回最大task分 | h_rt | 計画point | h_rt上限point |
|---|---:|---|---:|---:|---:|---|---:|---:|
| [e2-rnacentral-4f97fd3-20261021/rna_4_8kb](../../jobs/e2-rnacentral-4f97fd3-20261021/rna_4_8kb.sh) | 1 | cpu_40=1 | 10 | 8192 | 29.95 | 01:05:00 | 0.0853 | 0.1300 |
| [e2-rnacentral-4f97fd3-20261021/rna_8_20kb](../../jobs/e2-rnacentral-4f97fd3-20261021/rna_8_20kb.sh) | 1 | cpu_40=1 | 10 | 8192 | 58.34 | 02:05:00 | 0.1624 | 0.2500 |
| [e2-rnacentral-4f97fd3-20261021/rna_20_30kb](../../jobs/e2-rnacentral-4f97fd3-20261021/rna_20_30kb.sh) | 2 | cpu_40=1 | 10 | 8192 | 60.11 | 02:05:00 | 0.3325 | 0.5000 |
| [e2-synthetic-shapes-4f97fd3-20261021/synthetic_n250](../../jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n250.sh) | 1 | cpu_40=1 | 10 | 8192 | 3.60 | 00:15:00 | 0.0151 | 0.0300 |
| [e2-synthetic-shapes-4f97fd3-20261021/synthetic_n500](../../jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n500.sh) | 1 | cpu_40=1 | 10 | 8192 | 18.19 | 00:45:00 | 0.0545 | 0.0900 |
| [e2-synthetic-shapes-4f97fd3-20261021/synthetic_n1000](../../jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n1000.sh) | 5 | cpu_40=1 | 10 | 8192 | 15.14 | 00:35:00 | 0.2031 | 0.3500 |
| [e2-synthetic-shapes-4f97fd3-20261021/synthetic_n2000](../../jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n2000.sh) | 12 | cpu_40=1 | 10 | 8192 | 21.55 | 00:50:00 | 0.7576 | 1.2000 |
| [e2-synthetic-shapes-4f97fd3-20261021/synthetic_n4000](../../jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n4000.sh) | 36 | cpu_40=1 | 10 | 8192 | 28.42 | 01:05:00 | 2.9493 | 4.6800 |
| [e2-memory-breakdown-4f97fd3-20261021/rna_4_8kb](../../jobs/e2-memory-breakdown-4f97fd3-20261021/rna_4_8kb.sh) | 1 | cpu_40=1 | 10 | 8192 | 3.13 | 00:15:00 | 0.0141 | 0.0300 |
| [e2-memory-breakdown-4f97fd3-20261021/rna_8_20kb](../../jobs/e2-memory-breakdown-4f97fd3-20261021/rna_8_20kb.sh) | 1 | cpu_40=1 | 10 | 8192 | 6.37 | 00:20:00 | 0.0224 | 0.0400 |
| [e2-memory-breakdown-4f97fd3-20261021/rna_20_30kb](../../jobs/e2-memory-breakdown-4f97fd3-20261021/rna_20_30kb.sh) | 2 | cpu_40=1 | 10 | 8192 | 7.91 | 00:20:00 | 0.0516 | 0.0800 |
| [e2-parallel-impact-w1-4f97fd3-20261021/parallel_impact](../../jobs/e2-parallel-impact-w1-4f97fd3-20261021/parallel_impact.sh) | 1 | cpu_40=1 | 1 | 8192 | 75.54 | 02:40:00 | 0.2087 | 0.3200 |
| [e2-parallel-impact-w10-4f97fd3-20261021/parallel_impact](../../jobs/e2-parallel-impact-w10-4f97fd3-20261021/parallel_impact.sh) | 1 | cpu_40=1 | 10 | 8192 | 29.46 | 01:05:00 | 0.0842 | 0.1300 |
| [e3-neat1-bpp-4f97fd3-20261021/neat1_short](../../jobs/e3-neat1-bpp-4f97fd3-20261021/neat1_short.sh) | 1 | cpu_40=1 | 3 | 20480 | 8.55 | 00:25:00 | 0.0285 | 0.0500 |
| [e3-neat1-bpp-4f97fd3-20261021/neat1_long](../../jobs/e3-neat1-bpp-4f97fd3-20261021/neat1_long.sh) | 1 | cpu_40=1 | 3 | 20480 | 67.08 | 02:20:00 | 0.1852 | 0.2800 |
| [e3-neat1-capr-4f97fd3-20261021/neat1_short](../../jobs/e3-neat1-capr-4f97fd3-20261021/neat1_short.sh) | 1 | cpu_40=1 | 3 | 20480 | 5.43 | 00:15:00 | 0.0191 | 0.0300 |
| [e3-neat1-capr-4f97fd3-20261021/neat1_long](../../jobs/e3-neat1-capr-4f97fd3-20261021/neat1_long.sh) | 1 | cpu_40=1 | 3 | 20480 | 41.29 | 01:30:00 | 0.1163 | 0.1800 |
| [e3-neat1-raccess-4f97fd3-20261021/neat1_short](../../jobs/e3-neat1-raccess-4f97fd3-20261021/neat1_short.sh) | 1 | cpu_40=1 | 3 | 20480 | 5.47 | 00:15:00 | 0.0192 | 0.0300 |
| [e3-neat1-raccess-4f97fd3-20261021/neat1_long](../../jobs/e3-neat1-raccess-4f97fd3-20261021/neat1_long.sh) | 1 | cpu_40=1 | 3 | 20480 | 43.42 | 01:35:00 | 0.1222 | 0.1900 |
| [e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_short](../../jobs/e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_short.sh) | 1 | cpu_40=1 | 3 | 20480 | 5.50 | 00:15:00 | 0.0193 | 0.0300 |
| [e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_long](../../jobs/e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_long.sh) | 1 | cpu_40=1 | 3 | 20480 | 42.79 | 01:30:00 | 0.1196 | 0.1800 |
| [e3-neat1-sampling-4f97fd3-20261021/neat1_short](../../jobs/e3-neat1-sampling-4f97fd3-20261021/neat1_short.sh) | 1 | cpu_40=1 | 2 | 40960 | 2.86 | 00:10:00 | 0.0122 | 0.0200 |
| [e3-neat1-sampling-4f97fd3-20261021/neat1_long](../../jobs/e3-neat1-sampling-4f97fd3-20261021/neat1_long.sh) | 1 | cpu_40=1 | 2 | 40960 | 8.16 | 00:25:00 | 0.0276 | 0.0500 |
| [e3-neat1-centroid-4f97fd3-20261021/neat1_short](../../jobs/e3-neat1-centroid-4f97fd3-20261021/neat1_short.sh) | 1 | cpu_40=1 | 3 | 20480 | 10.22 | 00:25:00 | 0.0321 | 0.0500 |
| [e3-neat1-centroid-4f97fd3-20261021/neat1_long](../../jobs/e3-neat1-centroid-4f97fd3-20261021/neat1_long.sh) | 1 | cpu_40=1 | 3 | 20480 | 79.41 | 02:45:00 | 0.2185 | 0.3300 |
| [linearknot-pf-4f97fd3-20261021/n250](../../jobs/linearknot-pf-4f97fd3-20261021/n250.sh) | 1 | cpu_16=1 | 8 | 1024 | 2.30 | 00:10:00 | 0.0045 | 0.0080 |
| [linearknot-pf-4f97fd3-20261021/n500](../../jobs/linearknot-pf-4f97fd3-20261021/n500.sh) | 1 | cpu_16=1 | 8 | 1024 | 19.44 | 00:45:00 | 0.0229 | 0.0360 |
| [linearknot-pf-4f97fd3-20261021/n1000](../../jobs/linearknot-pf-4f97fd3-20261021/n1000.sh) | 1 | cpu_16=1 | 8 | 2048 | 191.14 | 06:30:00 | 0.2076 | 0.3120 |
| [linearknot-pf-4f97fd3-20261021/n2000](../../jobs/linearknot-pf-4f97fd3-20261021/n2000.sh) | 1 | cpu_16=1 | 8 | 2048 | 41.44 | 01:30:00 | 0.0467 | 0.0720 |
| [linearknot-pf-4f97fd3-20261021/n4000](../../jobs/linearknot-pf-4f97fd3-20261021/n4000.sh) | 1 | cpu_16=1 | 4 | 8192 | 220.36 | 07:25:00 | 0.2387 | 0.3560 |

合計：82 task、30 scripts、2768 CLI jobs。計画point 6.3811、h_rt上限point **10.0340**。

計画pointは旧版の割当実測に余裕を加えたもの。新版の実測値ではない。taskごとの点を合計しており、資源タイプの異なるtask-hoursを課金係数なしに足していない。
