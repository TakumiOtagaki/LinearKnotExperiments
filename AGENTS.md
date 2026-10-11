# LinearCPartyExperiments

This repository owns paper experiments: input preparation, configuration,
execution, aggregation and figure generation. It must work without the paper
repository. Keep software paths configurable and leave external implementation
working trees unchanged.

- Read README.md and docs/experiment_protocol.md before running experiments.
- Use development runs while the implementation is changing. Freeze source and
  build provenance before collecting manuscript measurements.
- Preserve failures, timeouts, commands, raw output, input hashes and versions.
- Keep large raw runs out of ordinary Git history. Commit small configurations,
  inputs, scripts and reviewed summaries.
- Do not replace maintained implementation tests with paper-level smoke checks.
- The shared /gs/bs project has a file (inode) limit used by the whole lab (`lfs quota -p 102571 /gs/bs`). Check it before large submissions and keep the planner's inode check enabled. Never write per-job files to shared storage; run jobs in node-local scratch and store one JSONL and one tar per task. Delete regenerable intermediates and pack retained raw outputs.
