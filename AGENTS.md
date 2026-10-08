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
