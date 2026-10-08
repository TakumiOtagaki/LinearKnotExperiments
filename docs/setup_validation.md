# Initial setup validation

On 2026-10-08 the development smoke ran two synthetic inputs, exact and beam
width 100, and five CLI programs: 20 of 20 measured jobs completed successfully.
The runner captured commands, stdout, stderr, the explicit parameter file,
input and executable hashes, elapsed time and per-process peak RSS. The
crossing-scaffold fixture produced crossing pairs in LinearCParty and MEA.

The implementation source state was
`c803e9c3802b1a9a470c541945912096fbfee779` with local changes. Source and
executable hashes did not change during the smoke. The existing binaries have
not been certified against this source revision; this is development validation
of the execution path, not a scientific baseline or a timing benchmark.

The runner integration test passed with one expected failed command, one
timeout and one success, each with a warm-up. The measured summary retains all
three statuses, records the exit code 7, preserves copied inputs and refuses
an existing output directory without changing its manifest.

Reproduce the checks with the commands in README.md. Raw setup runs are local
under ignored `results/runs/`; no manuscript result is inferred from them.
