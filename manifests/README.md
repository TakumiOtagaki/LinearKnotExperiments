# Software and input manifests

Formal experiment locks belong here. No changing implementation commit is
declared a paper baseline during this initial setup.

A formal software lock should record repository URL and commit for
LinearCParty and comparison methods, energy-parameter checksum, compiler,
build flags, binary checksums and verification that the build corresponds to
the pinned source. Each run additionally records the experiment commit and
its actual input and config hashes.

Development run manifests are generated inside `results/runs/RUN/run.json`.

`data-sources/` records acquired dataset snapshots separately from selected
experimental cohorts. The local bpRNA copy and the downloaded RNAcentral
panel have recorded hashes; neither acquisition manifest represents an
algorithm run or a completed pilot selection. Original source files are stored
under ignored `results/local/data-sources/` paths.
