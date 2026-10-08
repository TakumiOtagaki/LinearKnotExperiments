# Software and input manifests

Formal experiment locks belong here. No changing implementation commit is
declared a paper baseline during this initial setup.

A formal software lock should record repository URL and commit for
LinearCParty and comparison methods, energy-parameter checksum, compiler,
build flags, binary checksums and verification that the build corresponds to
the pinned source. Each run additionally records the experiment commit and
its actual input and config hashes.

Development run manifests are generated inside `results/runs/RUN/run.json`.
