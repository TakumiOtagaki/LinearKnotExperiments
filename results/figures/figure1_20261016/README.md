Draft Figure 1 (2026-10-16). Script: `scripts/make_figure1.py`.
Panel C data: exact outputs of build/linear-cparty-febd2bc1 (outputs identical to the adopted versions) for bpRNA_PDB_345 with G = nested pairs of its bpRNA reference structure (from results/local/pk-recovery-bprna-20261016/pk_sel.json).
Commands (CPARTY_FIXED_ENERGY_PARAM_FILE = params/rna_Turner04.par, one thread):
LinearCParty --exact --bpp-cutoff 0 --bpp bpp.tsv; LinearPKCapR --exact; LinearPKRaccess --exact --access-len 1,4,8 --probability;
LinearPKSampling --exact -k 1000 --seed 20261016 --unique; LinearPKMEA --exact --gamma 1; LinearPKSubopt --exact -k 5.
