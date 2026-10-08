# Experiment stages

| Stage | Question | Main outputs |
|---|---|---|
| E0 | Do exact calculations match the incorporated upstream CParty corrections and do derived quantities match independent checks? | Reference comparisons and small-case validation. |
| E1 | How do accuracy and resource use depend on beam width? | Errors in log Z/free energy, BPP, contexts, accessibility and decoding objectives. |
| E2 | How do runtime and peak RSS grow with sequence length and scaffold geometry? | Separate-process measurements on synthetic and biological inputs. |
| E3 | What do the different outputs reveal on the same long RNA? | Aligned BPP, contexts, accessibility, samples and representative structures. |
| E4 | How sensitive are results to the supplied scaffold? | Matched-sequence comparisons over predeclared scaffold choices. |

Prepared: a reproducible 18-input synthetic E1 cohort, a 324-job development
config and BPP/context/accessibility analysis. The matrix has not been executed.
Analysis was checked on the existing 20-job smoke and on an independent
two-structure probability distribution. Read [accuracy analysis](accuracy_analysis.md).

Next: prioritize a short biological E1 pilot and a biological E2 panel. Read the
[bpRNA audit](data_audit_bprna_20261008.md) and
[RNAcentral audit](data_audit_rnacentral_20261008.md) for candidate availability,
provenance and remaining checks. Preserve the synthetic cohort as a control
for scaffold geometry and local caps. Pin a reference version and a verified
implementation build, fix E0 enumeration inputs and the biological E1 cohort,
add sampling-frequency and
common-posterior decoder analysis, and add an external API stage driver. Set the E2
length tiers and scheduler resources after small pilots. Select the E3 RNA and
its biological question before choosing displayed regions. Record all
scaffolds independently; sums across scaffolds are not an unconditional
partition function.
