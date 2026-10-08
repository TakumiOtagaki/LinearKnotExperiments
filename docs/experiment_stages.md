# Experiment stages

| Stage | Question | Main outputs |
|---|---|---|
| E0 | Do exact calculations match the incorporated upstream CParty corrections and do derived quantities match independent checks? | Reference comparisons and small-case validation. |
| E1 | How do accuracy and resource use depend on beam width? | Errors in log Z/free energy, BPP, contexts, accessibility and decoding objectives. |
| E2 | How do runtime and peak RSS grow with sequence length and scaffold geometry? | Separate-process measurements on synthetic and biological inputs. |
| E3 | What do the different outputs reveal on the same long RNA? | Aligned BPP, contexts, accessibility, samples and representative structures. |
| E4 | How sensitive are results to the supplied scaffold? | Matched-sequence comparisons over predeclared scaffold choices. |

Next: pin a reference version and a verified implementation build, choose the
E0/E1 input set, and add parsers for numerical accuracy summaries. Set the E2
length tiers and scheduler resources after small pilots. Select the E3 RNA and
its biological question before choosing displayed regions. Record all
scaffolds independently; sums across scaffolds are not an unconditional
partition function.
