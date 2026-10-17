# GenomeSPOT Lens — ATB_SAMD00000344_longest20_stress_test

Exploratory model-input sensitivity audit. Predictions are not measured physiology; subset effects are not causal gene effects.

Full proteome: 2344 proteins. Selected subset: 468 proteins.

| Target | Full input | Selected input | Change | Full / selected flags |
|---|---:|---:|---:|---|
| oxygen | 0.9983 | 0.9993 | 0.001047 | none / none |
| temperature_optimum | 24.9 | 21.51 | -3.386 | none / none |
| ph_optimum | 5.286 | 5.61 | 0.3245 | none / none |
| salinity_optimum | 1.967 | 1.83 | -0.1363 | none / none |

Oxygen numbers are P(oxygen tolerant), not oxygen concentration. Temperature is °C. Other continuous targets retain upstream units.

| Target | Control | Finite / attempted | Unflagged | Control change 2.5–97.5% |
|---|---|---:|---:|---|
| oxygen | count_matched | 100/100 | 100 | -0.002143 to 0.0008887 |
| oxygen | length_matched | 100/100 | 100 | 0.0009184 to 0.001128 |
| temperature_optimum | count_matched | 100/100 | 79 | -7.377 to 11.71 |
| temperature_optimum | length_matched | 100/100 | 99 | -6.56 to -1.432 |
| ph_optimum | count_matched | 100/100 | 91 | -0.9876 to 0.8561 |
| ph_optimum | length_matched | 100/100 | 100 | -0.1195 to 0.5604 |
| salinity_optimum | count_matched | 100/100 | 82 | -1.967 to 2.704 |
| salinity_optimum | length_matched | 100/100 | 100 | -0.6203 to 0.2455 |

The control envelopes summarize finite predictions, including flagged ones. Missing predictions are counted and excluded from numerical summaries; consult raw records before interpreting results. Envelopes describe random subset variability, not confidence intervals for physiology. Matching uses protein count and fixed length bins; it does not control all biological confounders. Both controls draw from the whole proteome, including the selected group.

DNA is held constant: this is an input-selection audit, not a simulated genome deletion. A full-input prediction is a comparison baseline, not ground truth. Upstream errors for continuous traits are local cross-validation RMSE estimates, not per-genome confidence intervals. Novelty flags indicate unusual model inputs, not scientific discovery.

See provenance.json for checksums and software versions; predictions.json for every subset, protein ID, prediction and warning; summary.tsv for numeric comparisons.
