# GenomeSPOT Lens — M_jannaschii_current_reviewed_exact_matches

Exploratory model-input sensitivity audit. Predictions are not measured physiology; subset effects are not causal gene effects.

Full proteome: 1811 proteins. Selected subset: 1177 proteins.

| Target | Full input | Selected input | Change | Full / selected flags |
|---|---:|---:|---:|---|
| oxygen | 9.362e-06 | 1.076e-05 | 1.398e-06 | outside_training_distribution / outside_training_distribution |
| temperature_optimum | 78.66 | 77.75 | -0.9095 | none / none |
| ph_optimum | 5.654 | 5.965 | 0.3117 | none / none |
| salinity_optimum | 1.339 | 2.87 | 1.53 | outside_training_distribution / outside_training_distribution |

Oxygen numbers are P(oxygen tolerant), not oxygen concentration. Temperature is °C. Other continuous targets retain upstream units.

| Target | Control | Finite / attempted | Unflagged | Control change 2.5–97.5% |
|---|---|---:|---:|---|
| oxygen | count_matched | 100/100 | 0 | -3.288e-06 to 4.765e-06 |
| oxygen | length_matched | 100/100 | 2 | -3.004e-06 to 4.564e-06 |
| temperature_optimum | count_matched | 100/100 | 64 | -4.439 to 7.764 |
| temperature_optimum | length_matched | 100/100 | 78 | -4.868 to 5.925 |
| ph_optimum | count_matched | 100/100 | 81 | -1.076 to 0.6307 |
| ph_optimum | length_matched | 100/100 | 88 | -0.7949 to 1.037 |
| salinity_optimum | count_matched | 100/100 | 0 | -1.339 to 2.142 |
| salinity_optimum | length_matched | 100/100 | 0 | -1.339 to 2.325 |

The control envelopes summarize finite predictions, including flagged ones. Missing predictions are counted and excluded from numerical summaries; consult raw records before interpreting results. Envelopes describe random subset variability, not confidence intervals for physiology. Matching uses protein count and fixed length bins; it does not control all biological confounders. Both controls draw from the whole proteome, including the selected group.

DNA is held constant: this is an input-selection audit, not a simulated genome deletion. A full-input prediction is a comparison baseline, not ground truth. Upstream errors for continuous traits are local cross-validation RMSE estimates, not per-genome confidence intervals. Novelty flags indicate unusual model inputs, not scientific discovery.

See provenance.json for checksums and software versions; predictions.json for every subset, protein ID, prediction and warning; summary.tsv for numeric comparisons.
