# Execution and diagnostic changes

- Protocol committed as `bdb6752` before temperature–AMP model fitting.
- All current annotations were checked individually against each historical
  source protein by exact sequence; 15,786 taxon/unique-parent records survived.
- TEMPURA is a current table with 8,645 rows (549 archaeal), rather than silently
  treating its paper's older total as today's inventory. JSONP is parsed as JSON
  data without evaluating JavaScript. Its data license is CC BY-NC 3.0.
- Exact paired assembly matches are also checked against strain aliases.
  Uncorroborated labels are withheld, not “corrected” by guessing. Two missing
  aliases were verified independently in DSMZ/NCBI and recorded with evidence in
  `strain-aliases.json`; this happened before coefficient inspection.
- MMseqs2 v18 ran on all 19,710 historical protein records at prespecified
  30% and 50% identity, 80% bidirectional coverage, mode 2, four threads.
- The first baseline phylogenetic attempt used 1,000 tuning/1,000 retained draws
  per chain. Sampling warned about Rhat >1.01; saving also exposed the missing
  optional HDF5 backend. No biological result was reported from that attempt.
- Added the pinned h5py backend. Diagnostic repair increases all requested fits
  to 2,000 tuning/4,000 retained draws per chain with target_accept .98 and the
  same four chains/seed/model. This changes sampling effort, not the hypothesis,
  filtering or model to obtain a preferred direction/significance.
- The smaller/alternate homology cutoffs are sensitivities. Report both, including
  wide intervals and inconsistent conclusions; do not choose the smaller p-value.

- The 4,000-draw baseline reached adequate Rhat/ESS but retained one divergence.
  Its interval is withheld pending a dedicated baseline repair using 3,000 tuning
  draws and target_accept .995, retaining the same 4,000 draws/four chains/model.
- Initial 4,000-draw adjusted and >=100-parent fits also retained one divergence
  each. All remaining model variants are therefore rerun at the same stricter
  3,000-tuning/target_accept .995 setting. Final exports must explicitly identify
  the repaired run and pass the original diagnostics; earlier failed runs remain
  in the local audit and are not used to report a biological effect.

- Final baseline, adjusted and all three sensitivity models passed at .995:
  zero divergences, all-parameter Rhat <1.01 and bulk ESS >=400. Final exports use
  these stricter runs only; all five earlier summaries remain in the audit ledger.
