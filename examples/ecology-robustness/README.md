# Corrected real ecology comparison — 21 October 2026

This extends the real [Du et al. source data](https://doi.org/10.1186/s40168-025-02145-3)
(CC BY 4.0) retained in `../lactobacillaceae-ecology/`. Code is separately MIT.
The bare `s__` label in 20 genomes across six genera is unresolved species, not
conspecificity. Original notebook 05 incorrectly pooled it; those frozen fits
remain historical, and this analysis supersedes their species-control inference.

- `registry-with-taxonomy-status.csv`: all 1,964 screened genomes, original flags,
  explicit resolved status and **resolved_joint_model_eligible**. Original
  `count_model_eligible` is a legacy flag; do not use it without resolved status.
- `unresolved-species-exclusions.csv`: all 20 excluded assignments, retained in
  descriptive totals. No genus-level guessing or fabricated ANI clusters.
- `corrected-joint-effects.csv`: four comparisons after taxonomy correction.
  The common informative cohort is 1,292 genomes / 35 resolved species.
- `direct-effects.csv`: all 12 pair/endpoint/quality specifications; withheld and
  insufficient-support rows remain, without invented effects.
- `cohort-members.csv`: exact genome IDs for every specification, including
  cohorts that do not pass the ten-species inference gate.
- `species-influence.csv` and `influence-summary.csv`: all 35 joint omissions.
- `pair-influence/`: every species omission in all six direct count specifications,
  their support tables, summaries and provenance. These are post-lead checks.
- `provenance.json`: original input/code/output hashes, protocol commits and versions.

Quality restriction is completeness >=90%, contamination <5%; species overlap is
recomputed after filtering. Count ratios and presence odds ratios are separate
endpoints with different eligibility. All intervals are exploratory species-cluster
t intervals. Holm correction uses a three-contrast family for each endpoint/quality
combination; it is not a declaration of confirmatory discovery. One binomial fit
is withheld for extreme/nonconverged coefficients and two for insufficient support.

Species labels and genomic GC do not supply a full strain phylogeny or proteome
composition. Total screened small-ORF counts and source-study controls are absent.
These data cannot establish measured AMP potency, adaptation, peptide release,
or an ecological effect independent of **all** ancestry/composition/coverage factors.
See notebook 06, the correction notice and the updated research brief.
