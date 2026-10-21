# Exploratory robustness extension — 21 October 2026

Written after seeing the first ecology comparison; this is not a confirmatory
preregistration. Preserve all prior data, fits, code hashes and results. The
question is whether the crude habitat ranking persists beyond species sorting,
genome covariates, unequal overlap and influential species. Species sorting
itself can be genuine ecology: conditioning on species asks a narrower question
about differences within lineages and does not refute total habitat differences.

## Fixed analyses before inspecting new fits

1. **Direct pair support.** Fit each of the three habitat pairs separately,
   retaining only species observed in both members of that pair. Use the same
   count endpoint, species fixed effects, log genome bp, completeness,
   contamination and GC as the earlier quality-adjusted model. Condition out
   all-zero species, retain individual zero-hit genomes. Compare these direct
   pairs to the earlier joint three-habitat model, which can borrow indirect
   comparisons through connected species. Require at least ten species clusters
   and a full-rank design; otherwise report insufficient support.
2. **Genome quality restriction.** Repeat direct pairs after completeness >=90%
   and contamination <5%, defining species overlap after filtering. Recompute
   eligibility rather than transferring a previous species label. Keep the same
   covariates and minimum cluster count. Report excluded genomes and species.
   This is a sensitivity, not a guarantee of equal screening opportunity: counts
   of all small ORFs screened per genome are still unavailable.
3. **Presence endpoint.** Repeat direct-pair models for the binary indicator of
   any published Macrel hit using a species fixed-effect binomial logit model.
   Only species with variation in the binary outcome contribute; report the
   distinct cohort and odds ratios, never label these count rate ratios. Also
   apply the quality restriction. Withhold nonconverged/separated fits.
4. **Influential species.** Refit the original joint quality-adjusted count model
   after removing each of its 36 informative species in turn. Save all results
   for all three contrasts, not only favorable omissions. Describe effect range
   and interval stability without treating 36 refits as 36 new discoveries.

Use Poisson pseudo-maximum likelihood or binomial GLM mean models as specified,
species-cluster sandwich covariance with finite-sample correction and Student-t
intervals (clusters-1 df). Center/scale continuous covariates within each fit.
Report all unadjusted intervals and p-values, plus Holm adjustment across the
three contrasts within each endpoint/quality family, as exploratory summaries.
The original 1,964-genome registry supplies screened zeros; missing genomes are
not fabricated. Taxonomy is species membership, not a strain marker phylogeny.
Cross-species phylogenetic dependence, original-source effects and imperfect ORF
ascertainment remain limitations, even if every sensitivity agrees.

Preserve baseline inputs and source/hash checks. Export sample IDs, exclusions,
model diagnostics and source/code hashes. Notebook execution must recompute the
new real fits; it must not silently substitute simulated results. Do not claim
measured antimicrobial function or general ecological absence from prediction
counts. ATB encrypted fragments remain a separate endpoint and require their
original source/score/denominator export.
