# Actual published-data comparison: bacterial ecology and AMP candidate counts

Source: Du et al., Microbiome (2025), [article](https://doi.org/10.1186/s40168-025-02145-3),
original Additional File 1. Article and supplements are CC BY 4.0; the original
ZIP is retained unchanged under `sources/` and its SHA-256 is in provenance.json.
Derived data retain attribution and that data license, distinct from MIT code.

This is **not an independent replication of the archaeome temperature result**:
Macrel identifies small-ORF candidates, whereas the archaeome screened windows
inside reviewed parent proteins with APEX. Do not pool endpoints or translate
counts into experimentally measured activity or potency.

`genomes.csv` includes all 1,964 explicitly screened MAGs and their species,
habitat, genome size, completeness, contamination, GC and N50. There are 880
positive and 1,084 zero-hit genomes. Missing positives become zero **only within
that complete screened registry**. All 1,327 candidate records agree across the
supplement's S4 and S5, including identifiers, genome, species and habitat.
`hit-registry.csv` retains these joins without adding or rescoring sequences.

The same 1,312 genomes in 36 informative multi-habitat species are used by all
three matched-cohort count models. There are initially 40 multi-habitat species
(1,323 genomes); four all-zero species (11 genomes) cannot inform conditional
count contrasts and are excluded consistently from these three fits. Individual
zero-hit genomes in informative species remain. `species-habitat-overlap.csv`
shows which habitat comparisons each species can inform. Species fixed effects
do not fully control fine-scale phylogeny, cross-species dependence or source-study effects.

Poisson pseudo-maximum likelihood is used as a mean model, with species-cluster
sandwich covariance, finite sample correction, and t intervals (35 df in matched
fits). Genome size is an estimated covariate, **not an exact screened-ORF offset**:
that denominator is not supplied. All pairwise contrasts and omnibus F tests are
reported. Confidence intervals/p-values are exploratory and not multiplicity
adjusted. `permutation.json` is a separate unadjusted omnibus test assuming
within-species exchangeability; it does not adjust quality or source-study bias.

Read the executed `notebooks/05_real_ecology_species_comparison.ipynb` and
`study/ecology-comparison/results.md` for interpretation and reproduction.
