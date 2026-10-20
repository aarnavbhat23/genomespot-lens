# Prespecified analysis plan

These notebooks implement the three AMP research questions. The earlier
protein-selection audit remains available as a diagnostic tool. Freeze this
plan and source manifests before a real cohort run; changes after inspecting
results must be marked exploratory.

## Shared design

The unit is a versioned genome, with one independent full-proteome GenomeSPOT
prediction set and a complete screen registry. Use bacterial and archaeal cohorts
separately, maintaining their original screening rules. De-duplicate repeated
assemblies; account for repeated isolates within species. Exclude genomes with
completeness <95%, contamination >5%, incomplete screens or missing exposure.
The primary joint ecology models exclude any flagged/missing primary trait;
report inclusive and trait-specific coverage sensitivities rather than hiding
extremophile exclusions. Sampling must be independent of observed AMP richness.

Choose primary cutoffs before fitting: ATB's reported any-target <80 µM and the
archaeome mean-panel ≤100 µM outcome differ. Screening exposure and panel must be
homogeneous before cross-study comparisons. Publish coverage and exclusions,
including unassigned source proteins, ambiguous strain matches and failed calls.

## 1. Ecology and AMP richness

Primary endpoint: distinct qualifying sequences per genome. Fit NB2 counts with
log(distinct screened sequences /100,000) offset. Include standardized temperature,
pH, salinity and oxygen-tolerance probability, log full-proteome amino-acid count,
GC, independent habitat and family effects. Report species-clustered intervals,
NB overdispersion and observed/predicted zeros. The four trait coefficients form
the primary BH-FDR family. Check collinearity, nonlinear relationships and model
fit before claiming a null or positive association.

Family-held-out cross-validation compares numeric size/GC and size/GC/ecology
models with training-fold scaling. Repeat a stricter threshold, retain finite
flagged traits as sensitivity, and assess alternate transformations with a
recorded rationale. A tree-covariance NB2 mixed model is an explicit phylogenetic
sensitivity, not an unadjusted regression labeled phylogenetic.

## 2. Evolution and AMP richness

Prespecify the tested rank/clades and sample-size coverage threshold. Describe
rates by clade and fit clade contrasts with the same offset and size/GC/ecology
controls. Primary ancestry adjustment uses a rooted, branch-length covariance
matrix in a phylogenetic NB2 mixed model. Report exact tip matching, root choice,
prior sensitivity, posterior predictive checks, R-hat <1.01, effective samples
≥400 and zero divergent transitions. Retain failed diagnostics explicitly.

Clade contrasts and phylogenetic random effects can be weakly identified. Use
balanced subsets and alternate defensible trees/rootings. Test bacteria and
archaea separately; a bacterial marker tree does not provide archaeal ancestry.
The dense implementation is capped at 500 genomes. Larger collections need a
scalable model or a prespecified representative panel. The published top/bottom
twenty-genus table is descriptive context only, because it was outcome selected.

## 3. AMR and AMP opportunity

Primary burden: distinct AMR symbols among successful AMRFinderPlus calls. Keep
stress/virulence elements separate. Distinct drug classes and core-only or
curated acquired-only calls are sensitivities; the latter needs an external
classification. Successful no-hit genomes count as zero; failed/unrun status
stays missing. Control annotation/proteome size and fixed source release.

Endpoints: richness; conditional median log10 predicted MIC among ≥5 qualifying
sequences; expected cluster richness after rarefaction to ten sequences. Compare
across-cohort family-adjusted and within-species richness models. The latter
requires AMR variation within at least ten species. For conditional potency and
diversity, adjust candidate count and genome covariates. Treat within-species
richness, potency and diversity as an endpoint testing family with BH adjustment.
The phylogenetic model adds a relatedness sensitivity for richness.

## Evidence standards

Do not call an association co-evolution, a mechanism or clinical opportunity by
itself. Lineage, mobile elements, clinical sequencing and annotation coverage
can confound all three questions. Both GenomeSPOT and AMP predictors use protein
sequence information; composition can drive a circular association. Add
composition/length/coverage controls, independent annotation releases and measured
physiology/activity where possible before a manuscript claim. Retain assay
censoring and target identity. A whole-proteome prediction is not an environmental
measurement, and genetic AMR calls are not measured susceptibility.

The delivered notebooks execute the complete data preparation, count models,
cross-validation, sensitivity analyses, conditional endpoints and phylogenetic
model construction/prior prediction on an explicitly simulated 384-genome
fixture. Posterior MCMC is opt-in. These outputs validate software behavior;
real biological conclusions require the missing complete screen/source linkage.
