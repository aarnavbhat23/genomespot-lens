# Research brief for Vinay — 20 October 2026

## A concrete result and a defensible direction

**Ecological AMP associations can change substantially when the comparison
accounts for source coverage, sequence composition and species membership.**
We now have two separately labeled, real-data analyses supporting a focused
question: **which physiology/niche associations survive ancestry, composition
and screening ascertainment controls?** This is a useful comparative result and
reproducible resource; it is not a validated mechanism or a finished manuscript.

## 1. Archaeome: physiology and ancestry

The historical source reconstruction reproduces the paper's input/screen counts
and maps every published positive sequence to source proteins. GenomeSPOT ran
on 138 current complete source-taxon genomes (1,380 predictions). The stricter
actual comparative cohort contains 29 genomes with corroborated reported culture
temperatures, QC, exact marker-tree tips and 8,308 verified historical parents.
Counts and distinct-fragment exposure are recomputed for that subset.

The baseline ancestry/size/coverage model gives RR 1.34 per 10°C (95% posterior
interval 1.05–1.69). Adding screened-protein composition gives 1.17 (0.77–1.73).
Composition may be real thermoadaptation as well as a shared predictor input;
this change cannot identify causal mediation or prove a prediction artifact.

Five new fits vary coefficient shrinkage, GTDB release and covariance rooting.
All use identical 29-genome membership, pass the original sampling diagnostic
gates and have zero divergences. Every adjusted interval spans 1. The updated
GTDB r232/full-published-root estimate is 1.16 (0.77–1.72). A root audit found
that earlier pruning retained the selected MRCA instead of the full outgroup
root; the description is corrected and both definitions are now explicitly
fitted, without replacing the original results.

[Executed real-data temperature notebook](../notebooks/04_real_temperature_amp_study.ipynb)
— [results and complete sensitivity table](temperature/results.md).

## 2. Bacterial ecology: an independent methodological comparison

The original supplement of [Du et al.](https://doi.org/10.1186/s40168-025-02145-3)
supplies a complete 1,964-MAG screened registry, with 1,084 zero-hit genomes.
All 1,327 positive records and their species/habitat joins validate exactly.
Its Macrel small-ORF endpoint is different from archaeome APEX windows; the
analyses are not pooled or called independent temperature replication.

| Food-associated / human-associated count ratio | Estimate | 95% species-cluster interval |
|---|---:|---:|
| Full cohort, crude | 1.87 | 1.18–2.94 |
| Shared-species informative cohort, crude | 1.52 | 0.89–2.60 |
| Same cohort + species effects | 1.08 | 0.80–1.45 |
| Same cohort + genome size, quality and GC | 0.98 | 0.73–1.32 |

Only 40 of 137 species span habitats; after conditioning out four all-zero
species, all three matched fits use the same 1,312 genomes / 36 species.
A separate unadjusted within-species permutation test gives p=0.0414, while the
quality-adjusted species-cluster omnibus F test gives p=0.101. We retain both:
residual habitat effects remain possible. The supported finding is fragility of
the crude food advantage, not proof of ecological irrelevance. The screened-ORF
denominator, full strain phylogeny and source-study controls remain missing.

[Executed real-data ecology notebook](../notebooks/05_real_ecology_species_comparison.ipynb)
— [all contrasts, diagnostics and limitations](ecology-comparison/results.md).

## 3. What remains for the requested study

- **Resistance aim:** the full AllTheBacteria original screened-genome/parent/
  fragment-score export, including screened zeros, remains unresolved. Published
  threshold descriptions also conflict. The selected 24-peptide synthesis list
  cannot establish per-genome richness, potency distribution or AMR co-evolution.
  The analysis notebook/AMR adapter exist, but inferential examples in notebook 03
  remain explicitly synthetic. No real ATB resistance result is claimed.
- **Carbon-A:** reviewed as requested. Its documented eukaryotic CDS scope does
  not supply an AMP predictor or justify prokaryotic AMP inference here.
- **Novelty:** general AMP ecology and AMP/resistance surveys already exist,
  including [a 2023 global survey](https://doi.org/10.1021/acs.est.3c01664).
  The candidate contribution is exact source recovery plus a controlled,
  reusable audit of ecological associations. Standard statistical adjustment
  itself is not new. No 'first ever' or 'highest novelty' claim is justified.
- **Decisive next evidence:** a uniformly screened ATB cohort with ancestry,
  AMR and valid screening opportunities, plus independent temperature replication.
  Use the frozen ecology/evolution/AMR specifications, matched comparisons and
  source/threshold checks before inspecting favorable outcomes. New claims about
  antimicrobial function or potency need independent measured activity.

## Deliverables and verification

Public project: [aarnavbhat23/genomespot-lens](https://github.com/aarnavbhat23/genomespot-lens).
Five notebooks execute: 01–03 are explicitly synthetic method demonstrations;
04–05 are actual-data studies. Both real-data notebooks were executed locally,
50 tests pass, and 125 installed packages pass dependency checks. Source archives,
exact joins, exclusions, licenses, fit diagnostics, failed-attempt history and
figure-generating code are retained. Full posterior arrays stay local.

No message has been sent to Vinay or the paper authors automatically.
