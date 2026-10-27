# Strengthening the human/animal lead — 27 October 2026

**The lead reproduces exactly and is not an assembly-quality artifact. But
habitat is perfectly aliased with source cohort in this registry, a cohort
effect of comparable character is demonstrable with habitat held constant, and
few-cluster-valid inference puts the p-value near 0.065 rather than 0.002. The
signal is a weaker exploratory lead than the frozen interval implies, and it
cannot be separated from source-study effects inside this dataset.**

These are post-hoc checks added after the
[21 October robustness results](robustness-results-2026-10-21.md), chosen from
that document's own list of untested confounders. They are exploratory, not a
preregistered confirmatory analysis, and they do not establish measured
antimicrobial activity.

## 1. Exact reproduction

| Check | Genomes / species | Rate ratio | 95% cluster interval | p |
|---|---:|---:|---:|---:|
| Frozen specification, recomputed | 379 / 16 | 0.555 | 0.398–0.774 | 0.00182 |
| Published 21 October value | 379 / 16 | 0.555 | 0.398–0.774 | 0.00182 |

The pipeline reproduces to three decimal places, so the checks below act on the
same quantity the earlier document reported.

## 2. Assembly fragmentation does not explain it, and runs the wrong way

`contig_n50` and `assembly_status` are present in the registry but were absent
from the frozen control set. Small-ORF calling is sensitive to contig breaks, so
this was a live confounder.

| Added control | Rate ratio | 95% interval | p |
|---|---:|---:|---:|
| Frozen (size, completeness, contamination, GC) | 0.555 | 0.398–0.774 | 0.00182 |
| + log contig N50 | 0.552 | 0.404–0.755 | 0.00104 |
| + log N50 + quadratic log genome size | 0.550 | 0.400–0.755 | 0.00111 |
| + quadratic log genome size | 0.553 | 0.394–0.776 | 0.00202 |
| + assembly status | 0.555 | 0.398–0.774 | 0.00182 |
| Negative binomial, alpha 0.1 / 0.5 / 1.0 | 0.555 / 0.554 / 0.554 | — | 0.0019–0.0031 |

Pearson dispersion is 0.84, i.e. mild underdispersion, so the Poisson
pseudo-likelihood is not anticonservative here.

The direction of the quality difference matters more than the adjustment.
Within-species, human-habitat genomes are the *better* assemblies:

| Covariate | Human median | Animal median | Within-species mean difference (human − animal) |
|---|---:|---:|---:|
| contig N50 | 25,230 | 21,810 | +7,703 |
| genome size (bp) | 1,855,000 | 1,777,000 | +5,255 |
| completeness | 97.49 | 95.43 | +0.83 |
| contamination | 0.520 | 1.360 | −0.453 |
| GC fraction | 0.440 | 0.390 | −0.0004 |

Better assemblies should yield *more* called small ORFs and therefore more
candidates. Human genomes are better assembled and have fewer candidates, so
assembly quality cannot generate this association; it opposes it. This is a
genuine strengthening of the lead.

## 3. Habitat is perfectly aliased with source cohort

Genome identifiers in the source registry carry a cohort prefix. Cross-tabulated
against habitat over all resolved-species genomes:

| Prefix | animal | fermented food | human |
|---|---:|---:|---:|
| `Animal_` | 199 | 0 | 0 |
| `Chicken_` | 1,044 | 0 | 0 |
| `Rumin_` | 1 | 0 | 0 |
| `MAG_` | 0 | 197 | 0 |
| `Human_` | 0 | 0 | 523 |

No prefix spans more than one habitat; in the 379-genome analysis cohort, 0 of 3
prefixes appear in both arms. Habitat is therefore a deterministic function of
source cohort. No covariate adjustment within this registry can separate
host-associated biology from the sequencing, assembly and annotation batch that
produced those genomes. The earlier document listed source-study effects as
"not fully controlled"; they are in fact not identifiable here.

### A cohort effect exists with habitat held constant

The `animal` habitat is the one arm supplied by two cohorts, which permits a
negative control: compare `Chicken_` against `Animal_` genomes of the same
species, both labelled `animal`.

| Within-`animal` cohort contrast | Genomes / species | Rate ratio | 95% interval | p |
|---|---:|---:|---:|---:|
| Chicken / Animal, published quality | 1,018 / 19 | 1.131 | 0.922–1.386 | 0.221 |
| Chicken / Animal, high quality | 424 / 14 | 1.170 | 1.001–1.368 | 0.0485 |
| Chicken / Animal, high quality + log N50 | 424 / 14 | 1.189 | 1.042–1.357 | 0.0143 |

Zero-hit share is 50.7% for `Chicken_` and 64.4% for `Animal_` within the same
habitat and the same species set. So a cohort difference of 1.17–1.19× in the
same endpoint, in the same direction of magnitude, arises with habitat fixed.
That is smaller than the 1.80× implied by the habitat estimate, so batch effects
alone are unlikely to account for all of it — but they are not excluded, and the
design cannot exclude them.

### The direction does survive splitting the comparator

| High-quality decomposition | Genomes / species | Rate ratio | 95% interval | p |
|---|---:|---:|---:|---:|
| Human vs (Chicken + Animal) — frozen | 379 / 16 | 0.555 | 0.398–0.774 | 0.00182 |
| Human vs Chicken only | 298 / 12 | 0.544 | 0.379–0.781 | 0.00345 |
| Human vs Animal only | 124 / 9 | 0.369 | 0.134–1.019 | 0.0534 |

Both animal cohorts give the same direction, which argues against the result
being one anomalous cohort. The `Animal_`-only contrast falls below the
protocol's ten-cluster gate and is reported here only because the frozen
analysis pooled it; it is not an independent replication.

## 4. Few-cluster-valid uncertainty

The frozen interval uses a species-cluster sandwich with 16 clusters and t
critical values on 15 degrees of freedom. Cluster-robust sandwich inference is
known to under-cover at that cluster count. Two alternatives, both on the frozen
specification, 9,999 replicates, seed 23:

| Inference method | p-value |
|---|---:|
| Species-cluster sandwich t (as published) | 0.00182 |
| Wild cluster bootstrap-t, Rademacher weights, null imposed | 0.0646 |
| Within-species randomization of habitat labels | 0.0076 |

The two valid methods disagree by an order of magnitude because they assume
different things. The wild cluster bootstrap treats species as the resampling
unit, which is the conservative reading when the effect is carried by a few
clusters. The randomization test conditions on species composition and permutes
genome labels, which is anticonservative if genomes within a species-habitat
cell are clonal or otherwise non-exchangeable — a condition this registry cannot
rule out, since no ANI, MLST or strain metadata is supplied.

The defensible summary is that the evidence sits near conventional
significance, not two orders of magnitude inside it.

## 5. How little information the 379 genomes carry

Species-conditional contrasts use only within-species habitat variation. Arm
sizes for the 16 informative species:

| Smallest arm size | Species |
|---:|---:|
| 1 | 9 |
| 2–3 | 2 |
| 5–6 | 3 |
| 11, 19 | 2 |

Nine of 16 species contribute a contrast in which one habitat has a single
genome, including `Ligilactobacillus ruminis` (42 human against 1 animal) and
`Lacticaseibacillus paracasei` (20 human against 1 chicken). The sum of
smallest-arm sizes is 60 of 379 genomes, and the entire cohort contains 290 AMP
hit records. Leave-one-species-out estimates span 0.528–0.632, and the single
omission whose interval reaches 1 is `Ligilactobacillus sp945923665`, matching
the earlier influence analysis.

The effect is also substantially a zero-inflation difference: 63.4% of human
genomes are zero-hit against 45.1% of animal genomes, with means 0.44 and 0.96.
The separately reported presence endpoint, odds ratio 0.44 (0.15–1.28), remains
too imprecise to settle whether this is a prevalence or a per-genome-count
phenomenon.

## 6. What this changes

Retained:

- The association reproduces exactly and is stable to fragmentation, genome
  size, assembly status, dispersion family and single-species omission.
- Assembly quality confounding runs against the observed direction, so the
  simplest technical explanation is ruled out.
- The direction appears against both animal source cohorts.

Downgraded:

- The p-value should be read as roughly 0.065 under cluster-valid resampling,
  not 0.0018.
- Habitat and source cohort are perfectly aliased, and a 1.17–1.19× cohort
  effect is demonstrable within a single habitat. "Human-associated" and
  "produced by the `Human_` cohort" are the same variable in this registry.
- Effective information is far below 379 genomes: 60 genomes of paired
  within-species support and 290 total hits.

This remains a lead worth pursuing, not a finding. The decisive design is
unchanged from the earlier document and now better specified: same-species
strains from human and animal hosts, sequenced, assembled and annotated through
one pipeline, with the screened small-ORF denominator recorded per genome, and
isolate-level strain typing to establish independence.

## Reproduction

```bash
python scripts/audit_ecology_batch_structure.py        examples/lactobacillaceae-ecology/genomes.csv
python scripts/audit_ecology_cohort_aliasing.py        examples/lactobacillaceae-ecology/genomes.csv \
                                                       examples/lactobacillaceae-ecology/hit-registry.csv
python scripts/audit_ecology_cohort_decomposition.py   examples/lactobacillaceae-ecology/genomes.csv
python scripts/audit_ecology_few_cluster_inference.py  examples/lactobacillaceae-ecology/genomes.csv
```

Seed 23 throughout. Source data are Du et al., *Microbiome* (2025),
[article](https://doi.org/10.1186/s40168-025-02145-3), CC BY 4.0; the cohort
prefixes analysed in section 3 are the identifiers as published in that
supplement, and their mapping to specific originating studies has not been
confirmed against those studies' own accessions.
