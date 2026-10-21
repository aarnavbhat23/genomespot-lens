# Temperature adaptation and encrypted-AMP density

Protocol frozen before fitting temperature–AMP associations on 20 October 2026.
This extends the three notebooks with an actual-data study, not a replacement
for the still unresolved ATB resistance analysis. It is an exploratory secondary
analysis of a published prediction catalog. No clinical or causal claim is made.

## Question and plausible contribution

Does experimentally curated optimal growth temperature associate with encrypted-
AMP candidate density after controlling screening exposure, which homologous
parent proteins were screened, sequence composition and shared ancestry?

Habitat/AMP associations and extremophile AMP discovery already exist. The
candidate contribution is their separation from reviewed-protein sampling and
shared sequence-based predictors using an independently curated temperature,
versioned source mapping, a real marker tree and within-homology comparisons.
A literature search can support a gap; it cannot establish “highest novelty.”

## Data and eligibility fixed without association results

Use the reconstructed UniProt 2023_03 archaeome and its complete positive catalog
(mean-panel MIC <=100 micromolar). Do not rescore or optimize peptides. Match
historical proteins individually by exact full sequence to the current complete
proteome of the already resolved exact-taxon RefSeq assembly. Keep one record
per taxon/identical parent sequence; report every exclusion. Recompute distinct
screened fragment exposure and candidate unions for the retained subset. Never
transfer all historical hits onto an incomplete current match.

Current assemblies must pass recorded CheckM completeness >=95%, contamination
<=5%. Use independently curated TEMPURA optimal growth temperatures **only on an
exact versioned assembly or its documented GenBank/RefSeq pair**. Reject conflicting
matched temperatures or incompatible strain attributions; no species-name fuzzy
joins or assembly-version stripping. Retain reported temperature ranges and
source entry IDs. TEMPURA values are curated reported physiology, not newly
measured here; they are not natural-habitat measurements.

The primary genome comparison requires >=20 verified distinct reviewed parent
proteins, a unique TEMPURA match and an exact tip in the GTDB r226 archaeal marker
tree. Sensitivities use >=100 verified parents and the broader verified cohort.
Phylogeny comes from the published rooted ar53 tree, preserving branch lengths;
record the original tip and exact assembly mapping. No invented tree.

## Outcomes and models

1. Genome richness is the union of catalog-positive sequence IDs in verified
   proteins. Exposure is the exact union of canonical 8–50 aa fragments in those
   proteins. Fit a phylogenetic NB2 model with log(exposure/100000) offset, reported
   temperature per 10 C, three PCs of screened parent amino-acid composition,
   log(full-proteome size) and log(reviewed fraction). Scale non-temperature
   covariates. Compare a baseline without composition to the adjusted model on
   **the same genomes**; this is not a causal mediation estimate. Use four chains,
   1000 tuning/1000 retained draws initially, target_accept .95, seed 23. Require
   all-parameter Rhat <1.01, bulk ESS >=400, zero divergences; otherwise withhold
   posterior claims and repair diagnostics, with changes recorded.
2. Gene-level sensitivity counts distinct positive fragments per verified unique
   parent sequence, with that parent's exact screened-fragment offset. Cluster
   all archived proteins using MMseqs2 v18 (30% identity, 80% bidirectional coverage,
   greedy incremental mode); these are homologous groups, not proven orthologs.
   Retain families represented in >=5 eligible genomes. Compare temperature-only
   and composition-adjusted Poisson models with homology-family fixed effects,
   adding log(parent length). Use two-way cluster-robust covariance for source
   genomes and homologous groups, Student-t with min(cluster counts)-1 degrees of
   freedom. Report diagnostics and sensitivity to 50% identity clustering. Do not
   label cluster-robust covariance a phylogenetic mixed model. Families with zero
   total candidates cannot inform a conditional family contrast and are reported.
3. Compare GenomeSPOT temperature with independently curated values and repeat
   the genome model with unflagged GenomeSPOT temperature on its observed subset.
   The external table may overlap GenomeSPOT training strains; this checks label
   agreement, not independent predictive generalization. Check published holdouts
   before claiming an external model validation.

Primary estimand: adjusted phylogenetic rate ratio per 10 C with a 95% posterior
interval. Secondary family contrasts and alternate cutoffs are explicitly
sensitivity/exploratory results, not additional discovery claims. Save all fits,
exclusions and diagnostics. Do not select the reported result by significance.

## Interpretation and stopping rules

Review family/temperature coverage and covariate collinearity. If fewer than 20
eligible genomes or fewer than 10 informative homolog groups remain, report a
feasibility result instead of fitting the corresponding primary model. Sparse
phylogenetic coverage and large intervals limit conclusions. A residual
association does not prove thermoadaptive function; an attenuated association
may reflect real compositional adaptation, model bias, or both. Counts refer to
published *predicted* candidates in reviewed proteins, not the complete cryptic
repertoire, peptide expression, release, measured potency or organismal immunity.
Wet-lab validation remains necessary for mechanistic/therapeutic claims.

ATB replication/AMR analysis requires the unresolved complete source/score
export and threshold clarification. Do not manufacture it from synthesized hits.
