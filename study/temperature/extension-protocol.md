# Exploratory extension, 20 October 2026

This addendum is written after the original 29-genome pilot. These analyses are
exploratory robustness checks and an independent *methodological comparison*,
not an untouched confirmatory test or a new significance search. Preserve the
original protocol, scripts, cohorts, fits and failed-attempt ledger.

## Archaeome sensitivity

Keep the primary 29 genomes, exact fragment exposure, reported temperature,
composition adjustment, genome size and reviewed coverage fixed. Vary all
non-intercept coefficient prior standard deviations from 1.5 to 0.75 and 3;
retain the intercept Normal(0,1.5), phylogenetic SD HalfNormal(1), and NB shape
Exponential(1). This isolates fixed-effect shrinkage without moving the
intercept or variance priors. Compare the GTDB r226 tree with r232 using only
identical versioned assembly tips present in both releases. Report the common
cohort size and both release results; if the cohort changes refit r226 on that
same subset. Preserve published rooting and path lengths. Do not substitute
species-level or unversioned tips.

Four chains, 4000 retained draws, 3000 tuning steps, target_accept .995, seed 23;
require all-parameter Rhat <1.01, bulk ESS >=400, and zero divergences. Failed
models remain visible and withhold scientific claims pending repair. No
selection by interval exclusion of 1. Export source/code/data hashes and
cohort IDs. Raw posterior arrays stay local; public summary tables suffice for
lightweight notebook execution.

## Independent bacterial ecology comparison

Audit the original supplementary tables of Du et al. (2025),
https://doi.org/10.1186/s40168-025-02145-3. Their Macrel small-ORF endpoint is
biologically and computationally different from APEX windows encrypted in
longer proteins. Do not pool outcomes, call this archaeome-temperature
replication, or infer measured potency from candidate counts.

If complete MAG membership (including no-hit genomes), species assignments,
habitat and per-genome candidates can be joined exactly, reproduce the
published habitat totals and quantify habitat/species overlap. Compare crude
habitat count rates with a species fixed-effect count model. Use all documented
habitats; retain zero outcomes only for genomes explicitly in the screened
registry. Report species that contribute within-species habitat contrasts,
excluded unknown/missing annotations, and genomes outside common support.
Use genome size and assembly quality covariates if the tables contain them;
otherwise explicitly state their absence and avoid a causal habitat claim.
Use robust species-clustered uncertainty (small-sample t correction), and a
within-species habitat-label permutation test (9999 permutations, seed 23)
for the omnibus habitat association. All-zero species contribute to descriptive
totals but not conditional count contrasts. Report every habitat coefficient,
not a favorable pair only. If joins/denominators cannot be validated, report the
specific missing fields rather than fabricate a usable dataset.

The useful contribution sought is a reusable way to distinguish ecological
associations from ancestry, sequence composition, and ascertainment. A null,
attenuated, reversed, or unidentifiable result is reportable. This does not
establish absence of ecology or prove predictor bias.

## Rooting audit added before two further fits

Tree inspection after the initial sensitivity fits found that the original
pruning routine collapses unary root paths and resets the induced root branch
to zero. It preserves patristic distances and published clade orientation, but
its effective covariance root is the retained genomes' MRCA, not the full
published outgroup root. The discarded shared stems are 0.04572 (r226) and
0.14953935 (r232) substitutions/site. This corrects earlier prose describing the
root as unchanged; original data and fits are retained for audit.

Fit two additional composition-adjusted models on the same 29 genomes, one for
each release, retaining a unary published root and its full shared stem. Use the
original priors and the same sampling/diagnostic settings. Check original versus
pruned root-to-tip and root-to-MRCA depths, not only patristic distances. Report
both induced-MRCA and published-root results regardless of direction. This is
an exploratory correction/sensitivity, not a retroactive protocol amendment.
