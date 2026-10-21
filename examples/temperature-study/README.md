# Real temperature/AMP comparative pilot

Actual source inputs for the protocol in `study/temperature/`. These are not the
simulated notebook demonstrations. The cohort is frozen before model fitting.

- `genomes.csv`: all 138 current resolved genomes, explicit primary eligibility
  and exclusion reasons, verified-subset exposure/counts and composition. Curated
  temperatures are accepted only after assembly-pair and strain checks.
- `proteins.csv`: 15,786 taxon/unique-parent records whose complete historical
  sequence occurs in the current proteome; exposure, catalog-positive count,
  composition and homologous-group assignments. No new AMP scoring is performed.
- `phenotypes.csv`: 82 exact versioned assembly-pair TEMPURA matches, including
  11 uncorroborated culture labels withheld from model temperatures. Two external
  culture aliases were verified with cited evidence in the protocol directory.
- `exclusions.csv`: source mapping/phenotype exclusions. Primary exclusions also
  appear explicitly in `genomes.csv`; noneligible rows are not converted to zeros.
- `phylogeny.nwk`: published GTDB r226 ar53 tree pruned to 97 exact versioned
  current assembly tips. Unary branches were collapsed preserving patristic
  distances; the effective covariance root is the retained MRCA. Earlier
  provenance/protocol prose calling it the unchanged published root was
  imprecise; see the explicit root audit and full-published-root fits below.
- `sources/`: the original compressed marker tree and frozen MMseqs2 v18 cluster
  memberships. Clusters are homologous groups, not guaranteed orthologs.
- `homology-results/`: within-family Poisson sensitivity estimates and two-way
  cluster-robust uncertainty. The two identity cutoffs use different parent sets.
- `posterior-results/`: final genuine phylogenetic NB2 fits, diagnostic summaries,
  checks and a sampling-attempt ledger. Full posterior draws remain local.
- `sensitivity-results/`: exploratory coefficient-prior and r232 release checks.
- `rooting-results/`: two further fits preserving full published r226/r232 roots,
  shared stems and root-to-tip depths. Original fits are retained for audit.
- `provenance.json`: raw source and derived input checksums, protocol and versions.

The primary cohort comprises 29 genomes, 8,308 verified parents, and reported
optimal temperatures 22.5–100.5 C. Candidate counts concern historically reviewed
proteins and predicted activity, not a uniform full-proteome AMP screen, measured
activity or peptide expression. This is a small exploratory pilot; it does not
complete the ATB resistance analysis or establish a mechanism.

## Attribution and data licenses

The code is MIT; upstream data are **not relicensed by MIT**. Retain the source
components' notices and attribution when reusing the collection:

- UniProt 2023_03 / archaeome source-derived protein records: **CC BY 4.0**;
  cite [the archaeome paper/data](https://doi.org/10.1038/s41564-025-02061-0).
- TEMPURA-derived phenotype-containing tables/figures: **CC BY-NC 3.0**, Yu Sato
  and colleagues; [database](https://togodb.org/db/tempura),
  [paper](https://doi.org/10.1264/jsme2.ME20074). Preserve noncommercial restrictions.
- Original/adapted GTDB r226 and r232 trees: **CC BY-SA 4.0**, GTDB;
  [source release](https://data.gtdb.ecogenomic.org/releases/release226/226.0/),
  [license](https://gtdb.ecogenomic.org/licenses). The pruned tree is an adaptation
  and retains ShareAlike. Do not apply the software license to it.

[Reproduction and interpretation](../../study/temperature/README.md).
