# Source search and availability record

Checked 20 October 2026. This is a scoped, reproducible search of primary papers,
their deposits and code. “Not located” records what these routes yielded; it is
not a claim that no author has the files or that they do not exist elsewhere.
Checksums and actual table inventories are in `examples/amp-study/source-audit.json`.

## Author/publication route

The author's [publication page](https://cesarscience.wordpress.com/welcome/publications/)
links to [Google Scholar profile Wyggy_QAAAAJ](https://scholar.google.com/citations?user=Wyggy_QAAAAJ&hl=en).
Direct Scholar requests returned HTTP 403. I followed the author's publication
links, DOI/PMC versions, data-availability statements, supplemental assets,
Mendeley records and the paper-linked lab code instead of assuming profile
access. GitHub searches for archaeasins, AllTheBacteria and APEX were also checked.

## Archaeome study

- [Published article](https://doi.org/10.1038/s41564-025-02061-0).
- Nature-hosted Supplementary Data 1: downloaded CSV has **12,623 distinct
  positive peptide sequences** and eleven target MIC predictions plus a mean.
  The inspected table has no source-genome or parent-protein columns.
- Nature-hosted Supplementary Data 2: **265 selected sequences** with source and
  parent-protein information. Some rows name multiple sources. This is not the
  full positive catalog's mapping table or an inventory of all original inputs.
- [Mendeley version 3](https://doi.org/10.17632/d8yzgtdrcp.3): enumerated all six
  deposited files. Data S1's filename ends in `.csv`, but the downloaded content
  is FASTA with 12,623 records and generic EP identifiers. It does not restore
  parent/genome mapping. Data S2 has the same SHA-256 as Nature Data 2. Other
  deposited files concern GO summaries, statistics and synthesis certificates.
  Some published labels say 12,632; the inspected Nature CSV and Mendeley FASTA
  each actually contain 12,623 records. This report uses the observed file counts.
- Inspected the Nature source workbook for Fig. 1. Its sheets were `Fig 1b`,
  `Fig 1c`, `Fig 1d`, `Fig 1e`; these describe peptide feature/embedding data and
  comparisons, not a complete genome-screen registry.
- Inspected the supplementary PDF. Table 1 provides active-EP and NCBI gene
  totals for twenty **selected genera**. That table is retained as descriptive
  context, with provenance. Whole-genus NCBI gene totals do not establish the
  screened reviewed-protein set or per-genome fragment denominators.
- The paper's code link is [APEX pathogen on GitLab](https://gitlab.com/machine-biology-group-public/apex-pathogen),
  rather than an inferred personal GitHub repository. Inspected its file tree:
  inference code, model files and example sequences are present; a complete
  archaeome source-genome screening manifest was not found in that tree.

**Still required for the primary richness analysis:** original screened protein
accessions/releases; exact source taxa/strains/genome mappings; per-genome
screening counts after the original filters; all qualifying hit mappings,
including the ability to identify completed zero-hit screens. Global totals
alone do not supply these records. A current UniProt reanalysis is possible
but would define a new cohort and must be labeled accordingly.

## AllTheBacteria study

- [Version 8 preprint](https://doi.org/10.1101/2024.03.08.584059).
- Retrieved its complete supplemental archive via [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13370402/supplementaryFiles)
  after the direct PMC download returned HTML and the bioRxiv route rate-limited.
  Inspected the 22-page PDF. Table 1 lists 24 synthesized peptide IDs/sequences;
  it does not provide the full genome-level screen registry needed here.
- Inspected [AllTheBacteria/atb-paper](https://github.com/AllTheBacteria/atb-paper)
  at `98181f8429bddbb1bf98d3a70c1ce8651633f869`, including a file-content search
  for peptide/APEX/encrypted terms. A full linked AMP screening export was not
  located in that snapshot.
- Checked the [official documentation](https://allthebacteria.org/docs/),
  OSF project categories and the [master ATB file index](https://osf.io/download/r6gcp/).
  The downloaded index contained **3,009 file entries**. A search across its
  fields for `peptide`, `apex` and `encrypted` returned no matching entries.
  That filename/index search does not inspect every member of every archive.
- Checked the [Hugging Face ATB dataset](https://huggingface.co/datasets/AllTheBacteria/ATB):
  public genome-annotation sequence tables are available. Those annotations alone
  do not establish which fragments the AMP study screened or its score mappings.
- Located [current AMRFinderPlus output and status files](https://allthebacteria.org/docs/amr/)
  for tool v4.2.5/database 2025-12.03.1. Read the aggregate result header by streaming
  a small prefix. The complete results table is approximately 1.25 GB compressed;
  it was not downloaded wholesale. A version-specific, streaming subset adapter
  is provided in `scripts/prepare_amr_subset.py`.
- Located the [ATB bacterial marker tree and taxonomy](https://allthebacteria.org/docs/phylogeny/).
  The documented tree release and quality filtering must be matched to the cohort.

**Still required for the primary ATB models:** the screened source-genome manifest,
parent-protein/peptide mappings, complete qualifying score calls and per-genome
screening exposure. Genomes, AMR annotations and trees have accessible routes;
the complete AMP-to-genome link was not established by the checked releases.

## Carbon-A

Checked [HuggingFaceBio/Carbon-A-1.2B](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B)
at revision `07f7d86c96f27c3769220205beb0916e67f23448`. It predicts eukaryotic
protein-coding regions; it is not an AMP predictor. Bacterial/archaeal gene calling
is outside the documented training scope. Its model and annotation requirements
are separate from the lightweight GenomeSPOT analysis environment. No weights or
model-supplied remote code were executed. A useful later application would be a
separately validated eukaryotic gene-annotation comparison, not replacement of
the prokaryotic source annotations in these notebooks.

## What is already executed

GenomeSPOT has been run on current, exact-taxon RefSeq genome/proteome pairs for
six organisms represented in the selected archaeome source table. All sixty
trait predictions, warnings, source URLs and input/model hashes are archived in
`examples/amp-study/source-genomes/`. These are current-reference physiology
outputs, not a recreation of the original reviewed-protein screen. They have
not been joined to invented per-genome AMP counts. These six taxa were selected
from the most represented unambiguous taxa in the selected hit table solely
for technical validation; that outcome-derived subset is not an inference cohort.

All three notebooks have a fully executable, explicitly simulated demonstration
and a real-data mode that fails on absent/incompatible input instead of silently
substituting data. The demonstrations test models and joins; no AMR/ecology/
evolution association has yet been established on the complete original cohorts.
