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

### Historical reconstruction recovered in this search

The Methods give a UniProt access date of **24 August 2023**, with reviewed
canonical and isoform sequences for taxon 2157. That date falls between UniProt
2023_03 (28 June) and 2023_04 (13 September). The
[official archived release statistics](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_03/knowledgebase/UniProtKB_SwissProt-relstat.html)
report 19,708 reviewed archaeal entries. The
[release archive and published checksum](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_03/knowledgebase/RELEASE.metalink)
make this a reproducible route, rather than relying on a live UniProt query.

Downloaded the Swiss-Prot-only archive (1,601,187,312 bytes), verified its
published MD5 `11e3f6e03ad6ee60219974e2ddb0cd54`, and streamed the text and isoform
members. All four reported input/screen counts are reproduced exactly:

| Check | Paper | Reconstruction |
|---|---:|---:|
| Canonical plus isoform protein records | 19,710 | 19,708 canonical + 2 isoforms |
| Distinct protein sequences | 18,677 | 18,677 |
| Source taxa | 233 | 233 |
| Distinct canonical 8–50 aa fragments | 193,331,608 | 193,331,608 |

Exact substring mapping recovers **all 12,623 published positives** in these
proteins, yielding 14,237 parent occurrences. In addition, all 285 parent
accession mentions in the 265-row selected source table are present in this
archive, contain the stated peptide and match a reported taxon. This independently
checks the published source attribution. Summing distinct per-taxon positives
also reproduces all twenty published genus totals, with an explicit
Lokiarchaeum/Candidatus Lokiarchaeum label crosswalk. Those outcome-selected
genera remain validation context only. No AMP model was rerun for these joins.

The committed [reconstruction audit and tables](../examples/amp-study/archaeome-reconstruction/)
include the full protein registry, all positive-to-parent mappings and all
233 source taxa with distinct-fragment exposure. There are 105 taxa with catalog
hits and 128 with no catalog hits under the published mean-panel ≤100 µM rule.
Shared positives occur in multiple taxa (658 sequences); all exact occurrences
are retained. Per-taxon exposures overlap, so their sum is greater than the
global distinct-fragment count. Coordinates are one-based inclusive; catalog
row numbers are zero-based and checksum pinned.

This is a **count-validated historical reconstruction**, not an author-supplied
original accession manifest: the paper supplies no original input checksum for
an identity comparison. Reviewed-source coverage is extremely uneven: one to
1,787 protein records per taxon, with a median of four. Only five taxa have
at least 500 reviewed protein records. Taxon catalog counts are
not whole-proteome discovery rates. Exposure normalization does not eliminate
selection of which proteins were reviewed. A current UniProt query now returns
19,908 reviewed entries in release 2026_03, so silently substituting today's
proteins would not reproduce this cohort.

**Still required for the primary genome models:** resolve the reconstructed taxa
and source accessions to exact versioned whole-genome assemblies, full proteomes,
quality metrics, independent habitats and an appropriate rooted archaeal tree.
A current-reference resolution of all 233 recovered taxa finds 138 unique
exact-taxon complete RefSeq assemblies, leaves 94 unresolved and records one
bounded query failure for the broad uncultured-archaeon taxon. These statuses
and provenance are retained in `examples/amp-study/reconstructed-genomes/`.
Preserve unresolved mappings and duplicate/one-to-many strain assignments.
The original protein/taxon screen and its positive linkage are now available;
the whole-genome model bundle is not yet complete. Reproduction commands are in
[study setup](README.md#reconstruct-the-historical-archaeome-input).

## AllTheBacteria study

- [Version 8 preprint](https://doi.org/10.1101/2024.03.08.584059) and its
  [full PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13370402/).
  **Unresolved protocol discrepancy:** Results describe 1,867 candidates selected
  with predicted MIC <80 µM against at least one target. Methods instead describe
  a median-panel MIC ≤80 µM rule. The 3,919,096-fragment screen therefore needs its
  original export/provenance before choosing a primary real-data endpoint. The
  demonstration uses the Results-style rule solely to exercise the software;
  the implemented aggregator supports mean/any, not median.
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

The initial GenomeSPOT pilot used six exact-taxon RefSeq genome/proteome pairs
from the selected source table, solely for technical validation. Those older
results remain in `examples/amp-study/source-genomes/`.

After the reconstruction, all 233 recovered taxa were checked independently of
AMP outcomes. Downloaded all 138 uniquely resolved current complete RefSeq genome
and proteome pairs and ran GenomeSPOT successfully: **1,380 actual trait
predictions**. Every source input/model hash, feature set, prediction and warning
is retained in `examples/amp-study/reconstructed-genomes/`.

Historical source proteins are all contained in 52 current proteomes, and all
match as full proteins in 47. Current CheckM QC passes 103 genomes; 89 have all
four primary traits unflagged. These intersect in **24 preliminary historical
join candidates**. Other records are retained with their limitations. Source
identity, quality and trait reliability are separate checks; exact taxon labels
alone did not establish protein identity. The 24 candidates are a feasibility
result, not a powered, finalized association cohort. No ecological or resistance
association has been fitted on this reconstructed/current-reference cohort.

All three notebooks have a fully executable, explicitly simulated demonstration
and a real-data mode that fails on absent/incompatible input instead of silently
substituting data. The demonstrations test models and joins; no AMR/ecology/
evolution association has yet been established on the complete original cohorts.
The historical archaeome reconstruction and source joins above are real data
recovery results, separate from those simulated statistical demonstrations.
