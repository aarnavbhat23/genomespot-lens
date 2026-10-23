# ATB original-export search and conditional screen-scope audit

**This directory does not contain the original AMP screen export.**
The deeper public search on 21 October 2026 did not locate that export in the
checked routes. It contains inventories, source hashes and a released-input
counting check, with its limits explicitly recorded in `search-manifest.json`.

The original OSF index is retained (3,009 entries), alongside inventories of
11 public lab GitLab projects and two alternate APEX branches, ten ATB Hugging
Face repositories and the main GitHub tree. Publication metadata confirm eight
bioRxiv versions, latest v8 dated 7 July 2026. Metadata, filenames and documentation
are not a complete inspection of the terabyte-scale archives or Parquet contents.
Private, restricted, unindexed or generically titled deposits can remain unseen.
The unrelated anti-EV71 releases route returned HTTP 403 and is not called checked.

## Released Bakta input, not original-screen membership

We downloaded the smallest indexed Bakta archive with positive recorded size, independently of AMP
outcomes: `atb.bakta.r0.2.batch.174.tar.xz`, 8,640,348 bytes. Its official index MD5
`93c4c3a54584a463a60a64bee4e16ce1` was verified. The first annotation member is
`atb.bakta.r0.2.batch.174/SAMN34537047.bakta.json`; its SHA-256 is in `audit.json`.
The derived `released-proteins.fasta.gz` has 2,773 coding/small-ORF protein records,
2,772 unique sequences. The source archive/annotation remain local; the compact
FASTA and registry permit exact recounting. ATB source-data attribution/license
are distinct from the analysis code; see the [ATB resource](https://allthebacteria.org/).

Exact canonical fragment unions contain 33,330,305 distinct 8–50-aa fragments or
32,612,241 for 8–49 aa. The former is 8.50 times the paper's reported 3,919,096.
**Original screen membership is unknown.** If a screen includes every gene in
this sample and every specified window, even global deduplication cannot reduce
the entire screen below that member's count. This is a conditional scope question,
not proof the published count is wrong and not an original-screen reconstruction.
An unreported source restriction, different fragmentation procedure or description
correction could resolve it. No AMP scores, per-genome hit counts or zeros are
inferred here, and no substitute APEX screen was run.

[Source-search report and required data](../../study/atb-export-search-2026-10-21.md).

To reproduce: download the URL in `audit.json`, verify the index MD5, extract only
the named JSON member, then run `scripts/audit_atb_screen_scope.py` with the saved
annotation and the source-provenance object from `audit.json` saved as JSON.
Raw archive extraction must preserve the member contents for checksum matching.
