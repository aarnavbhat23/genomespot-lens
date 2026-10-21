# Historical archaeome reconstruction

Count-validated reconstruction from UniProt release 2023_03 and the published
[archaeome study](https://doi.org/10.1038/s41564-025-02061-0). The original author
input checksum is not published; this is not an author-supplied accession list.

- `audit.json`: archive/catalog checksums, all four input/screen count checks,
  all-positive mapping, 285 published parent-attribution checks and 20 genus checks.
- `protein_registry.csv`: 19,710 canonical/isoform records, source taxon IDs,
  sequence hashes and archived RefSeq/proteome cross-references.
- `source_taxa.csv`: all 233 taxa, reviewed-protein coverage, exact distinct
  canonical 8–50 aa screened-fragment exposure and distinct catalog-positive counts.
- `positive_parent_mappings.csv`: all 14,237 exact occurrences for the 12,623
  published positive sequences. `catalog_row` is **zero-based** in the checksum-
  pinned Nature Data 1 CSV. Coordinates are **one-based inclusive** in the parent.
  Sequence SHA-256 provides another join key. Multiple valid sources are retained.

Taxon counts describe the historical **reviewed-protein screen**. They are not a
uniform screen of complete proteomes. Shared fragments/positives belong to each
source taxon, so summed taxon counts need not equal a globally deduplicated count.
The 128 taxa with zero catalog hits are retained. Median reviewed coverage is
four records per taxon; only five taxa have at least 500 reviewed records.

Reproduction: [study instructions](../../../study/README.md#reconstruct-the-historical-archaeome-input).
Raw source proteins and the 1.60-GB archive stay local under `runs/`.

Source attribution: UniProt Consortium, release 2023_03, **CC BY 4.0**; retain its
[license and archive provenance](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_03/knowledgebase/RELEASE.metalink)
and the linked archaeome paper/data attribution when reusing these tables.
The repository's MIT software license does not replace upstream data licenses.
