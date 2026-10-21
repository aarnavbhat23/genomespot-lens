# Source interpretation and reuse

This analysis uses real published sequences and predictions, current versioned
assemblies, independently curated phenotype records and a published marker tree.
No new AMP design, optimization, expression experiment or antimicrobial assay is
performed. GenomeSPOT is run on full current proteomes; historical AMP outcomes
are restricted to individually verified, exactly matching source proteins.

TEMPURA's live table returned 8,645 records (549 archaeal); its paper describes an
older 8,639-record snapshot. Raw response SHA-256 is frozen. JSONP wrapping is
parsed as data without evaluating JavaScript. Assembly matching initially found
82 records, but culture labels must also be corroborated. Eleven assembly-linked
records remain uncorroborated after the documented alias checks. These may be
metadata errors, missing synonyms or different strains; they are not automatically
classified as proven errors. Two additional aliases were verified using DSMZ and
NCBI taxonomy records before coefficient inspection. No fuzzy organism-name join
is used. Only reported strain phenotypes with corroborated aliases enter models.

Pruning the GTDB r226 ar53 tree preserves patristic distances and exact versioned
RefSeq tip IDs. A subsequent audit found that collapsing unary root paths made
the effective covariance root the retained MRCA, rather than the full published
root described in earlier prose. Both releases were additionally fitted with the
full published root and shared stem retained; see extension-protocol.md and the
separate rooting-results provenance. Original inputs and fits remain unchanged. It is a marker phylogeny, not a tree inferred from
AMP counts. A full audit retains genomes missing this release's representative
tree; missing tips are not grafted by species name.

## Licenses and attribution

Software is covered by the repository's MIT license. Upstream data licenses are
separate and must be respected; MIT does not relicense them:

- [UniProt](https://www.uniprot.org/help/license), 2023_03: CC BY 4.0. The original
  archaeome paper/deposit attribution is preserved in the reconstruction audit.
- [TEMPURA](https://togodb.org/db/tempura): CC BY-NC 3.0, as displayed by the source
  download notice. Phenotype-containing tables/figures retain that source notice
  and are provided for attributed noncommercial research reuse.
- [GTDB](https://gtdb.ecogenomic.org/licenses), r226/r232 ar53 trees: CC BY-SA 4.0. The
  pruned tree is an adaptation (tips selected/renamed, unary branches collapsed)
  and remains CC BY-SA 4.0 with GTDB attribution.

Keep source components and their notices separate when reusing this collection.
Raw source downloads and full posterior draws remain local under `runs/`; scripts,
source hashes, derived study tables, diagnostic summaries and figures provide the
reproducible public analysis. Metadata joins do not create experimental evidence.
