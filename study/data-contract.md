# Input contract for the three AMP notebooks

Version 1. All files use UTF-8 and comma delimiters unless specified otherwise.
One analysis bundle represents one domain, screen protocol, score model/version,
target panel and threshold rule. Freeze releases and source URLs/checksums before
analysis. Never force organism-name-only joins or silently mix ATB sample IDs,
RefSeq accessions and UniProt organisms.

`bundle.json` contains:

```json
{
  "schema_version": 1,
  "data_kind": "real",
  "threshold_mode": "any",
  "threshold_comparison": "lt",
  "threshold_um": 80,
  "targets": ["EXACT_TARGET_ID_1", "EXACT_TARGET_ID_2"],
  "tree_rooted": true,
  "tree_rooting_method": "document the actual method here",
  "sources": {"genomes": "versioned source URL and processing record"}
}
```

The threshold is illustrative too: ATB Results say any-target <80 µM, whereas
Methods say median-panel ≤80 µM. Resolve that discrepancy before preparing a real
ATB bundle. A verified median-panel protocol would require extending and testing
the current aggregator, which currently supports mean and any only.

The target IDs and rooting method above illustrate field names; they must be
replaced with the actual panel and provenance. `threshold_mode` is `any` (minimum
MIC across the panel) or `mean` (arithmetic mean across the panel).
`threshold_comparison` is `lt` (<) or `le` (≤). MIC units are µmol/L. Every peptide
must have all declared panel scores. The archived demonstration has
`data_kind: demo`; real mode rejects that bundle.

## Required tables

| File | Grain and required fields |
|---|---|
| `genomes.csv` | One versioned assembly per row. `genome_id,cohort,assembly_accession,domain,phylum,family,genus,species,protein_count,proteome_aa,gc_fraction,completeness,contamination,habitat`. |
| `screening.csv` | One row per genome, **including zero-hit and failed screens**. `genome_id,status,hits_complete,n_screened,screen_protocol,exposure_unit,score_model,score_version,target_panel_id`. |
| `hits.csv` | One score per genome/peptide/target, allowing repeated parent copies. `genome_id,peptide_id,parent_protein_id,target,predicted_mic_um,mapping_status,peptide_cluster_id`. |
| `genomespot.csv` | Long-format full-proteome predictions: `genome_id,target,value,error,is_novel,warning`. All four optimum/tolerance targets are required for the primary models. Preserve all ten raw predictions separately. |
| `amr_status.csv` | `genome_id,status`; status is `PASS`, `FAIL` or `NOT DONE`. |
| `amr_hits.csv` | `genome_id,element_type,gene_symbol,drug_class,locus_id`. Type is canonical `AMR`, `STRESS` or `VIRULENCE`; empty/unknown types are errors. |
| `phylogeny.nwk` | Rooted, branch-length Newick tree. Tip names match `genome_id` exactly. Rooting method and provenance must be documented. |

`cohort` is an explicit dataset selection such as `atb` or `archaeome`. Use ATB
sample accessions as ATB genome IDs so they join to AMR calls, and keep assembly
versions in their own field. Harmonize GTDB/NCBI historical names using a recorded
mapping. Current source-taxon assemblies are not automatically the genomes whose
proteins the original paper screened.

Genome completeness and contamination are percentages (0–100), GC is a fraction
(0–1), and `proteome_aa` is total amino-acid residues in the matched full proteome.
`habitat` is independent source metadata. Use a recorded `unknown` category where
appropriate; do not invent habitat from the GenomeSPOT predictions. Record the
missingness rate and repeat known-habitat-only analyses.

`screening.status` is `COMPLETE`, `FAIL` or `NOT DONE`. `hits_complete=true` attests
that **all qualifying hits**, not merely a shortlist for synthesis, are supplied.
The exposure unit is `unique_peptide_sequences`: `n_screened` is the number of
distinct sequences actually scored per genome after the same eligibility filters.
Do not substitute whole-genome gene counts, protein counts, the global number of
peptides screened, or an assumed sliding-window count. If original per-genome
exposure cannot be reconstructed, that screen is not eligible for the primary
richness model. A separate protein-normalized descriptive endpoint may be useful,
but must not be labeled the same estimand.

`peptide_id` is a stable sequence identifier shared across genomes, ideally the
SHA-256 of the peptide sequence. Gene copies and multiple target scores are
collapsed to one qualifying sequence per genome. Preserve source coordinates and
all parent mappings outside the model table for audit. Mapping status is `exact`
only when source attribution is established; unresolved genome mappings are
excluded and counted. A peptide occurring exactly in proteins from several
verified sources belongs to each source and is not forced to just one. The
reconstruction preserves these multiple-source occurrences. Missing peptide rows never establish a zero without a completed registry.

Cluster IDs must be assigned consistently across the full candidate cohort using
a recorded sequence-clustering method, threshold and coverage rule. Missing IDs
are allowed and make cluster diversity unavailable. Identity to a known catalog
and distinctness among candidates are different endpoints. Do not claim novelty
from a cluster ID. No peptide design or optimization is performed by the notebooks.

Potency here means conditional **predicted** MIC among qualifying candidates.
If experimental results are available, store them separately in long format with
`measured_mic_um,target,censoring,assay_id,replicate` and assay conditions. Right-
or left-censored MICs require an appropriate censored model; they are not exact
measurements at the reporting limit. The current notebooks do not infer measured
potency or silently replace predictions with assay values.

For ecology/evolution cohorts without AMR data, supply empty AMR tables with the
required headers. Missing status remains missing and does not block those models.
The resistance notebook requires an ATB cohort with actual successful AMR calls.

## Importing observed ATB AMR output

The checked v4.2.5 aggregate result header uses `Name`, `Protein id`, `Contig id`,
`Start`, `Stop`, `Strand`, `Element symbol`, `Type` and `Class`; do not apply a
v3-specific rename blindly. The version-specific adapter streams the large table:

```bash
.venv/bin/python scripts/prepare_amr_subset.py \
  --results AMRFP_results.tsv.gz --status AMRFP_status.tsv.gz \
  --genomes data/amp-study/genomes.csv --out runs/amr-import
```

Copy the two normalized CSVs to the real bundle and retain the import report.
`--core-only` performs a documented `Scope=core` sensitivity. The primary burden
counts distinct reported symbols, including reported resistance-associated
mutation symbols; it is not gene copy number or phenotypic resistance. Missing
AMR status always stays missing. Database/tool versions and any filters must
match across the cohort. See [ATB AMR documentation](https://allthebacteria.org/docs/amr/)
and [AMRFinderPlus interpretation](https://github.com/ncbi/amr/wiki/Interpreting-results).
