# Three AMP research notebooks

The notebooks correspond directly to the ecology, evolution and resistance
questions. They share validated data preparation, and each contains its own
question, model, diagnostics, sensitivity analyses, interpretation and sources.

| Notebook | Primary analysis | Additional checks |
|---|---|---|
| [01 Ecology](../notebooks/01_ecology_amp_richness.ipynb) | Negative-binomial richness model with screened-sequence offset and size/GC/habitat/taxonomic controls | Family-held-out prediction, stricter threshold, flagged-trait sensitivity, phylogenetic mixed-model option |
| [02 Evolution](../notebooks/02_evolution_amp_richness.ipynb) | Clade coverage/contrasts and a branch-length phylogenetic negative-binomial mixed model | Exact tree/sample alignment, covariance checks, prior prediction, optional posterior diagnostics |
| [03 Resistance](../notebooks/03_resistance_amp_opportunity.ipynb) | AMR burden versus richness, conditional predicted potency and equal-effort cluster diversity | Within-species comparison, alternative AMR burden, endpoint FDR, phylogenetic richness sensitivity |

**Notebooks 01–03 use simulated inferential demonstrations.** A separate
[real temperature pilot](temperature/README.md) now uses verified protein subsets,
curated strain temperatures and a real marker tree; its
[results](temperature/results.md) do not yet establish an independent temperature
effect. The original notebooks remain software demonstrations. Real published
source audits, a selected published genus table and real GenomeSPOT predictions
for 138 current source genomes are displayed separately, with the earlier
six-taxon pilot retained. A count-validated
archaeome reconstruction now supplies the 233 source taxa, per-taxon screened
fragment exposure, all 12,623 positive peptide mappings and zero-catalog-hit taxa.
Historical protein identity and whole-genome QC remain explicit; ATB’s original
full screen export remains
unresolved. The [source search](source-search.md) documents the evidence. No
ecology/AMR/co-evolution finding is claimed from the demonstration.

## Install and execute

From the repository root, using Python 3.11 and `uv`:

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e '.[test,plot]' -r study/requirements.txt
.venv/bin/python -m pytest tests/test_amp_study.py tests/test_archaeome_archive.py -q
.venv/bin/python scripts/execute_study_notebooks.py
```

The runner uses the current Python interpreter explicitly and executes all three
notebooks, saving their outputs in place. It performs no network downloads or
paid compute. The delivered notebook metadata uses a generic Python 3 kernel;
select this environment if opening them in Jupyter or an editor. Shared code
lives in `src/genomespot_lens/amp_study.py`; scientific decisions are explained
in the notebook text and [analysis plan](analysis-plan.md).

The fixture has 384 explicitly synthetic genomes, fake peptide identifiers,
fake AMR calls, simulated traits and a simulated tree. Its purpose is testing
joins, offsets, filters and models. `scripts/make_study_demo.py` regenerates it
in a new directory with seed 23. Raw published inputs stay outside Git in `runs/`.

## Reconstruct the historical archaeome input

The checked release archive is 1.60 GB compressed. Download it locally; only the
compact source registry, parent mappings and audit are committed. This command
requires the analysis dependencies installed above:

```bash
mkdir -p runs/archive-inputs
curl -fL --retry 2 -o runs/archive-inputs/uniprot_sprot-only2023_03.tar.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_03/knowledgebase/uniprot_sprot-only2023_03.tar.gz
curl -fL --retry 2 -o runs/archive-inputs/archaeome_scores.csv \
  'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_MOESM3_ESM.csv'
.venv/bin/python scripts/reconstruct_archaeome.py \
  --archive runs/archive-inputs/uniprot_sprot-only2023_03.tar.gz \
  --catalog runs/archive-inputs/archaeome_scores.csv \
  --genus-table examples/amp-study/published_genus_table.csv --out runs/archaeome-reconstruction
```

The script verifies pinned checksums, canonical+isoform counts, distinct protein
and taxon counts, the exact distinct-fragment total and complete positive mapping.
Pass `--selected PATH_TO_MOESM4_ESM.csv` to independently check all 285 parent
accession/sequence/taxon mentions in the published 265-row selected table. The
genus-table option checks all twenty published genus totals, using an explicit
Lokiarchaeum/Candidatus Lokiarchaeum label crosswalk.
It counts distinct fragments exactly through suffix arrays, without allocating
193 million strings. Small exhaustive tests check boundaries, noncanonical
residues and shared fragments. The output FASTA stays local. Tables use taxon IDs;
they must not be passed to the genome models as if they were complete assemblies.
Median reviewed-protein coverage is only four records per taxon, so source
coverage is a major scientific confound even after exposure normalization.

## Resolve and run the reconstructed source taxa

These steps process all 233 recovered taxa independently of AMP outcomes. The
resolver retains ambiguous/unmatched taxa. The download step compares historical
source proteins with each current full proteome, keeping exact and containment
coverage separately. Current-reference technical predictions do not by themselves
establish a valid historical screen join.

```bash
.venv/bin/python scripts/resolve_archaeome_genomes.py \
  --taxa examples/amp-study/archaeome-reconstruction/source_taxa.csv \
  --out runs/archaeome-genome-resolution
.venv/bin/python scripts/download_archaeome_genomes.py \
  --mapping runs/archaeome-genome-resolution/source_genome_mappings.csv \
  --reconstruction runs/archaeome-reconstruction --out runs/archaeome-current-inputs
.venv/bin/python scripts/run_genomespot_manifest.py \
  --manifest runs/archaeome-current-inputs/genome_inputs.csv \
  --models vendor/GenomeSPOT/models --out runs/archaeome-current-predictions
.venv/bin/python scripts/summarize_archaeome_batch.py \
  --mapping runs/archaeome-genome-resolution/source_genome_mappings.csv \
  --inputs runs/archaeome-current-inputs --predictions runs/archaeome-current-predictions \
  --out runs/archaeome-batch-export
```

The checked current-reference resolution finds 138 unique exact-taxon complete
RefSeq assemblies, leaves 94 taxa unresolved, and records one bounded query
failure for the broad “uncultured archaeon” taxon (>1,000 descendants). All 233
statuses are retained in `examples/amp-study/reconstructed-genomes/`. Reference
assembly choice and original protein identity are separate checks. Current NCBI
CheckM values are preserved where present; absent values are not filled in.

All 138 downloads and GenomeSPOT runs completed, yielding 1,380 trait predictions.
All historical source proteins are contained in 52 current full proteomes (47
match as full proteins); 103 genomes pass the quality thresholds and 89 have all
four primary traits unflagged. The intersection contains 24 preliminary join
candidates. These checks are distinct: their counts must not be added together.
The actual results, per-genome features/model provenance and exclusions are in
`examples/amp-study/reconstructed-genomes/`. Twenty-four candidates are a
feasibility audit, not a completed association cohort; final modeling also needs
taxonomy, independent habitat, a defensible tree and an adequately powered design.

## Real-data run

1. Prepare a release-pinned, outcome-independent source cohort and the original
   screen registry. Resolve genome, parent protein and peptide IDs. Follow the
   [data contract](data-contract.md); do not use selected synthesis candidates
   as a complete richness table.
2. Obtain source contigs and full matching proteomes. Run GenomeSPOT via the
   manifest runner below, using the pinned upstream source/models from the
   repository's main installation instructions.
3. Join the same ATB sample IDs to the release-matched AMRFinderPlus result and
   status tables, using `scripts/prepare_amr_subset.py`. Complete only genuinely
   successful no-hit AMR records with zero; other statuses remain missing.
4. Supply a rooted, appropriate marker tree with branch lengths. Document the
   tip mapping and rooting method. The full tree may need to be pruned outside
   the notebooks to a prespecified representative cohort.
5. Point the notebooks at the real bundle:

```bash
AMP_STUDY_MODE=real AMP_STUDY_COHORT=atb \
AMP_STUDY_DATA="$PWD/data/amp-study" \
.venv/bin/python scripts/execute_study_notebooks.py
```

This intentionally errors if real inputs are missing or the supplied bundle
identifies itself as simulated. There is no automatic demonstration fallback.
Study outputs and provenance are written under `runs/amp-study/<analysis>/<mode>/`.

## Full-proteome GenomeSPOT batch

Create `genome_inputs.csv` with `genome_id,contigs,proteins`. File paths are
relative to that manifest. The runner validates sequence IDs, records input/model
checksums and package versions, saves genome features/predictions, resumes only
matching completed jobs and retains failures explicitly:

```bash
.venv/bin/python scripts/run_genomespot_manifest.py \
  --manifest data/genome_inputs.csv --models vendor/GenomeSPOT/models \
  --out runs/source-predictions
```

Copy the successful long-format `genomespot.csv` into the real analysis bundle
and retain the status/provenance files. A failed prediction does not establish
a zero trait value. Confirm that protein and nucleotide inputs are from the
same strain/assembly; sequence parsing alone cannot establish that relationship.

To repeat the already executed, bounded source-taxon pilot:

```bash
.venv/bin/python scripts/source_genome_pilot.py \
  --models vendor/GenomeSPOT/models --out runs/source-genomes --limit 6
```

This resolves current complete RefSeq assemblies for reported source taxa and
is explicitly not a reconstruction of the original reviewed-protein input set.

## Phylogenetic posterior sampling

Every notebook builds and checks its real phylogenetic likelihood and runs a
prior-predictive check. Full four-chain MCMC is opt-in:

```bash
AMP_STUDY_MODE=real AMP_STUDY_COHORT=atb \
AMP_STUDY_DATA="$PWD/data/amp-study" AMP_STUDY_PHYLO_MCMC=1 \
.venv/bin/python scripts/execute_study_notebooks.py
```

The default posterior run uses 1,000 tuning and 1,000 retained draws per chain,
four chains and two local CPU workers. Runtime depends on cohort and hardware;
the dense implementation refuses more than 500 genomes. It checks convergence,
effective sample size and divergences before exporting posterior estimates.
The standard demonstration executes the model likelihood/prior prediction and
frequentist analyses; it does not imply posterior MCMC was completed.

## Carbon-A assessment

The requested [Carbon-A model](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B)
was reviewed and its revision recorded. Its documented task is eukaryotic CDS
annotation, outside the bacterial/archaeal scope of these source cohorts. It
has not been inserted into the prokaryotic pipeline or used as an AMP scorer.
See the source-search note for an appropriate separate evaluation direction.
