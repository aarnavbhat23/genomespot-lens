# Three AMP research notebooks

The notebooks correspond directly to the ecology, evolution and resistance
questions. They share validated data preparation, and each contains its own
question, model, diagnostics, sensitivity analyses, interpretation and sources.

| Notebook | Primary analysis | Additional checks |
|---|---|---|
| [01 Ecology](../notebooks/01_ecology_amp_richness.ipynb) | Negative-binomial richness model with screened-sequence offset and size/GC/habitat/taxonomic controls | Family-held-out prediction, stricter threshold, flagged-trait sensitivity, phylogenetic mixed-model option |
| [02 Evolution](../notebooks/02_evolution_amp_richness.ipynb) | Clade coverage/contrasts and a branch-length phylogenetic negative-binomial mixed model | Exact tree/sample alignment, covariance checks, prior prediction, optional posterior diagnostics |
| [03 Resistance](../notebooks/03_resistance_amp_opportunity.ipynb) | AMR burden versus richness, conditional predicted potency and equal-effort cluster diversity | Within-species comparison, alternative AMR burden, endpoint FDR, phylogenetic richness sensitivity |

**Saved inferential outputs are simulated demonstrations.** Real published
source audits, a selected published genus table and real GenomeSPOT predictions
for six archaeome-linked source taxa are displayed separately. The complete
original AMP-to-genome screening linkage is not yet established; the
[source search](source-search.md) documents what was actually checked. No
ecology/AMR/co-evolution finding is claimed from the demonstration.

## Install and execute

From the repository root, using Python 3.11 and `uv`:

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e '.[test,plot]' -r study/requirements.txt
.venv/bin/python -m pytest tests/test_amp_study.py -q
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
