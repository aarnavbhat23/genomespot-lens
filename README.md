# GenomeSPOT Lens

**Explore microbial physiology, ancestry, AMR burden and encrypted-AMP richness.**

Three executable research notebooks implement ecology, evolution and resistance
analyses with explicit screening denominators and source-mapping checks:

1. [Ecology and AMP richness](notebooks/01_ecology_amp_richness.ipynb)
2. [Evolution and AMP richness](notebooks/02_evolution_amp_richness.ipynb)
3. [Resistance and AMP opportunity](notebooks/03_resistance_amp_opportunity.ipynb)

A separate [real-data temperature study](notebooks/04_real_temperature_amp_study.ipynb)
now follows the ecology/evolution question using 29 genomes, independently curated
temperatures, a real archaeal marker phylogeny and verified parent proteins.
Its composition-adjusted temperature interval remains uncertain; see the
[actual results and limits](study/temperature/results.md). The protocol was frozen
before fitting. Novelty and mechanism remain unestablished.

Read [study setup and execution](study/README.md) first. The saved statistical
outputs in notebooks **01–03** remain explicitly simulated software demonstrations. Actual GenomeSPOT results for 138 current exact-taxon source genomes
(1,380 trait predictions) and a thorough public-data audit are included.
**New source reconstruction:** the UniProt 2023_03 archive reproduces all four
archaeome input/screen counts, and all 12,623 published positive sequences map
back to its proteins. The [233-taxon registry and full parent mappings](examples/amp-study/archaeome-reconstruction/)
are included with checksums. All 233 source taxa were checked; 138 had unique current complete RefSeq
matches. Protein identity, genome QC and trait checks are retained separately,
with 24 preliminary historical-join candidates. Those original all-trait checks are distinct from the new per-protein
verified temperature study; ATB’s full AMP screen
export has not been located. See
[source search](study/source-search.md). The notebooks do not turn selected hit
lists into zero-hit genomes or claim biological results from the demonstration.

## Protein-selection audit

Does the protein set supplied to a physiology model change its answer?

GenomeSPOT Lens audits how selecting a protein group changes predicted oxygen
tolerance and growth temperature. It compares the full input with a selected
subset and two random controls: equal protein count, and equal counts within
protein-length bins. It retains GenomeSPOT's uncertainty estimates, novelty
flags, missing predictions and provenance.

Useful for researchers evaluating curated protein collections, gene-family
panels, and core/accessory proteome selections before interpreting microbial
physiology. You supply the group membership; Lens does not annotate gene
functions or construct a pangenome.

**Status: working research prototype with tested notebooks and genome pilots.**
Biological hypotheses still require a complete real cohort. GenomeSPOT supplies the predictive
models; Lens supplies the selection audit and matched controls. A fuller
literature review is needed before claiming methodological novelty.

## What is included

- Strict FASTA and protein-ID validation; optional exact-sequence mapping.
- A deterministic, seeded audit with both random-control schemes.
- Oxygen probabilities consistently expressed as P(oxygen tolerant).
- JSON, TSV and readable Markdown reports, with per-run warnings and IDs.
- Input/model SHA-256 checksums and software versions.
- Unit tests and upstream parity tests, plus GitHub Actions CI.
- A bounded, reproducible AllTheBacteria/archaeal pilot.
- [Research protocol](docs/research-protocol.md), [paper guide](docs/paper-guide.md),
  [novelty assessment](docs/novelty.md) and [pilot findings](docs/pilot-results.md).

## Install

Use Python 3.11. The upstream serialized models require scikit-learn 1.2.2;
changing this or upgrading NumPy without checking compatibility can break inference.
Commands below assume `uv` is installed and that you are in this repository.

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e '.[test,pilot]'
git clone https://github.com/cultivarium/GenomeSPOT.git vendor/GenomeSPOT
git -C vendor/GenomeSPOT checkout f79f7d6fdece71b766ff47f4fb7114daae5287d7
uv pip install --python .venv/bin/python --no-deps -e vendor/GenomeSPOT
```

The source checkout must remain unmodified. Models are read from it rather than
redistributed here. For the exact pilot package set, install
`requirements-pilot-lock.txt` before the editable installs.

## Run your own audit

Provide nucleotide contigs and their corresponding full protein FASTA from the
**same assembly and annotation**, plus a text file with one selected protein ID
per line. FASTA IDs are the first whitespace-delimited token. Do not supply short
peptide hits as a whole-organism proteome. GenomeSPOT does not perform gene calling.

```bash
.venv/bin/python -m genomespot_lens.cli audit \
  --contigs genome.fna.gz --proteins proteins.faa.gz \
  --subset-ids selected_ids.txt --models vendor/GenomeSPOT/models \
  --out runs/my-audit --replicates 100 --seed 23
```

The default report focuses on oxygen and temperature. For the later ecology
analysis, add `--targets oxygen,temperature_optimum,ph_optimum,salinity_optimum`.
All ten upstream predictions are preserved in raw JSON regardless of report targets.

If your subset has different sequence IDs:

```bash
.venv/bin/python -m genomespot_lens.cli map-subset \
  --proteins proteins.faa.gz --subset-fasta curated.faa \
  --out runs/mapping
```

Use `runs/mapping/subset_ids.txt` in the audit and inspect `mapping.json`.
Ambiguous, unmatched and partial sequences are never silently assigned.
This checks sequence correspondence, not whether the nucleotide/protein files
belong to the same assembly; you must establish that through accession metadata.

Use the `python -m` invocation shown here, including for pytest. The upstream
package installs a `genome_spot.py` script beside Python executables; executable
launchers can otherwise resolve that file instead of the importable package.

## Reproduce the real-data pilot

```bash
.venv/bin/python scripts/pilot.py --models vendor/GenomeSPOT/models \
  --out runs/pilot --replicates 100
GENOMESPOT_ROOT="$PWD/vendor/GenomeSPOT" .venv/bin/python -m pytest -q
```

This downloads small inputs directly from official repositories, calls proteins
for one ATB assembly with Pyrodigal, and tests a current reviewed-protein subset
for *Methanocaldococcus jannaschii*. It performs 404 model evaluations: two
baselines, two selected inputs and 400 controls. Output directories must be new.
Public databases change; checksums capture the retrieved inputs, and later runs
may differ. It does not reconstruct the archaeome paper's historical input release.

To regenerate the archived figure, install `.[plot]` and run
`.venv/bin/python scripts/plot_pilot.py`.

## Interpretation

These are model-input perturbations, **not gene knockouts or causal genetic
effects**. DNA remains unchanged. Whole-proteome predictions are comparison
baselines, not ground truth. Control envelopes describe sampling variation,
not measurement confidence intervals. Length matching uses fixed 100-aa bins;
inspect residue totals and repeat with other bin widths before research claims.

An upstream novelty flag means unusual input relative to the training data,
not a novel organism or discovery. Predicted oxygen tolerance is not habitat
oxygen concentration. Predicted pH preference is not measured soil pH.

## Credit and license

Lens code is MIT licensed. GenomeSPOT is by Barnum and colleagues at Cultivarium
and collaborators: [paper](https://doi.org/10.1101/2024.03.22.586313),
[upstream repository](https://github.com/cultivarium/GenomeSPOT).
See [third-party notice](THIRD_PARTY_NOTICES.md). Data retain their source terms;
this repository's software license does not relicense third-party datasets.
