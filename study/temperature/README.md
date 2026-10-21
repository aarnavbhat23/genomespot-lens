# Temperature adaptation and encrypted-AMP density

This is an **actual-data comparative pilot** following the ecology and evolution
questions, with a protocol committed before model fitting. It does not complete
the unresolved AllTheBacteria resistance study or validate antimicrobial function.

Question: does optimal growth temperature associate with predicted encrypted-AMP
density after screening exposure, source-protein coverage, sequence composition
and ancestry controls? A within-homology comparison provides a separate check on
which parent protein functions were sampled.

- [Actual results](results.md), [real-data notebook](../../notebooks/04_real_temperature_amp_study.ipynb),
  [protocol](protocol.md), [literature comparison](literature.md),
  [source/licensing notes](data-notes.md), [execution/diagnostic record](run-notes.md).
- The derived real data and diagnostic summaries are in
  `examples/temperature-study/`; full posterior draws stay in `runs/`.
- The first three requested notebooks retain their explicitly simulated software
  demonstrations. A separate real-data notebook documents this pilot.

## Reproduce source preparation

Install the project and `study/requirements.txt` using Python 3.11. First run the
historical reconstruction and all-source genome download steps in
[the main study instructions](../README.md). MMseqs2 v18 is an external executable;
record its version and use the following frozen clustering settings:

```bash
mkdir -p runs/temperature-sources
curl -fL 'https://togodb.org/togodb/fetch/1376?rp=10000&page=1' \
  -o runs/temperature-sources/tempura.json
curl -fL https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_r226.tree.gz \
  -o runs/temperature-sources/gtdb-ar53.tree.gz
mmseqs easy-cluster runs/archaeome-reconstruction/reviewed_proteins.fasta \
  runs/temperature-sources/homology30 runs/temperature-sources/tmp30 \
  --min-seq-id 0.3 -c 0.8 --cov-mode 0 --cluster-mode 2 --threads 4 -s 7.5
mmseqs easy-cluster runs/archaeome-reconstruction/reviewed_proteins.fasta \
  runs/temperature-sources/homology50 runs/temperature-sources/tmp50 \
  --min-seq-id 0.5 -c 0.8 --cov-mode 0 --cluster-mode 2 --threads 4 -s 7.5
.venv/bin/python scripts/prepare_temperature_study.py \
  --raw runs/temperature-sources --reconstruction runs/archaeome-reconstruction \
  --out runs/temperature-data
```

The source scripts resolve current NCBI genomes, so later annotation downloads may
change. TEMPURA is live and can also change. Compare raw and derived hashes in the
committed provenance rather than claiming identical historical reproduction from
an unpinned new query. Frozen derived inputs are available in the repository.

## Fit the real models

The repaired four-chain runs use 3,000 tuning and 4,000 retained draws per chain,
with target_accept .995. They take local CPU time; no paid API is used.

```bash
.venv/bin/python scripts/analyze_temperature_study.py \
  --data examples/temperature-study --out runs/temperature-posterior \
  --draws 4000 --tune 3000 --target-accept .995
.venv/bin/python scripts/analyze_temperature_homology.py \
  --data examples/temperature-study --out runs/temperature-homology
```

To repeat the label-agreement and training-overlap audit:

```bash
.venv/bin/python scripts/compare_temperature_labels.py \
  --data examples/temperature-study --upstream vendor/GenomeSPOT
```

The model runner saves every diagnostic result. Any failed model is explicitly
marked `WITHHOLD_DIAGNOSTICS_FAILED`; it must not support a biological claim.
The protocol defines prior/cohort choices and the interpretation limits. Formal
publication additionally requires prior/tree sensitivity, adequate power and
external replication; statistical output alone does not establish a new mechanism.
