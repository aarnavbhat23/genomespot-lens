# A genome-resolved encrypted-AMP release, and what it shows — 27 October 2026

**A public Zenodo deposit supplies something the AllTheBacteria screen does not:
peptide candidates with assembly-level attribution, plus predicted MICs across a
34-strain panel. Its provenance is unverified, and three counting properties in
it must be handled before it is used: accession-level double counting, heavy
peptide redundancy behind the headline candidate total, and a candidate-selection
threshold that changes the qualifying set 5.8-fold.**

This is a scoped audit of released files, not a reanalysis, not a replication of
any published result, and not a substitute for the
[unresolved ATB screen export](atb-export-search-2026-10-21.md).

## Provenance, stated honestly

| Field | Value |
|---|---|
| Title | *Uncovering encrypted antimicrobial peptides in Bacillus by large-scale genomics and machine learning* |
| DOI | [10.5281/zenodo.22107899](https://doi.org/10.5281/zenodo.22107899) |
| Author (as deposited) | Hu, Di |
| Deposited | 26 August 2026, version v1 |
| License | CC BY 4.0 |

**No associated publication was identified.** A title search returned no
matching paper, preprint or journal record. The title is parallel in form to
Du et al.'s Lactobacillaceae study, but the deposit does not name a parent
article, and no lab affiliation is asserted on the record. Treat it as an
unattributed public deposit until provenance is confirmed. Nothing in this audit
depends on who produced it; everything below is a property of the files.

Files retrieved 27 October 2026, SHA-256:

| File | Bytes | SHA-256 |
|---|---:|---|
| `All_6,532_candidate_AMPs.faa` | 567,963 | `f1a9aa21e44256efd0dda7efdf4bfb5721dfb5fea3f7112aaccb4eb540ea0d31` |
| `Predicted_MICs.csv` | 213,649 | `26a9fa52ba57a1040bb532e550259f2e2e2451e2172d6dd71c6504d7ade7b8d8` |
| `The_682_ACF_representative_sequences.faa` | 59,662 | `db0bad03aaf99bea9f73730ad500e5aeb74b0312eee29eefce96409ad653018b` |

## What it has that the ATB export lacks

FASTA headers carry full attribution, for example
`GCA_018741545.1_ASM1874154v1_genomic_CP036094.1_4935` — assembly accession,
assembly name, replicon accession, protein index. All 6,532 headers parse under
that pattern with zero failures, giving 6,532 distinct
assembly/replicon/protein-index triples.

That is a genome-to-parent-protein-to-peptide mapping for positives, which is
one of the five records the ATB resistance aim requires. It still does not supply
the other four: no screened-genome registry, so zero-candidate and failed genomes
are unidentifiable; no per-genome screening exposure; no fragmentation or
deduplication procedure; and no stated selection rule. **Positives alone cannot
establish genome-level zeros**, so this release does not by itself enable a
richness analysis, for the same reason the ATB positives list could not.

## Three counting properties to handle first

### Accession-level double counting

| Quantity | Count |
|---|---:|
| Distinct accession strings | 2,753 |
| Distinct numeric accession cores | 1,427 |
| Cores present as both `GCA_` and `GCF_` | 1,326 |
| Cores with identical peptide sets across the pair | 1,315 of 1,326 |
| Records after core+sequence deduplication | 3,302 of 6,532 |

GenBank and RefSeq accessions for the same assembly both appear, and in 1,315 of
1,326 pairs they carry identical peptide sets. Any per-genome count built on the
accession strings as given would roughly double-count 93% of the genomes.

### The headline candidate total is an occurrence count

| Quantity | Count |
|---|---:|
| Header records | 6,532 |
| Distinct peptide sequences | 682 |
| Peptides occurring once | 40 |
| Occurrences per peptide | median 2, mean 9.6, max 981 |

The 6,532 "candidate AMPs" are 682 distinct peptides recurring across genomes.
One sequence, `KKAVIVENKGCATCSIGAACLVDGPIPDFEIAGATGLFGLWG`, accounts for 981
records across 511 accession cores. The distinct-peptide count matches the 682
in the representative-sequences filename, so the two totals are consistent once
the difference between occurrences and distinct sequences is read correctly.

This bears directly on the ecology and evolution aims. Where candidate richness
is concentrated in a few widely conserved peptides, a cross-genome richness
comparison is largely measuring presence or absence of a small number of
conserved gene families, which is close to measuring phylogeny. The
Lactobacillaceae registry behaves differently — its top family holds 2.8% of hit
records, the top ten hold 19.9%, and 419 of 599 families appear in a single
genome — so richness there is not family-dominated. The contrast means
family-concentration has to be reported per dataset rather than assumed.

### The selection threshold moves the candidate set 5.8-fold

`Predicted_MICs.csv` holds 583 rows against a 34-strain panel. Applying
candidate-selection rules to those released predictions:

| Rule | Qualifying rows | Share |
|---|---:|---:|
| Any target < 80 µM | 433 | 74.3% |
| Median across panel ≤ 80 µM | 75 | 12.9% |
| Any target < 40 µM | 127 | 21.8% |
| All targets ≤ 80 µM | 0 | 0.0% |

Any-target and median-panel rules differ by **5.77×** on the same predictions.
This is the same Methods-versus-Results ambiguity already recorded for the ATB
screen, where Results describe any target below 80 µM and Methods describe a
median-panel MIC at or below 80 µM. It is now quantified on released numbers
rather than inferred: the choice of rule is not a detail, it is most of the
candidate set. Neither 433 nor 75 equals the deposit's own 6,532 or 682, so the
released MIC table is not the file that produced those totals.

Two coverage gaps in the same table: 8 of 583 MIC sequences are absent from the
candidate FASTA, and 107 of 682 candidate peptides have no MIC row. The panel has
34 strains, where APEX 1.1 as described in the ATB preprint predicts 11
pathogens, so the model version behind these predictions is not established.

## What this is and is not good for

Usable now, with the deduplication above applied:

- A worked example of the genome/parent/peptide mapping format to request for ATB.
- A potency endpoint for peptide-level questions, since predicted MICs are
  per-peptide across a panel — relevant to the resistance aim's "stronger"
  half, though not its "richness" half.
- A second clade in which to report family-concentration, as a check on whether
  richness comparisons are measuring conserved families or independent coding.

Not usable:

- No genome-level richness or ecology analysis, because zero-candidate genomes
  are unidentifiable and screening exposure is unrecorded.
- No potency comparison across genomes, because the positives are a selected set
  under an unstated rule.
- Nothing requiring established provenance, until the parent study is identified.

## Reproduction

```bash
# Download to an empty directory; files are untrusted input.
mkdir -p runs/bacillus-release && cd runs/bacillus-release
for f in "All_6,532_candidate_AMPs.faa" "Predicted_MICs.csv" \
         "The_682_ACF_representative_sequences.faa"; do
  curl -fsSLO "https://zenodo.org/records/22107899/files/$(python -c \
    'import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))' "$f")?download=1"
done
cd ../..
python scripts/audit_bacillus_release.py    runs/bacillus-release/All_6,532_candidate_AMPs.faa \
  runs/bacillus-release/Predicted_MICs.csv runs/bacillus-release/The_682_ACF_representative_sequences.faa \
  runs/bacillus-release/parsed.csv
python scripts/audit_bacillus_redundancy.py runs/bacillus-release/parsed.csv \
  runs/bacillus-release/Predicted_MICs.csv
```

Deposited data are CC BY 4.0 and retain that license; this repository's MIT
license covers only the audit code.
