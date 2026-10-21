# Research brief for Vinay — 21 October 2026

## What we can now report

**The crude food-associated AMP-count advantage is not robust to within-species
controls. A different human/animal count contrast remains a promising exploratory
lead after genome-quality restriction.** Both statements concern predicted
Macrel small-ORF candidates, not measured activity or archaeome APEX fragments.

A taxonomy audit found that our earlier comparison pooled 20 genomes carrying
only the unresolved label `s__`, across six genera. We corrected that interpretation,
kept the original records/fits visible, and refit after excluding unresolved
assignments from species-based inference. The corrected common cohort is
1,292 genomes / 35 resolved species. Its food/human ratio is **1.02 (95% interval
0.75–1.39)**; the earlier qualitative conclusion remains.

## The surviving exploratory lead

Directly compare human- and animal-associated genomes only within shared resolved
species, require completeness >=90% and contamination <5%, then adjust species,
log genome size, completeness, contamination and GC:

- **379 genomes / 16 species:** human/animal candidate-count ratio **0.56**,
  95% species-cluster t interval **0.40–0.77**; Holm p≈0.0055 within its specified
  three-contrast family.
- Omitting each species gives point estimates **0.53–0.63**; 15 of 16 omission
  intervals remain below 1. One influential-species omission crosses 1.
- Binary presence is imprecise: OR **0.44 (0.15–1.28)**, using a different cohort.
  Insufficient-support and unstable fits are explicitly withheld, not called null.

This survives the tested **species and genome-covariate controls**, not yet all
ancestry/composition/coverage controls. Full strain phylogeny, proteome amino-acid
composition, clonality, source-study effects and exact ORF screening denominators
remain unresolved. Ecological sorting among species can also be genuine biology;
within-species adjustment asks a narrower question. No functional or causal
adaptation claim is justified.

[Current executed notebook 06](../notebooks/06_real_ecology_robustness.ipynb) —
[all estimates, correction and limitations](ecology-comparison/robustness-results-2026-10-21.md).
The archaeome temperature results and prior/root checks remain unchanged and
uncertain; see [notebook 04](../notebooks/04_real_temperature_amp_study.ipynb).

## Original ATB export: stronger search, still unresolved

The deeper search covered the current OSF index (3,009 entries), 15 top-level
components, official GitHub trees/history/issues/releases, 11 public lab GitLab
projects and two alternate APEX branches, ten ATB Hugging Face repositories,
Mendeley, DataCite, Zenodo and all eight bioRxiv versions. The original fragment /
parent / genome / score export was **not located in the checked public routes**.
These scoped checks do not prove no private or unindexed deposit exists.

We also checked the screen's stated scope using one independently selected
released Bakta annotation. Its protein set contains **33,330,305 distinct 8–50-aa
fragments**, versus **3,919,096** reported for the screen. The archived input MD5
and annotation/FASTA/code hashes are retained. **Its original screen membership
is unknown**: this is a conditional scope question, not proof the paper's total
is wrong or a reconstruction of its original screen. Global deduplication could
not reduce the full screen below a member proteome's own union, if included.

The original source restriction/fragmentation procedure and the conflicting
any-target-versus-median MIC threshold require clarification. Genome annotations
and AMR calls alone cannot recover those missing screen assignments or zeros.

[Audited search and counting evidence](atb-export-search-2026-10-21.md) —
[prepared author request, not sent](atb-export-request-2026-10-21.md).
The request specifies the screened registry/status, exact parent/genome links,
denominators, eleven-target scores and intended selection rule. It is ready for
a forward through Vinay's author connection. No outreach was sent automatically.

## Direction and deliverables

Pursue **which ecological AMP associations survive ancestry, composition and
screening ascertainment controls**, with the human/animal count signal as a
candidate for independent verification. A uniformly screened, closely matched
strain cohort and a real strain tree are the decisive next evidence. This is not
a completed mechanism paper or a claim of highest novelty. The ATB resistance
aim remains unfinished; notebook 03's inferential demonstrations are synthetic.

Public project: https://github.com/aarnavbhat23/genomespot-lens

The new correction, exclusions, model statuses, all contrasts, influence fits,
source inventories and data request are included. Notebook 05 remains labeled
historical/superseded; 04 and 06 are the current actual-data studies.


Verification: 59 local tests pass; the current real-data notebooks execute and
refit the specified checks. Source/data/code hashes and separate upstream data
licenses are retained. Nothing has been sent to the mentor or authors.
