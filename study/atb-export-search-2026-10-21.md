# Original AllTheBacteria AMP export: deeper search, 21 October 2026

**The original screened-genome/parent/fragment-score export was not located in
the checked public routes.** Released genomes, annotation corpora and AMR results
are accessible; they do not identify which inputs generated the AMP study's
reported screen, scores and qualifying calls. We have not inferred absent calls
as zeros or labeled a newly generated screen as the original one.

[23 October follow-up: exact public post and linked source checks](atb-social-source-check-2026-10-23.md)
confirm the linked paper version and record public-reply/timeline access limits.

## What was checked beyond the previous search

- The [bioRxiv version API](https://api.biorxiv.org/details/biorxiv/10.1101/2024.03.08.584059)
  returns eight versions; the latest is v8, 7 July 2026. The current Europe PMC
  full-text XML was retrieved and its Methods and data/code links inspected.
- The [OSF root](https://osf.io/xv7q9/) has 15 public top-level components. Its
  current master index again has 3,009 entries and no peptide/APEX/encrypted
  filename match. The Media file route exposes a Logos folder. This does not
  inspect every archive member.
- Five official GitHub repositories were enumerated. The main resource repository
  was cloned, its current file contents searched and historical paths inspected;
  the paper repository and its historical paths were also checked. Its figure
  material is not the full AMP screen. All 111 main-repository issue/PR titles and
  bodies were checked for the relevant AMP terms, and public release assets were
  examined. A public code-index query produced no match; code indexing is limited.
- [The lab's public GitLab group](https://gitlab.com/machine-biology-group-public)
  has 11 public projects. Recursive default-branch inventories and release/branch
  routes were examined, plus both alternative APEX branch trees: 933 blob entries
  across those inventories. An unrelated anti-EV71 release endpoint returned 403;
  that failure is retained. Public APEXDB.csv has 1,736 database rows, without the
  original genome/parent screen mapping; db.csv has 200 demonstration rows.
  Neither is the reported full ATB screen. Runtime user/job databases were not read.
- All ten [ATB Hugging Face repositories](https://huggingface.co/AllTheBacteria)
  were enumerated with exact revisions, file inventories and README snapshots.
  These describe sequence/annotation corpora and representative IDs; no original
  score/mapping export was identified in the inspected material. The multi-terabyte
  Parquet contents were not scanned wholesale.
- Correct Mendeley `query=AllTheBacteria` search returns zero indexed records;
  `ATB20` also returns zero. DataCite returns 20 AllTheBacteria metadata records
  and Zenodo five records, none identifying this export in the reviewed metadata.
  Broad peptide-title queries are nonspecific and are not evidence of absence.

Inventories, exact revisions, API URLs, response hashes and coverage limits are
public under `examples/atb-screen-audit/`. These observations are scoped negative
search evidence, not proof the authors lack a file or that no deposit exists.
The earlier supplementary PDF still supplies only the synthesis shortlist, not
the full screen registry.

## A quantitative scope issue that must be clarified

The [paper](https://doi.org/10.1101/2024.03.08.584059) describes screening all
8–50-aa substrings of Bakta-annotated proteins, but reports 3,919,096 screened
fragments. We performed an independently selected released-input audit:

| Released input check | Exact result |
|---|---:|
| Selected Bakta sample | SAMN34537047 |
| Protein records / unique sequences | 2,773 / 2,772 |
| Distinct canonical 8–50-aa fragments | 33,330,305 |
| Distinct canonical 8–49-aa fragments | 32,612,241 |
| Single input / reported total (8–50 aa) | 8.50× |

The smallest public annotation archive with positive recorded size was selected by size, not AMP
outcomes, downloaded in full and checked against its official MD5. Its first
annotation member supplies the protein sequences. The suffix-array exact-union
counter is independently tested against explicit window enumeration on small
inputs. Source/download/member/code hashes and a compact recountable FASTA are
retained. No AMP predictor was run.

**This sample's original screen membership is unknown.** If every gene/window
from this sample were included, the global union could not be smaller than its
own fragment union. The result therefore asks for the actual source restriction
or fragmentation procedure; it does not establish that the reported total is
incorrect. The difference also survives the inclusive/exclusive 50-aa boundary.

A separate unresolved discrepancy is the endpoint: Results use any target below
80 µM, whereas Methods describe median-panel MIC <=80 µM. Without the original
export or author clarification, reproducing a current arbitrary cohort cannot
resolve either provenance or the primary selection rule.

## Minimum next evidence

The required data are the original screened-source registry,
fragment/parent/genome links, eleven-target predictions,
fragmentation/normalization/deduplication details and the intended threshold.
A positive list alone cannot establish genome-level zeros or full potency
variation. Genome annotations or AMR calls alone do not fill that gap.
