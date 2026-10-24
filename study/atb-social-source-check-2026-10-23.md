# AllTheBacteria social-post and source-link check — 23 October 2026

**The exact public post was recovered and its paper link followed. It did not
resolve the missing original APEX screen export.** This is a scoped negative
search result, not evidence that no such file exists.

## Exact post and discussion

The [20 July 2026 post](https://x.com/delafuentelab/status/2079170015056871790)
announces AllTheBacteria and its ATB20 experimental demonstration. Its sole
non-media external URL is [version 8 of the preprint](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v8).
There is no direct AMP export or supplementary dataset URL in the post text.
The complete post was read on the public X page and independently retrieved
through the [FxTwitter public text mirror](https://api.fxtwitter.com/status/2079170015056871790).

Three replies were visible without login. None supplied a screen-export link.
X required login to see all replies; the remaining discussion was not inspected.
The live profile loaded its identity but its post feed returned an error/login
prompt. A search-indexed [public timeline mirror](https://www.sotwe.com/delafuentelab)
showed the same pinned announcement and a partial timeline. No additional original
ATB export link was identified in that visible material. This does not constitute
an exhaustive inspection of the account's timeline. A direct download of the
mirror returned a challenge page, which was not treated as timeline data.

## Following the paper's data links

The [Europe PMC full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13370402/fullTextXML)
was retrieved again. Its availability sections point to the ATB website, OSF,
GitHub, Hugging Face annotation tables, AWS assemblies and the ATB CLI. Its
supplementary-material declarations identify one PDF, consistent with the earlier
supplement audit. No separate original fragment-score export URL was identified
in the availability text or declared supplement links.

The old ReadTheDocs URL in the paper currently returns 404. Current entry points:

- [Resource directory](https://allthebacteria.org/resources/).
- [Bakta annotation downloads](https://allthebacteria.org/docs/annotation/).
- [AMRFinderPlus results and matched status files](https://allthebacteria.org/docs/amr/).
- [Documentation-to-OSF links](https://allthebacteria.org/docs/osf_links/).

These routes provide genuine genome and AMR resources. They do not establish
which genomes and fragments entered the original AMP screen or which genomes
were successfully screened with zero candidates.

## Checks against the previous search

The newly downloaded [OSF index](https://osf.io/download/r6gcp/) is byte-identical
to the 21 October index: SHA-256
`b36aa796435185a8ea606e63ed5031f40b6e327853e3ef4aab8d07802fd26d3f`.
It still has no peptide/APEX/encrypted filename match.

The recursive GitHub trees remain unchanged and untruncated:

| Repository | Tree revision | Entries |
|---|---|---:|
| AllTheBacteria/AllTheBacteria | `d6fc56cc2c20a7e18268806d1ac5ec58dcafc9d1` | 142 |
| AllTheBacteria/atb-paper | `98181f8429bddbb1bf98d3a70c1ce8651633f869` | 50 |

Neither tree identifies a newly released AMP export. This check does not inspect
every archive member, unindexed deposit or restricted resource. The broader
[21 October search and scope audit](atb-export-search-2026-10-21.md) remains applicable.

## Remaining data requirement

The resistance analysis still needs the original successfully screened genome
registry (including zero-hit genomes and failed/not-processed status), complete
fragment-to-parent-to-genome mappings, per-genome screening exposure, APEX scores
and the intended candidate-selection rule. A targeted author clarification or a
specific file URL is the next actionable evidence; the announcement alone does
not supply those records. No outreach was sent.

[Machine-readable check and provenance](../examples/atb-screen-audit/social-source-audit-2026-10-23.json)
records the inspected routes, response hashes and access limitations. Full raw
responses are retained locally under `runs/social-source-audit-2026-10-23/`.
