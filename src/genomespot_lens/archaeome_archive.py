"""Recover reviewed archaeal proteins and exact fragment/source joins from UniProt."""
import gzip
import io
import re
import tarfile

import ahocorasick
import numpy as np
from pydivsufsort import divsufsort, kasai

ARCHIVE_URL = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2023_03/knowledgebase/uniprot_sprot-only2023_03.tar.gz"
ARCHIVE_MD5 = "11e3f6e03ad6ee60219974e2ddb0cd54"
CATALOG_SHA256 = "be8f62397e6afd97245fb72ebc08e62a1d5afeec42a185d5950ddf625792447a"
EXPECTED = {"protein_records": 19710, "distinct_proteins": 18677, "source_taxa": 233,
            "distinct_screened_fragments": 193331608, "positive_catalog_rows": 12623}


def read_archaeal_records(archive):
    """Stream only the Swiss-Prot text and isoform FASTA members; never extract paths."""
    records, isoforms = [], []
    with tarfile.open(archive, "r:gz") as tf:
        for member in tf:
            if member.name == "uniprot_sprot.dat.gz":
                block = []
                with io.TextIOWrapper(gzip.GzipFile(fileobj=tf.extractfile(member))) as fh:
                    for line in fh:
                        block.append(line)
                        if not line.startswith("//"):
                            continue
                        if any(x.startswith("OC   Archaea") for x in block):
                            text = "".join(block)
                            accession = re.search(r"^AC   ([^;]+);", text, re.M).group(1)
                            taxid = re.search(r"NCBI_TaxID=(\d+)", text).group(1)
                            sequence_block = text.split("\nSQ   ", 1)[1]
                            sequence = "".join(re.findall(r"^     (.*)", sequence_block, re.M)).replace(" ", "")
                            records.append({"accession": accession, "taxid": taxid,
                                "organism": " ".join(x[5:].strip() for x in block if x.startswith("OS   ")),
                                "sequence": sequence, "isoform": False,
                                "refseq": [x[5:].strip() for x in block if x.startswith("DR   RefSeq;")],
                                "proteomes": [x[5:].strip() for x in block if x.startswith("DR   Proteomes;")]})
                        block = []
            elif member.name == "uniprot_sprot_varsplic.fasta.gz":
                canonical = {r["accession"]: r for r in records}
                if not canonical:
                    raise ValueError("Archive must provide canonical records before isoforms")
                def retain(header, parts):
                    if header is None:
                        return
                    accession = header.split("|")[1]
                    base = accession.split("-")[0]
                    if base in canonical:
                        isoforms.append({**canonical[base], "accession": accession,
                                         "sequence": "".join(parts), "isoform": True})
                header, parts = None, []
                with io.TextIOWrapper(gzip.GzipFile(fileobj=tf.extractfile(member))) as fh:
                    for line in fh:
                        line = line.strip()
                        if line.startswith(">"):
                            retain(header, parts)
                            header, parts = line, []
                        else:
                            parts.append(line)
                    retain(header, parts)
                break
    if not records:
        raise ValueError("No archaeal records in the expected archive members")
    return records + isoforms


def distinct_fragment_count(sequences, minimum=8, maximum=50):
    """Count distinct canonical substrings exactly, without materializing all windows.

    For each suffix, lengths beyond its preceding suffix's LCP are new strings.
    Cap at the next noncanonical residue/protein boundary and the length cutoff.
    pydivsufsort's LCP compares to the NEXT suffix, so shift it one position.
    """
    if not 1 <= minimum <= maximum:
        raise ValueError("Invalid fragment length interval")
    sequences = sorted(set(sequences))
    if not sequences:
        return 0
    joined = "#".join(sequences).encode("ascii")
    if not joined:
        return 0
    chars = np.frombuffer(joined, dtype=np.uint8)
    valid = np.isin(chars, np.frombuffer(b"ACDEFGHIKLMNPQRSTVWY", dtype=np.uint8))
    ix = np.arange(len(joined), dtype=np.int32)
    available = np.minimum.accumulate(np.where(valid, len(joined), ix)[::-1])[::-1] - ix
    suffixes = divsufsort(joined)
    lcp = kasai(joined, suffixes)
    preceding_lcp = np.r_[0, lcp[:-1]]
    return int(np.maximum(0, np.minimum(maximum, available[suffixes]) -
                          np.maximum(minimum - 1, preceding_lcp)).sum())


def map_catalog(records, sequences):
    """Return EVERY exact occurrence, retaining repeated parents and ambiguous taxa.

    Coordinates are one-based inclusive. Catalog rows are zero-based and must be
    used with the checksum-pinned catalog; sequence hashes provide another key.
    """
    if len(sequences) != len(set(sequences)) or any(not x for x in sequences):
        raise ValueError("Catalog sequences must be unique and nonempty")
    automaton = ahocorasick.Automaton()
    for i, sequence in enumerate(sequences):
        automaton.add_word(sequence, i)
    automaton.make_automaton()
    mappings = []
    for record in records:
        for end, i in automaton.iter(record["sequence"]):
            mappings.append({"catalog_row": i, "accession": record["accession"],
                "taxid": record["taxid"], "start": end - len(sequences[i]) + 2, "end": end + 1})
    return mappings
