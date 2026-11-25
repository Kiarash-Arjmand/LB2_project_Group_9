"""Shared inputs and cleavage coordinates for the analysis and von Heijne model."""

from pathlib import Path

import pandas as pd
from Bio import SeqIO

AA = "ACDEFGHIKLMNPQRSTVWY"
# Expasy Swiss-Prot release 2026_04, release statistics section 6.1.
# https://web.expasy.org/docs/relnotes/relstat.html
SWISS_PROT_PERCENT = dict(zip(AA, (
    8.25, 1.38, 5.46, 6.71, 3.86, 7.07, 2.27, 5.90, 5.79, 9.64,
    2.41, 4.06, 4.75, 3.93, 5.52, 6.66, 5.36, 6.85, 1.10, 2.92,
)))


def load_split(name, data_dir=Path("data/splits")):
    """Return split metadata and a sequence map, with exact matching IDs."""
    frame = pd.read_csv(data_dir / f"{name}.tsv", sep="\t", dtype={"accession": str})
    if frame["accession"].duplicated().any() or not set(frame["label"]) <= {0, 1}:
        raise ValueError(f"{name}: duplicate accession or invalid label")
    seqs = {}
    for label, suffix in ((1, "pos"), (0, "neg")):
        records = {}
        for record in SeqIO.parse(data_dir / f"{name}_{suffix}.fasta", "fasta"):
            parts = record.id.split("|")
            accession = parts[1] if len(parts) >= 3 else parts[0]
            if accession in records:
                raise ValueError(f"{name}: duplicate FASTA accession {accession}")
            records[accession] = str(record.seq).upper()
        expected = frame.loc[frame.label == label, ["accession", "protein_length"]]
        if set(records) != set(expected.accession):
            raise ValueError(f"{name}: FASTA IDs do not match metadata for class {label}")
        for row in expected.itertuples(index=False):
            if len(records[row.accession]) != row.protein_length:
                raise ValueError(f"{name}: length mismatch for {row.accession}")
        seqs.update(records)
    return frame, seqs


def cleavage_window(seq, site):
    """Return positions -13 through -1 and +1 through +2 around the cut.

    UniProt's SIGNAL end is the number of residues before cleavage; it is
    therefore also the zero-based boundary between the two sides.
    """
    site = int(site)
    if site < 13 or site + 2 > len(seq):
        raise ValueError(f"Cannot extract a complete 15-residue window at {site}")
    return seq[site - 13:site + 2]
