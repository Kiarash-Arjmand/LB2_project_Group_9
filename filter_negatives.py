"""Build negative metadata from the course-specific UniProt search export.

The experimental-location restriction belongs in the UniProt query. Do not
filter its results a second time using substrings of the displayed comments.
This script requires the collection record produced with the corrected files.
Python's standard library is sufficient.
"""

import csv
import hashlib
import json
import re
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
LOCATIONS = {
    "cytosol": "SL-0091", "nucleus": "SL-0191",
    "mitochondrion": "SL-0173", "plastid": "SL-0209",
    "peroxisome": "SL-0204", "cell membrane": "SL-0039",
}
QUERY = (
    "reviewed:true AND taxonomy_id:2759 AND fragment:false "
    "AND length:[40 TO *] AND existence:1 AND NOT ft_signal:* AND ("
    + " OR ".join("cc_scl_term_exp:" + term for term in LOCATIONS.values()) + ")"
)
EXPORT_FIELDS = (
    "accession,reviewed,id,protein_name,gene_names,organism_name,length,"
    "ft_signal,lineage,cc_subcellular_location,ft_transmem,protein_existence,fragment"
)
COLUMNS = ["accession", "organism", "kingdom", "protein_length", "tm_helix_first_90"]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def accession_digest(accessions):
    return hashlib.sha256(("\n".join(sorted(accessions)) + "\n").encode()).hexdigest()


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        return reader.fieldnames, list(reader)


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    records = {}
    accession = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                token = line[1:].split()[0]
                accession = token.split("|")[1] if token.startswith(("sp|", "tr|")) else token
                if accession in records:
                    raise ValueError(f"Duplicate FASTA accession: {accession}")
                records[accession] = [line, ""]
            elif line:
                if accession is None:
                    raise ValueError("Sequence before first FASTA header")
                records[accession][1] += line
    if any(not sequence for _, sequence in records.values()):
        raise ValueError("Empty FASTA sequence")
    return records


def write_fasta(path, records, accessions):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in accessions:
            header, sequence = records[accession]
            handle.write(header + "\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def kingdom(lineage):
    # UniProt includes taxonomic ranks, for example "Metazoa (kingdom)".
    terms = {re.sub(r"\s+\([^)]*\)$", "", term.strip())
             for term in re.split(r"[,;]", lineage)}
    for term, group in [("Metazoa", "Metazoa"), ("Fungi", "Fungi"), ("Viridiplantae", "Plants")]:
        if term in terms:
            return group
    return "Other"


def tm_first_90(annotation):
    for feature in re.split(r"(?=TRANSMEM )", annotation):
        start = re.match(r"TRANSMEM (\d+)\.\.", feature)
        note = re.search(r'/note="([^"]+)"', feature)
        if start and note and 1 <= int(start.group(1)) <= 90:
            description = note.group(1).lower()
            if "helical" in description or "signal-anchor" in description:
                return True
    return False


def validate_candidates(fields, rows, records):
    required = {"Entry", "Reviewed", "Organism", "Length", "Signal peptide",
                "Taxonomic lineage", "Subcellular location [CC]", "Transmembrane",
                "Protein existence", "Fragment"}
    if not required.issubset(fields):
        raise ValueError(f"Missing source columns: {sorted(required - set(fields))}")
    accessions = [row["Entry"] for row in rows]
    if not rows or len(accessions) != len(set(accessions)):
        raise ValueError("Empty dataset or duplicate TSV accession")
    if set(accessions) != set(records):
        raise ValueError("TSV and FASTA accession sets differ")
    for row in rows:
        accession = row["Entry"]
        if (row["Reviewed"] != "reviewed" or int(row["Length"]) < 40
                or "Eukaryota" not in row["Taxonomic lineage"]
                or row["Signal peptide"].strip() or row["Fragment"].strip()
                or row["Protein existence"] != "Evidence at protein level"):
            raise ValueError(f"Common selection criteria failed: {accession}")
        if len(records[accession][1]) != int(row["Length"]):
            raise ValueError(f"Sequence length mismatch: {accession}")
        if not row["Organism"].strip():
            raise ValueError(f"Missing organism: {accession}")
    return accessions


def make_metadata(rows):
    return [dict(zip(COLUMNS, [row["Entry"], row["Organism"],
                 kingdom(row["Taxonomic lineage"]), row["Length"],
                 str(tm_first_90(row["Transmembrane"])).lower()])) for row in rows]


def positive_overlap(root, accessions):
    path = root / "data/processed/positive.tsv"
    if path.exists():
        _, positives = read_tsv(path)
        overlap = set(accessions) & {row["accession"] for row in positives}
        if overlap:
            raise ValueError(f"Positive/negative accession overlap: {sorted(overlap)[:5]}")


def build(root=ROOT):
    raw = root / "data/raw"
    tsv, fasta = raw / "negative_candidates.tsv", raw / "negative_candidates.fasta"
    manifest = json.loads((raw / "negative_collection.json").read_text(encoding="utf-8"))
    if manifest["query"] != QUERY:
        raise ValueError("Collection query does not match the course-specific negative query")
    for name in (tsv.name, fasta.name):
        if sha256(raw / name) != manifest["sha256"][name]:
            raise ValueError(f"Input changed since verified collection: {name}")
    fields, rows = read_tsv(tsv)
    records = read_fasta(fasta)
    accessions = validate_candidates(fields, rows, records)
    if (len(accessions) != manifest["search_result_count"]
            or accession_digest(accessions) != manifest["accession_set_sha256"]):
        raise ValueError("Input accessions differ from the verified UniProt search results")
    positive_overlap(root, accessions)
    metadata = make_metadata(rows)
    processed = root / "data/processed"
    processed.mkdir(parents=True, exist_ok=True)
    names = ["negative.tsv", "negative_filtered.fasta",
             "negative_filtered_annotations.tsv", "negative_accessions.txt"]
    with tempfile.TemporaryDirectory(prefix="negative_prepare_", dir=processed) as temp:
        stage = Path(temp)
        write_tsv(stage / names[0], COLUMNS, metadata)
        write_fasta(stage / names[1], records, accessions)
        shutil.copyfile(tsv, stage / names[2])
        (stage / names[3]).write_text("\n".join(accessions) + "\n", encoding="utf-8")
        _, saved_metadata = read_tsv(stage / names[0])
        saved_records = read_fasta(stage / names[1])
        if (saved_metadata != metadata or set(saved_records) != set(accessions)
                or sha256(stage / names[2]) != sha256(tsv)
                or (stage / names[3]).read_text().splitlines() != accessions):
            raise ValueError("Generated output verification failed")
        if any(saved_records[a][1] != records[a][1] for a in accessions):
            raise ValueError("Generated FASTA sequences differ from raw sequences")
        previous = [processed / name for name in names if (processed / name).exists()]
        if previous:
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            backup = root / "data/archive" / (stamp + "_negative_processed")
            backup.mkdir(parents=True)
            for path in previous:
                shutil.copy2(path, backup / path.name)
        for name in names:
            (stage / name).replace(processed / name)
    print(f"Raw negative candidates: {len(rows):,}")
    print(f"Processed negatives: {len(metadata):,}")
    print("Kingdoms:", dict(Counter(row["kingdom"] for row in metadata)))
    print("TM helix starts in residues 1-90:", dict(Counter(row["tm_helix_first_90"] for row in metadata)))
    print("Validated: unique matching IDs, sequence lengths, and common selection criteria.")
    return metadata


if __name__ == "__main__":
    try:
        build()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Negative preparation stopped: {error}")
