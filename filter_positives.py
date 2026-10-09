import csv
from pathlib import Path

# Locate the TSV relative to this script.
project_folder = Path(__file__).resolve().parent
tsv_file = project_folder / "data" / "raw" / "positive_candidates.tsv"

# Read the table. TSV columns are separated by tabs.
with open(tsv_file, encoding="utf-8", newline="") as file:
    reader = csv.DictReader(file, delimiter="\t")
    proteins = list(reader)

print("Number of proteins:", len(proteins))
print("First protein:", proteins[0]["Entry"])
print("Signal peptide:", proteins[0]["Signal peptide"])

import re

kept_proteins = []
excluded_proteins = []

for protein in proteins:
    annotation = protein["Signal peptide"]
    cleavage_site = None

    # Check each signal peptide feature separately.
    features = re.split(r"(?=SIGNAL )", annotation)

    for feature in features:
        match = re.match(r"SIGNAL 1\.\.(\d+)(?:;|$)", feature)

        if match is None:
            continue

        end_position = int(match.group(1))
        experimental = "ECO:0000269" in feature
        not_cleaved = "not cleaved" in feature.lower()

        if experimental and end_position >= 14 and not not_cleaved:
            cleavage_site = end_position
            break

    if cleavage_site is not None:
        protein["cleavage_site"] = cleavage_site
        kept_proteins.append(protein)
    else:
        excluded_proteins.append(protein)

print("Proteins kept:", len(kept_proteins))
print("Proteins excluded:", len(excluded_proteins))

# Create a folder for the filtered output.
output_folder = project_folder / "data" / "processed"
output_folder.mkdir(parents=True, exist_ok=True)

# Read the FASTA and index its records by UniProt accession.
fasta_file = project_folder / "data" / "raw" / "positive_candidates.fasta"
fasta_records = {}

with open(fasta_file, encoding="utf-8") as file:
    for record in file.read().split(">"):
        if not record.strip():
            continue

        header = record.splitlines()[0]
        accession = header.split("|")[1]
        fasta_records[accession] = ">" + record.rstrip() + "\n"

# Check that every selected protein has a sequence.
for protein in kept_proteins:
    if protein["Entry"] not in fasta_records:
        raise ValueError("Missing sequence: " + protein["Entry"])

# Save the annotations, including the extracted cleavage site.
columns = list(reader.fieldnames) + ["cleavage_site"]

for filename, records in [
    ("positive_filtered_annotations.tsv", kept_proteins),
    ("positive_excluded_annotations.tsv", excluded_proteins),
]:
    with open(output_folder / filename, "w",
              encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

# Save the full protein sequences for the selected proteins.
with open(output_folder / "positive_filtered.fasta", "w",
          encoding="utf-8") as file:
    for protein in kept_proteins:
        file.write(fasta_records[protein["Entry"]])

print("Saved filtered TSV and FASTA.")
print("Saved excluded entries for inspection.")
print("Output folder:", output_folder)

import re

kept_proteins = []
excluded_proteins = []

for protein in proteins:
    annotation = protein["Signal peptide"]
    cleavage_site = None

    # Check each signal peptide feature separately.
    features = re.split(r"(?=SIGNAL )", annotation)

    for feature in features:
        match = re.match(r"SIGNAL 1\.\.(\d+)(?:;|$)", feature)

        if match is None:
            continue

        end_position = int(match.group(1))
        experimental = "ECO:0000269" in feature
        not_cleaved = "not cleaved" in feature.lower()

        if experimental and end_position >= 14 and not not_cleaved:
            cleavage_site = end_position
            break

    if cleavage_site is not None:
        protein["cleavage_site"] = cleavage_site
        kept_proteins.append(protein)
    else:
        excluded_proteins.append(protein)

print("Proteins kept:", len(kept_proteins))
print("Proteins excluded:", len(excluded_proteins))

# Read the lineage information.
lineage_file = project_folder / "data" / "raw" / "positive_candidates_lineage.tsv"

with open(lineage_file, encoding="utf-8", newline="") as file:
    lineage_reader = csv.DictReader(file, delimiter="\t")
    lineage_by_id = {
        row["Entry"]: row["Taxonomic lineage"]
        for row in lineage_reader
    }

# Assign the course's kingdom labels.
final_proteins = []

for protein in kept_proteins:
    lineage = lineage_by_id[protein["Entry"]]

    if not lineage.strip():
        raise ValueError("Missing lineage: " + protein["Entry"])

    if "Metazoa" in lineage:
        kingdom = "Metazoa"
    elif "Fungi" in lineage:
        kingdom = "Fungi"
    elif "Viridiplantae" in lineage:
        kingdom = "Plants"
    else:
        kingdom = "Other"

    final_proteins.append({
        "accession": protein["Entry"],
        "organism": protein["Organism"],
        "kingdom": kingdom,
        "protein_length": protein["Length"],
        "cleavage_site": protein["cleavage_site"],
    })

# Save the five-column positive table.
final_columns = [
    "accession", "organism", "kingdom",
    "protein_length", "cleavage_site"
]

with open(output_folder / "positive.tsv", "w",
          encoding="utf-8", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=final_columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(final_proteins)

print("Saved positive.tsv:", len(final_proteins), "proteins")

for kingdom in ["Metazoa", "Fungi", "Plants", "Other"]:
    count = sum(p["kingdom"] == kingdom for p in final_proteins)
    print(kingdom + ":", count)