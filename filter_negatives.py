import csv
import re
from pathlib import Path

project_folder = Path(__file__).resolve().parent
input_file = project_folder / "data" / "raw" / "negative_candidates.tsv"
output_folder = project_folder / "data" / "processed"
output_folder.mkdir(parents=True, exist_ok=True)

with open(input_file, encoding="utf-8", newline="") as file:
    reader = csv.DictReader(file, delimiter="\t")
    columns = list(reader.fieldnames)
    proteins = list(reader)

# Allowed locations, including their subcompartments.
location_pattern = (
    r"^(?:\[[^\]]+\]:\s*)?"
    r"(?:Cytoplasm,\s*cytosol|Nucleus|Mitochondrion|"
    r"Plastid|Peroxisome|Cell membrane)(?=\s|,|$)"
)

kept_proteins = []

for protein in proteins:
    if protein["Signal peptide"].strip():
        continue

    annotation = protein["Subcellular location [CC]"]
    valid_location = False

    # Handle separate location annotations, including isoform annotations.
    for block in annotation.split("SUBCELLULAR LOCATION:")[1:]:
        # Free-text notes are not location evidence.
        block = block.split("Note=", 1)[0]

        # Pair each location statement with its own evidence.
        statements = re.findall(r"([^{}.;]+)\{([^{}]+)\}", block)

        for location, evidence in statements:
            allowed = re.match(location_pattern, location.strip())
            experimental = "ECO:0000269" in evidence

            if allowed and experimental:
                valid_location = True

    if valid_location:
        kept_proteins.append(protein)

# Save the selected annotations.
with open(output_folder / "negative_filtered_annotations.tsv", "w",
          encoding="utf-8", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(kept_proteins)

# Save accessions for downloading only the selected sequences.
with open(output_folder / "negative_accessions.txt", "w",
          encoding="utf-8") as file:
    for protein in kept_proteins:
        file.write(protein["Entry"] + "\n")

print("Candidates:", len(proteins))
print("Proteins kept:", len(kept_proteins))
print("Proteins excluded:", len(proteins) - len(kept_proteins))
print("Saved selected annotations and accession list.")

def get_kingdom(lineage):
    if "Metazoa" in lineage:
        return "Metazoa"
    if "Fungi" in lineage:
        return "Fungi"
    if "Viridiplantae" in lineage:
        return "Plants"
    return "Other"


def has_early_tm_helix(annotation):
    features = re.split(r"(?=TRANSMEM )", annotation)

    for feature in features:
        start = re.match(r"TRANSMEM (\d+)\.\.", feature)
        note = re.search(r'/note="([^"]+)"', feature)

        if start is None or note is None:
            continue

        position = int(start.group(1))
        description = note.group(1).lower()

        is_helix = (
            "helical" in description
            or "signal-anchor" in description
        )

        if 1 <= position <= 90 and is_helix:
            return True

    return False


final_columns = [
    "accession",
    "organism",
    "kingdom",
    "protein_length",
    "tm_helix_first_90",
]

helix_count = 0

with open(
    output_folder / "negative.tsv",
    "w",
    encoding="utf-8",
    newline=""
) as file:
    writer = csv.writer(file, delimiter="\t")
    writer.writerow(final_columns)

    for protein in kept_proteins:
        has_helix = has_early_tm_helix(protein["Transmembrane"])
        helix_count += int(has_helix)

        writer.writerow([
            protein["Entry"],
            protein["Organism"],
            get_kingdom(protein["Taxonomic lineage"]),
            protein["Length"],
            str(has_helix).lower(),
        ])

print("Saved negative.tsv:", len(kept_proteins), "proteins")
print("TM helix starting within residues 1–90:", helix_count)
print("Without such an annotated helix:", len(kept_proteins) - helix_count)