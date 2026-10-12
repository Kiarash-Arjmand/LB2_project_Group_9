import os
import re
import pandas as pd
from Bio import SeqIO
from sklearn.model_selection import StratifiedKFold, train_test_split

def extract_uniprot_accession(header_str):
    """
    Extracts a standard UniProt accession (e.g., P12345 or Q9Y6X9) 
    from various FASTA header formats.
    """
    # Try matching pipe-separated format: >sp|P12345|... or >tr|P12345|...
    pipe_match = re.search(r'(?:sp|tr)\|([A-Z0-9]{6,10})\|', header_str)
    if pipe_match:
        return pipe_match.group(1)
    
    # Try matching standard standalone accession at the start of string
    acc_match = re.search(r'([A-Z0-9]{6,10})', header_str)
    if acc_match:
        return acc_match.group(1)
        
    return header_str.split()[0]

def get_fasta_ids(fasta_path):
    """Extract clean Accession IDs from a FASTA file using BioPython."""
    accessions = set()
    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = extract_uniprot_accession(record.description)
        accessions.add(acc)
    return accessions

# 1. Load representative IDs from MMseqs2 output files
pos_rep_ids = get_fasta_ids("data/processed/pos_clustered_rep_seq.fasta")
neg_rep_ids = get_fasta_ids("data/processed/neg_clustered_rep_seq.fasta")

print(f"Extracted {len(pos_rep_ids)} positive representative IDs from FASTA.")
print(f"Extracted {len(neg_rep_ids)} negative representative IDs from FASTA.")

# 2. Load full metadata TSVs
pos_df = pd.read_csv("data/processed/positive.tsv", sep="\t")
neg_df = pd.read_csv("data/processed/negative.tsv", sep="\t")

# Clean accession column in metadata TSVs (just in case of whitespace/formatting differences)
pos_df["accession"] = pos_df["accession"].astype(str).str.strip()
neg_df["accession"] = neg_df["accession"].astype(str).str.strip()

# 3. Filter metadata TSVs to retain ONLY representative entries
pos_nr = pos_df[pos_df["accession"].isin(pos_rep_ids)].copy()
neg_nr = neg_df[neg_df["accession"].isin(neg_rep_ids)].copy()

print(f"Non-redundant Positives matched: {len(pos_nr)} (from {len(pos_df)})")
print(f"Non-redundant Negatives matched: {len(neg_nr)} (from {len(neg_df)})")

# Sanity check: Ensure non-zero rows
if len(pos_nr) == 0 or len(neg_nr) == 0:
    raise ValueError(
        "Filtering resulted in 0 entries. Check sample accession IDs between TSV and FASTA:\n"
        f"Sample TSV pos accession: {pos_df['accession'].iloc[0] if len(pos_df)>0 else 'None'}\n"
        f"Sample FASTA pos accession: {list(pos_rep_ids)[0] if pos_rep_ids else 'None'}"
    )

# 4. Add binary class targets (1 = Positive, 0 = Negative)
pos_nr["label"] = 1
neg_nr["label"] = 0

# 5. Combine Positives and Negatives
combined_df = pd.concat([pos_nr, neg_nr], ignore_index=True)

# 6. Perform 80% Train / 20% Benchmarking Split (Stratified)
train_df, benchmark_df = train_test_split(
    combined_df, test_size=0.20, random_state=42, stratify=combined_df["label"]
)

# 7. Set up 5-Fold Stratified Cross-Validation on the Training Set
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
train_df["cv_fold"] = -1

for fold, (_, val_idx) in enumerate(skf.split(train_df, train_df["label"])):
    train_df.iloc[val_idx, train_df.columns.get_loc("cv_fold")] = fold

# 8. Save output metadata files to data/splits/
os.makedirs("data/splits", exist_ok=True)
train_df.to_csv("data/splits/train.tsv", sep="\t", index=False)
benchmark_df.to_csv("data/splits/benchmark.tsv", sep="\t", index=False)

# 9. Save matching representative FASTA files to data/splits/
def save_filtered_fasta(input_fasta, valid_ids, output_fasta):
    records = []
    for rec in SeqIO.parse(input_fasta, "fasta"):
        acc = extract_uniprot_accession(rec.description)
        if acc in valid_ids:
            records.append(rec)
    SeqIO.write(records, output_fasta, "fasta")

save_filtered_fasta("data/processed/pos_clustered_rep_seq.fasta", set(train_df[train_df["label"]==1]["accession"]), "data/splits/train_pos.fasta")
save_filtered_fasta("data/processed/neg_clustered_rep_seq.fasta", set(train_df[train_df["label"]==0]["accession"]), "data/splits/train_neg.fasta")
save_filtered_fasta("data/processed/pos_clustered_rep_seq.fasta", set(benchmark_df[benchmark_df["label"]==1]["accession"]), "data/splits/benchmark_pos.fasta")
save_filtered_fasta("data/processed/neg_clustered_rep_seq.fasta", set(benchmark_df[benchmark_df["label"]==0]["accession"]), "data/splits/benchmark_neg.fasta")

print("Part B complete! All train/benchmark split files created in data/splits/")