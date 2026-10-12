import pandas as pd
from Bio import SeqIO

def extract_cleavage_windows():
    train_df = pd.read_csv("data/splits/train.tsv", sep="\t")
    pos_train_df = train_df[train_df["label"] == 1]

    # Load sequences from FASTA
    seq_dict = {}
    for rec in SeqIO.parse("data/splits/train_pos.fasta", "fasta"):
        # Handle various accession formats in headers
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id.split()[0]
        seq_dict[acc] = str(rec.seq)

    windows = []
    for _, row in pos_train_df.iterrows():
        acc = str(row["accession"]).strip()
        cs = int(row["cleavage_site"])  # 1-based cleavage site position

        if acc in seq_dict:
            seq = seq_dict[acc]
            # Extract window [-13, +2] around cleavage site (0-based slicing)
            start_idx = cs - 13
            end_idx = cs + 2

            if start_idx >= 0 and end_idx <= len(seq):
                window = seq[start_idx:end_idx].upper()
                # Ensure exactly 15 characters and valid amino acid characters
                if len(window) == 15 and window.isalpha():
                    windows.append(window)

    output_path = "figures/cleavage_windows_13_2.txt"
    with open(output_path, "w") as f:
        for w in windows:
            f.write(w + "\n")

    print(f"Extracted {len(windows)} valid 15-aa windows to {output_path}")

if __name__ == "__main__":
    extract_cleavage_windows()