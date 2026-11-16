import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.preprocessing import StandardScaler

# Standard 20 amino acids
AMINO_ACIDS = sorted(list("ACDEFGHIKLMNPQRSTVWY"))
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# Kyte-Doolittle hydropathy index values
KYTE_DOOLITTLE = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Polar and Apolar amino acid groups
POLAR_AAS = set("STNQCHY")
APOLAR_AAS = set("AVLIMFWP")

def load_fasta_to_dict(fasta_paths):
    """Loads multiple FASTA files and creates a lookup dictionary for sequences."""
    seq_dict = {}
    for path in fasta_paths:
        if os.path.exists(path):
            for record in SeqIO.parse(path, "fasta"):
                full_id = record.id
                seq_dict[full_id] = str(record.seq)
                parts = full_id.split('|')
                if len(parts) >= 2:
                    seq_dict[parts[1]] = str(record.seq)
                seq_dict[full_id.split()[0]] = str(record.seq)
    return seq_dict

def extract_aac(sequence, n_term_len=50):
    """
    1. Global Sequence Feature: Amino Acid Composition (AAC)
    Computes a 20-dimensional frequency vector across the N-terminal region.
    """
    n_term = sequence[:n_term_len]
    counts = np.zeros(20, dtype=float)
    total = len(n_term)
    if total == 0:
        return counts
    for aa in n_term:
        if aa in AA_TO_IDX:
            counts[AA_TO_IDX[aa]] += 1.0
    return counts / total

def extract_physicochemical(sequence, n_term_len=50):
    """
    1. Global Sequence Feature: Physicochemical Properties
    Extracts net charge, mean Kyte-Doolittle hydropathy index, and polar/apolar ratios.
    """
    n_term = sequence[:n_term_len]
    total = len(n_term)
    if total == 0:
        return [0.0, 0.0, 0.0, 0.0]
        
    net_charge = 0.0
    hydropathy_sum = 0.0
    polar_count = 0
    apolar_count = 0
    
    for aa in n_term:
        # Net charge at physiological pH: K, R (+1), D, E (-1)
        if aa in ('K', 'R'):
            net_charge += 1.0
        elif aa in ('D', 'E'):
            net_charge -= 1.0
            
        # Hydropathy sum
        hydropathy_sum += KYTE_DOOLITTLE.get(aa, 0.0)
        
        # Polar and Apolar counts
        if aa in POLAR_AAS:
            polar_count += 1
        elif aa in APOLAR_AAS:
            apolar_count += 1
            
    mean_hydropathy = hydropathy_sum / total
    polar_ratio = polar_count / total
    apolar_ratio = apolar_count / total
    
    return [net_charge, mean_hydropathy, polar_ratio, apolar_ratio]

def extract_local_window_features(sequence, max_len=50):
    """
    2. Local / Window Feature: One-hot encoded matrix for the N-terminal region.
    """
    n_term = sequence[:max_len]
    matrix = np.zeros((max_len, len(AMINO_ACIDS)), dtype=float)
    for i, aa in enumerate(n_term):
        if aa in AA_TO_IDX:
            matrix[i, AA_TO_IDX[aa]] = 1.0
    return matrix.flatten()

def extract_features_for_dataset(tsv_path, pos_fasta_path, neg_fasta_path, n_term_len=50):
    """Extracts all feature components for all records in a given dataset split."""
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"Dataset split not found at {tsv_path}")
        
    df = pd.read_csv(tsv_path, sep='\t')
    fasta_dict = load_fasta_to_dict([pos_fasta_path, neg_fasta_path])
    
    features_list = []
    labels = []
    
    for _, row in df.iterrows():
        acc = str(row['accession'])
        seq = fasta_dict.get(acc, None)
        if seq is None:
            for k, v in fasta_dict.items():
                if acc in k:
                    seq = v
                    break
        if seq is None:
            continue
            
        aac = extract_aac(seq, n_term_len=n_term_len)
        phys_chem = extract_physicochemical(seq, n_term_len=n_term_len)
        local_feat = extract_local_window_features(seq, max_len=n_term_len)
        
        # Concatenate global and local features into a single feature vector
        combined = np.concatenate([aac, phys_chem, local_feat])
        features_list.append(combined)
        labels.append(int(row['label']))
        
    return np.array(features_list), np.array(labels), df

def prepare_and_scale_features():
    """
    3. Feature Standardization:
    Fits StandardScaler ONLY on training data and transforms both train and benchmark sets.
    """
    train_tsv = "data/splits/train.tsv"
    train_pos = "data/splits/train_pos.fasta"
    train_neg = "data/splits/train_neg.fasta"
    
    bench_tsv = "data/splits/benchmark.tsv"
    bench_pos = "data/splits/benchmark_pos.fasta"
    bench_neg = "data/splits/benchmark_neg.fasta"
    
    print("Extracting features for training set...")
    X_train, y_train, df_train = extract_features_for_dataset(train_tsv, train_pos, train_neg, n_term_len=50)
    
    print("Extracting features for benchmark set...")
    X_bench, y_bench, df_bench = extract_features_for_dataset(bench_tsv, bench_pos, bench_neg, n_term_len=50)
    
    # Feature Standardization (preventing data leakage by fitting ONLY on train)
    print("Fitting StandardScaler on training features and transforming benchmark...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_bench_scaled = scaler.transform(X_bench)
    
    print(f"Training feature matrix shape: {X_train_scaled.shape}")
    print(f"Benchmark feature matrix shape: {X_bench_scaled.shape}")
    
    return X_train_scaled, y_train, X_bench_scaled, y_bench, df_train, df_bench

if __name__ == "__main__":
    prepare_and_scale_features()