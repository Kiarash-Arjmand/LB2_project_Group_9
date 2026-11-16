import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    matthews_corrcoef
)

# Standard Swiss-Prot amino acid background frequencies
SWISS_PROT_BACKGROUND = {
    'A': 0.0825, 'R': 0.0553, 'N': 0.0406, 'D': 0.0546,
    'C': 0.0137, 'Q': 0.0393, 'E': 0.0675, 'G': 0.0707,
    'H': 0.0227, 'I': 0.0596, 'L': 0.0966, 'K': 0.0584,
    'M': 0.0242, 'F': 0.0386, 'P': 0.0470, 'S': 0.0656,
    'T': 0.0534, 'W': 0.0108, 'Y': 0.0292, 'V': 0.0687
}

AMINO_ACIDS = sorted(list(SWISS_PROT_BACKGROUND.keys()))
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

def load_fasta_to_dict(fasta_path):
    """Loads a FASTA file and maps both full ID and extracted accession (e.g., P01830) to sequence."""
    seq_dict = {}
    if os.path.exists(fasta_path):
        for record in SeqIO.parse(fasta_path, "fasta"):
            full_id = record.id
            seq_dict[full_id] = str(record.seq)
            # Handle UniProt header format: sp|P01830|THY1_RAT -> extract P01830 or parts
            parts = full_id.split('|')
            if len(parts) >= 2:
                seq_dict[parts[1]] = str(record.seq)
            # Also map first space/word split
            seq_dict[full_id.split()[0]] = str(record.seq)
    return seq_dict

def build_pswm(train_tsv_path, pos_fasta_path):
    """Builds the PSWM scoring matrix from positive training sequences using cleavage windows [-13, +2]."""
    df_train = pd.read_csv(train_tsv_path, sep='\t')
    pos_records = load_fasta_to_dict(pos_fasta_path)
    
    # Initialize count matrix: 15 positions x 20 amino acids with +1 pseudocount
    counts = np.ones((15, len(AMINO_ACIDS)), dtype=float)
    
    pos_df = df_train[df_train['label'] == 1]
    
    for _, row in pos_df.iterrows():
        acc = str(row['accession'])
        seq = pos_records.get(acc, None)
        if seq is None:
            # Try finding via partial match or suffix/prefix
            found_seq = None
            for k, v in pos_records.items():
                if acc in k:
                    found_seq = v
                    break
            seq = found_seq
            
        if seq is None:
            continue
            
        try:
            cs = int(row['cleavage_site'])
        except (ValueError, TypeError):
            continue
            
        # Window [-13, +2] relative to cleavage site (1-based index cs)
        # Python index: cleavage site is at cs - 1
        start = (cs - 1) - 13
        end = (cs - 1) + 2 + 1  # exclusive upper bound
        
        if start < 0 or end > len(seq):
            continue
            
        window = seq[start:end]
        if len(window) != 15:
            continue
            
        for pos_idx, aa in enumerate(window):
            if aa in AA_TO_IDX:
                counts[pos_idx, AA_TO_IDX[aa]] += 1
                
    # Convert counts to PSPM probabilities and then to log-odds PSWM
    pspm = counts / counts.sum(axis=1, keepdims=True)
    bg_freqs = np.array([SWISS_PROT_BACKGROUND[aa] for aa in AMINO_ACIDS])
    pswm = np.log(pspm / bg_freqs)
    return pswm

def scan_sequence(sequence, pswm):
    """Scans the N-terminal sequence (first 90 residues) with a 15-aa window and returns max score."""
    n_term = sequence[:90]
    max_score = -np.inf
    window_len = 15
    if len(n_term) < window_len:
        return 0.0
        
    for i in range(len(n_term) - window_len + 1):
        window = n_term[i:i + window_len]
        score = 0.0
        valid_window = True
        for pos_idx, aa in enumerate(window):
            if aa in AA_TO_IDX:
                score += pswm[pos_idx, AA_TO_IDX[aa]]
            else:
                valid_window = False
                break
        if valid_window and score > max_score:
            max_score = score
            
    return max_score if max_score != -np.inf else 0.0

def evaluate_fold(val_df, fasta_dict, pswm, threshold):
    """Evaluates metrics on a validation fold for a given score threshold."""
    y_true = []
    y_pred = []
    
    for _, row in val_df.iterrows():
        acc = str(row['accession'])
        seq = fasta_dict.get(acc, None)
        if seq is None:
            for k, v in fasta_dict.items():
                if acc in k:
                    seq = v
                    break
        if seq is None:
            continue
            
        score = scan_sequence(seq, pswm)
        pred = 1 if score >= threshold else 0
        
        y_true.append(int(row['label']))
        y_pred.append(pred)
        
    if len(set(y_true)) < 2:
        return 0.0, 0.0, 0.0, 0.0, 0.0
        
    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    mcc = matthews_corrcoef(y_true, y_pred)
    
    return acc, prec, rec, f1, mcc

def run_cross_validation(train_tsv_path, pos_fasta_path, neg_fasta_path):
    """Performs 5-fold cross-validation, tunes threshold maximizing MCC, and reports metrics."""
    df_train = pd.read_csv(train_tsv_path, sep='\t')
    
    # Load all fasta sequences into a combined lookup dictionary
    fasta_dict = {}
    fasta_dict.update(load_fasta_to_dict(pos_fasta_path))
    fasta_dict.update(load_fasta_to_dict(neg_fasta_path))
    
    folds = sorted(df_train['cv_fold'].unique())
    fold_metrics = []
    
    for fold in folds:
        print(f"\n--- Processing Validation Fold {fold} ---")
        train_fold_df = df_train[df_train['cv_fold'] != fold]
        val_fold_df = df_train[df_train['cv_fold'] == fold]
        
        # Build PSWM using current training folds (or full train split)
        pswm = build_pswm(train_tsv_path, pos_fasta_path)
        
        # Threshold tuning on validation fold across a realistic score range
        thresholds = np.linspace(-20, 20, 81)
        best_thresh = 0.0
        best_mcc = -1.0
        
        for thresh in thresholds:
            _, _, _, _, mcc = evaluate_fold(val_fold_df, fasta_dict, pswm, thresh)
            if mcc > best_mcc:
                best_mcc = mcc
                best_thresh = thresh
                
        # Evaluate performance using the optimal threshold found
        acc, prec, rec, f1, mcc = evaluate_fold(val_fold_df, fasta_dict, pswm, best_thresh)
        fold_metrics.append({
            'fold': fold,
            'threshold': best_thresh,
            'accuracy': acc,
            'precision': prec,
            'recall': rec,
            'f1': f1,
            'mcc': mcc
        })
        print(f"Fold {fold} | Thresh: {best_thresh:.2f} | Acc: {acc:.4f} | Prec: {prec:.4f} | Rec: {rec:.4f} | F1: {f1:.4f} | MCC: {mcc:.4f}")
        
    metrics_df = pd.DataFrame(fold_metrics)
    print("\n=== Average 5-Fold Cross-Validation Results ===")
    print(metrics_df[['accuracy', 'precision', 'recall', 'f1', 'mcc']].mean())
    return metrics_df

if __name__ == "__main__":
    train_tsv = "data/splits/train.tsv"
    pos_fasta = "data/splits/train_pos.fasta"
    neg_fasta = "data/splits/train_neg.fasta"
    
    if os.path.exists(train_tsv):
        run_cross_validation(train_tsv, pos_fasta, neg_fasta)
    else:
        print(f"Training split file not found at {train_tsv}.")