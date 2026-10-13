# Laboratory of Bioinformatics II — Group 9

We are comparing two methods for predicting secretory signal peptides in eukaryotic proteins: the **von Heijne method (1986)** and a **Support Vector Machine (SVM)** classifier.

---

## Current Progress

Data collection, annotation filtering, sequence redundancy reduction (using MMseqs2), train/benchmark dataset splits, and statistical dataset analyses are **100% complete**.

1. **Part A (Data Collection)**: Filtered reviewed, non-fragment eukaryotic proteins with protein-level evidence from UniProtKB.
2. **Part B (Redundancy Reduction & Splits)**:
   - Clustered positive and negative FASTA sequences using **MMseqs2** (`--cluster-mode 1`, 30% sequence identity, 40% coverage).
   - Reduced candidate counts to **1,104 non-redundant positives** and **9,082 non-redundant negatives**.
   - Generated an 80% / 20% stratified train/benchmark split with **5-fold cross-validation fold assignments** (`cv_fold` 0–4).
3. **Part C (Data Analysis)**: Generated sequence length distributions, kingdom counts, and extracted 15-aa cleavage site windows `[-13, +2]` for WebLogo sequence logo generation (`figures/cleavage_site_logo.png`).

---

## Project File Structure

Paths are relative to the repository root.

| Path | Description |
| :--- | :--- |
| `filter_positives.py` | Positive candidate annotation filtering script. |
| `collect_negative_candidates.py` | Script to retrieve course-specific negative collection. |
| `filter_negatives.py` | Negative candidate annotation filtering script. |
| `src/preprocessing.py` | Parses MMseqs2 outputs, builds non-redundant TSVs, and generates 80/20 splits with 5-fold CV assignments. |
| `src/visualization.py` | Generates summary figures and extracts 15-aa cleavage windows `[-13, +2]` for WebLogo. |
| `src/features.py` | *(Pending)* Feature extraction routines for SVM (AAC, physicochemical features, sliding windows). |
| `src/models.py` | *(Pending)* Implementation of von Heijne PSWM scoring and SVM model wrappers. |
| `data/raw/` | Original candidate annotations and FASTA downloads from UniProtKB. |
| `data/processed/` | Filtered entries, annotation logs, and MMseqs2 clustering outputs (`pos_clustered_*`, `neg_clustered_*`). |
| `data/splits/train.tsv` | Non-redundant training metadata (80%, 8,149 proteins) with `cv_fold` column (values 0–4). |
| `data/splits/benchmark.tsv` | Non-redundant independent benchmarking metadata (20%, 2,037 proteins). |
| `data/splits/*.fasta` | Corresponding FASTA files for positive/negative training and benchmarking splits. |
| `figures/` | Saved statistical plots (`*_sp_length_distribution.png`, `*_protein_length_comparison.png`, `*_kingdom_distribution.png`) and `cleavage_site_logo.png`. |

---

## Dataset Statistics

### Filtering & Redundancy Reduction Counts

| Dataset | UniProt Candidates | Retained (Part A) | Non-Redundant Reps (MMseqs2) |
| :--- | :---: | :---: | :---: |
| **Positive** | 2,972 | 2,956 | **1,104** |
| **Negative** | 20,974 | 20,974 | **9,082** |

### Final Dataset Splits (80% Train / 20% Benchmark)

| Split File | Positives | Negatives | Total Entries |
| :--- | :---: | :---: | :---: |
| `data/splits/train.tsv` | 883 | 7,266 | **8,149** |
| `data/splits/benchmark.tsv` | 221 | 1,816 | **2,037** |
| **Total Non-Redundant** | **1,104** | **9,082** | **10,186** |

---

## Next Steps

### Task 1: Implement the von Heijne PSWM Algorithm (`src/models.py`)

1. **Build the PSWM Scoring Matrix**:
   - Load positive sequences from `data/splits/train_pos.fasta` and extract the 15-aa cleavage window `[-13, +2]` using `cleavage_site` coordinates in `data/splits/train.tsv`.
   - Calculate the Position-Specific Probability Matrix (PSPM) across the 15 positions ($15 \times 20$ amino acids). Apply a **pseudocount of +1** to every cell.
   - Convert to a Position-Specific Weight Matrix (PSWM) using log-odds against Swiss-Prot background frequencies ($b_k$):
     $$W_{k, j} = \log \left( \frac{M_{k, j}}{b_k} \right)$$

2. **Sequence Scanning Function**:
   - Write a function `scan_sequence(sequence, pswm)` that takes an N-terminal sequence (first 90 residues), slides the 15-aa window across positions 1 to 76, calculates the score for each window, and returns the **maximum score** found.

3. **5-Fold Cross-Validation & Threshold Tuning**:
   - Run 5-fold cross-validation using the `cv_fold` column (0–4) in `data/splits/train.tsv`.
   - On the training folds, construct the PSWM and evaluate on the validation fold across candidate decision thresholds ($\theta$).
   - Determine the optimal threshold $\theta^*$ that maximizes the **Matthews Correlation Coefficient (MCC)**.
   - Report average **Accuracy, Precision, Recall, F1-Score, and MCC** across the 5 validation folds.

---

### Task 2: Feature Extraction & Encoding (`src/features.py`)

1. **Extract Global Sequence Features**:
   - **Amino Acid Composition (AAC)**: Compute a 20-dimensional frequency vector across the N-terminal region (first 30–50 residues) for each protein in `train.tsv` and `benchmark.tsv`.
   - **Physicochemical Properties**: Extract net charge, Kyte-Doolittle hydropathy index, and polar/apolar ratios.

2. **Extract Local / Window Features**:
   - Generate one-hot encoded matrices or sliding-window features for the N-terminal region.

3. **Feature Standardization**:
   - Fit `sklearn.preprocessing.StandardScaler` **only on the training split** (`data/splits/train.tsv`) and transform both `train` and `benchmark` feature sets to prevent data leakage.

---

### Task 3: SVM Classifier Training & Benchmarking Evaluation

1. **Model Tuning**:
   - Train an `sklearn.svm.SVC` (testing `kernel='rbf'` and `kernel='linear'`, with `class_weight='balanced'`) on the training folds (`data/splits/train.tsv`).
   - Tune hyperparameters ($C, \gamma$) via cross-validation.

2. **Holdout Evaluation**:
   - Evaluate both the optimized **von Heijne model** and the **SVM classifier** on the unseen 20% holdout set (`data/splits/benchmark.tsv`).
   - Specifically measure and compare the **false-positive rate** on negative proteins containing N-terminal transmembrane helices (`tm_helix_first_90 == True`).

---

## 🛠️ How to Run Existing Code

```bash
# 1. Regenerate non-redundant splits and CV assignments
python src/preprocessing.py

# 2. Regenerate data analysis figures and cleavage site windows
python src/visualization.py