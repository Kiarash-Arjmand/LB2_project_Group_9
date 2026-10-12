# Laboratory of Bioinformatics II — Group 9

This project compares two methods for predicting secretory signal peptides in eukaryotic proteins: the von Heijne method and a support vector machine (SVM).

## Progress

Data collection, annotation filtering, sequence redundancy reduction (using MMseqs2), train/benchmark dataset splits, and statistical dataset analyses are complete.

The preliminary dataset contained 2,956 positive proteins and 16,852 negative proteins. Following sequence clustering with MMseqs2 (using 30% sequence identity and 40% alignment coverage), the non-redundant dataset comprises **1,104 positive proteins** and **8,113 negative proteins**.

The non-redundant data have been randomly split into an 80% training set (with 5-fold cross-validation fold assignments) and a 20% holdout benchmarking set. Descriptive visualizations (length distributions, taxonomy, and cleavage site logos) have been generated and saved to `figures/`. Predictor development (von Heijne PSWM and SVM feature extraction) is pending.

## Files

All paths are relative to the repository root.

| File or folder | Contents |
| :--- | :--- |
| `filter_positives.py` | Positive annotation filtering and output preparation. |
| `filter_negatives.py` | Negative annotation filtering and output preparation. |
| `src/preprocessing.py` | MMseqs2 parsing, non-redundant dataset construction, 80/20 train/benchmark splitting, and 5-fold cross-validation fold assignment. |
| `src/visualization.py` | Data analysis plotting script generating length distributions, taxonomic counts, and cleavage site window extractions. |
| `data/raw/` | Original candidate annotations and fasta downloads from UniProtKB. |
| `data/processed/` | Retained positive/negative proteins, annotation logs, and MMseqs2 clustering output files (`pos_clustered_*`, `neg_clustered_*`). |
| `data/splits/train.tsv` | Non-redundant training metadata (80%) including `cv_fold` (0–4) column assignments. |
| `data/splits/benchmark.tsv` | Non-redundant independent benchmarking metadata (20%). |
| `data/splits/*.fasta` | Corresponding positive and negative sequence fasta files for training and benchmark splits. |
| `figures/` | Saved statistical plots (`*_sp_length_distribution.png`, `*_protein_length_comparison.png`, `*_kingdom_distribution.png`) and the cleavage site sequence logo (`cleavage_site_logo.png`). |
| `notebooks/` | Reserved for subsequent exploratory code and model evaluation notebooks. |

## Data Collection

The candidate data were downloaded from UniProtKB/Swiss-Prot on 6 October 2026. Both searches select reviewed eukaryotic proteins (taxonomy ID 2759), exclude fragments, require a sequence length of at least 40 residues, and require evidence at protein level (existence:1).

### Filtering & Redundancy Reduction Results

| Dataset | Candidates | Retained (Part A) | Non-Redundant (MMseqs2 Reps) |
| :--- | :--- | :--- | :--- |
| **Positive** | 2,972 | 2,956 | **1,104** |
| **Negative** | 66,172 | 16,852 | **8,113** |

### Dataset Splits (80% Train / 20% Benchmark)

| Split | Positives | Negatives | Total |
| :--- | :--- | :--- | :--- |
| **Training Set (80%)** | 883 | 6,490 | 7,373 |
| **Benchmarking Set (20%)** | 221 | 1,623 | 1,844 |
| **Total Non-Redundant** | **1,104** | **8,113** | **9,217** |

## Next Steps

1. **Feature Extraction**: Extract amino acid composition, physicochemical properties, and position-specific sliding windows around the N-terminus.
2. **Model Implementation**:
   - Build von Heijne's Position-Specific Weight Matrix (PSWM) using cleavage site context $[-13, +2]$.
   - Train and tune a Support Vector Machine (SVM) classifier on the 5-fold cross-validation training set.
3. **Evaluation**: Compare cross-validation and benchmarking metrics (Accuracy, Precision, Recall, MCC, F1-score) between both predictors.