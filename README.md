# Laboratory of Bioinformatics II — Group 9

We compare a von Heijne position-specific weight matrix (PSWM) with a support vector machine (SVM) for detecting secretory signal peptides in reviewed eukaryotic UniProtKB proteins.

| Stage | Positive | Negative | Total |
|---|---:|---:|---:|
| Downloaded candidates | 2,972 | 20,974 | 23,946 |
| After annotation selection and preparation | 2,956 | 20,974 | 23,930 |
| MMseqs2 representatives | 1,104 | 9,082 | 10,186 |
| Training (five stratified folds) | 883 | 7,265 | 8,148 |
| Independent benchmark | 221 | 1,817 | 2,038 |

## Reproduce the analysis

Run the commands from the repository root. The data, cluster outputs and splits are checked in. Python 3.12 was used for this run.

```bash
python -m pip install -r requirements.txt
python -m src.visualization
python -m src.models
python -m src.svm
python -m src.benchmark
```

The SVM run includes nested cross-validation, so it takes longer than the other steps. It uses four workers by default; set the `LB2_JOBS` environment variable to reduce that number if needed. The exploratory `src/features.py` module is separate; `src/svm.py` contains the features used for the reported results.

To regenerate the two committed cleavage logos, install [Ghostscript](https://www.ghostscript.com/) and make sure it is on your PATH. Then run `python -m pip install weblogo==3.7.12 setuptools==80.9.0` followed by `python -m src.visualization --logos`.

### Dataset preparation

- Raw inputs are in `data/raw/`. Run `python filter_positives.py` and `python filter_negatives.py` to rebuild the processed TSVs and FASTAs from those inputs. The negatives were selected with an experimental compartment restriction in the UniProt search, so processing the export retains the same 20,974 entries. The query, release and checksums are in `data/raw/negative_collection.json`. `collect_negative_candidates.py` queries the live database again and can return different data in a later release; the committed inputs are used here.
- Clustered separately by class with MMseqs2 using `--min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`. The cluster memberships and representatives are in `data/processed/`. Run `python -m src.preprocessing` to regenerate the 80/20 split and five training folds (seed 42) from those representatives.
- `src.visualization` generates the length distributions, kingdom and species plots, and a comparison of amino acid composition between the annotated signal peptides and [Swiss-Prot release 2026_04, section 6.1](https://web.expasy.org/docs/relnotes/relstat.html). The species plots show the 12 most frequent species plus one pooled category. The cleavage logo inputs preserve 13 residues before the cut and two after it.

### Evaluation

- **von Heijne:** Count residues at the 15 cleavage-window positions with a pseudocount of 1, convert the matrix to log-odds against Swiss-Prot, then take the highest valid window score within the first 90 residues of each sequence. In each of five runs, three folds fit the matrix, the next fold selects an MCC threshold, and a separate fold measures performance. For the final run, choose one MCC cutoff from training-only out-of-fold scores and refit the matrix on all training positives before scoring the benchmark.
- 
- **SVM:** Use N-terminal amino acid composition, hydropathy and charge. Nested cross-validation tunes feature count and RBF parameters by MCC; the feature ranking and scaler are fitted inside the training portion of each fold. Tune on all training folds, refit on the full training set and then score the benchmark.
- The two CV protocols have different training sizes and tuning schemes, so the **shared independent benchmark** is the direct comparison. Its metrics, confusion counts, early transmembrane helix error rates and ROC/precision–recall curves are in `results/`.
- **SVM detector:** The SVM predicts whether a eukaryotic protein carries a secretory signal peptide (1 = SP, 0 = no SP) from its N-terminal sequence only. The cleavage site is never used as a feature: it is used only to compute the window length n from the training positives.
**Features**
N-terminal composition: the fraction of each of the 20 amino acids in the first n residues, where n is the median signal-peptide length of the current training set (22 on the full training set; recomputed in every fold).
Hydropathy: Kyte–Doolittle profile (window 5) over the first 40 residues: maximum, mean and position of the peak.
Charge: net charge (K, R positive; D, E negative) of residues 1–10.
**Model and protocol**
Features → Random-Forest importance ranking (500 trees, fixed seed) → top-k features → StandardScaler → RBF SVC(class_weight="balanced").
Nested cross-validation on the five predefined training folds. The outer loop evaluates (each fold is the test set once); the inner 4-fold loop chooses k, C and γ by mean MCC. The test fold is never used for a choice.
n, the RF ranking, the scaler and the SVM are fitted on training data only, separately in every fold.
Final model: 5-fold CV over the grid on the whole training set, refit on the whole training set, then scored once on the independent benchmark.
Results are reported as mean ± standard error over the five outer folds (SE = σ/√5).
**Results**
Nested cross-validation (mean ± SE over 5 outer folds): see results/ for the per-fold tables and chosen parameters. Benchmark metrics are in the table under Independent benchmark results.
**Error analysis**
False positives: an N-terminal transmembrane helix is the main confounder. Its hydrophobic stretch resembles the h-region of a signal peptide, so the false-positive rate is about ten times higher for proteins with a TM helix starting within the first 90 residues (10.3%) than for those without (1.2%).
False negatives: compared with true positives, the missed signal peptides tend to have less alanine and leucine in the N-terminus. With only a few dozen cases this is a tendency, not a result.
Taxonomy: the benchmark contains few positives outside Metazoa (roughly 16 plant, 17 fungal and 5 "other"), so per-kingdom error rates are very uncertain.
**Limitations**
The benchmark has 221 positives, so its metrics carry a wide uncertainty (a few hundredths of MCC). Differences smaller than that should not be over-interpreted.
After clustering, Fungi are about 28% of the negatives but only about 9% of the positives (Metazoa: 53% vs 79%). A composition-based model may partly learn taxonomic signal; the per-kingdom analysis is there to check this.
Some negatives are multi-localised proteins that also carry a secretory annotation but no signal peptide.
A few sequences contain non-standard residues (U, B, X); they are ignored when computing composition.
The Random-Forest selection kept all features in every fold, so in practice it did not reduce the feature set.

## Independent benchmark results

| Method | MCC | Precision | Recall | F1 | Accuracy | ROC AUC |
|---|---:|---:|---:|---:|---:|---:|
| von Heijne | 0.657 | 0.708 | 0.679 | 0.693 | 0.935 | 0.952 |
| SVM | 0.853 | 0.840 | 0.900 | 0.869 | 0.971 | 0.981 |

The von Heijne model called 28 of 185 negatives with a transmembrane helix starting within the first 90 residues positive (15.1%), compared with 34 of 1,632 negatives without one (2.1%). The SVM rates were 19 of 185 (10.3%) and 19 of 1,632 (1.2%), respectively. Both methods were less specific on the early transmembrane subset; the SVM made fewer false positives and fewer false negatives overall on this benchmark.

See [benchmark metrics and confusion counts](results/benchmark_metrics.tsv), [subgroup errors](results/benchmark_tm_false_positives.tsv) and [ROC and precision–recall curves](results/benchmark_curves.png).
