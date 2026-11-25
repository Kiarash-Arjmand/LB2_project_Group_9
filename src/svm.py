# Support Vector Machines with scikit-learn — signal-peptide detection
#   0. Imports
#   1. The dataset (train / benchmark, features)
#   2. A first SVM on the raw features
#   3. Standardizing the features
#   4. Looking at the trained model: support vectors, gamma, the decision function rebuilt by hand
#   5. Probabilities with predict_proba
#   6. What happens with a badly chosen gamma
#   7. Choosing the hyper-parameters: Random-Forest feature selection + nested cross-validation
#   8. Final model and benchmark (the benchmark set is used ONLY here, ONLY once)
#   9. Saving and reloading the model
#  10. Error analysis (false positives / false negatives)
#
# Run from the repository root:   python src/svm.py        (or open it cell by cell: "# %%" markers)
#
# Rule of the whole script: everything that is "learned" (n, the RF ranking, the scaler, the SVM)
# is fitted on TRAINING data only.

# %% 0. Imports
import glob
import gzip
import os
import pickle

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, f1_score, matthews_corrcoef,
                             precision_score, recall_score)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
# ------CONFIGURATION--------------------------------
DATA_DIR = "data/splits"
TRAIN_TSV = f"{DATA_DIR}/train.tsv"
BENCH_TSV = f"{DATA_DIR}/benchmark.tsv"
ID_COL, LABEL_COL, CLEAVAGE_COL, FOLD_COL = "accession", "label", "cleavage_site", "cv_fold"
TM_COL, KINGDOM_COL = "tm_helix_first_90", "kingdom"
FIG_DIR = "results/svm"
MODEL_FILE = f"{FIG_DIR}/svm_final_model.pkl.gz"
N_JOBS = int(os.environ.get("LB2_JOBS", "4"))   # parallel SVM fits
RUN_NESTED_CV = True      # section 7 takes a few minutes: set False to skip it
SHOW_PLOTS = False        # True: also open the figures on screen
os.makedirs(FIG_DIR, exist_ok=True)
if not SHOW_PLOTS:
    matplotlib.use("Agg")


def finish_plot(name):
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/{name}.png", dpi=150)
    if SHOW_PLOTS:
        plt.show()
    plt.close()

# %% 1. The dataset
# The TSV files hold the metadata, the sequences are in the FASTA files of data/splits/.
def read_fasta(paths):
    """FASTA -> {id: sequence}; handles UniProt headers like >sp|P12345|NAME_HUMAN ..."""
    seqs = {}
    for p in paths:
        name, buf = None, []
        for line in open(p):
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                if "|" in name:
                    name = name.split("|")[1]
                buf = []
            elif line:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf)
    return seqs


FASTA = read_fasta(glob.glob(f"{DATA_DIR}/*.fasta"))


def load_split(path):
    df = pd.read_csv(path, sep="\t")
    df["seq"] = df[ID_COL].map(FASTA).astype(str).str.upper()
    df["y"] = df[LABEL_COL].astype(int)                       # 1 = has a signal peptide, 0 = no SP
    df["sp_len"] = pd.to_numeric(df[CLEAVAGE_COL], errors="coerce")   # NaN for negatives
    if FOLD_COL in df.columns:
        df["fold"] = df[FOLD_COL].astype(int)
    assert df["seq"].ne("NAN").all(), "some sequences are missing from the FASTA files"
    return df


train = load_split(TRAIN_TSV)
bench = load_split(BENCH_TSV)
print("train :", train.shape, "| positives:", train.y.sum(), "| negatives:", (train.y == 0).sum())
print("bench :", bench.shape, "| positives:", bench.y.sum(), "| negatives:", (bench.y == 0).sum())
print("classes: 0 = no SP, 1 = SP  (always the positive class = 1 = signal peptide)")
print("positive fraction per CV fold:", train.groupby("fold").y.mean().round(3).to_dict())

# ---- Features -------------------------------------------------------------------------------
# Each protein -> a point in a D-dimensional space:
#   * amino-acid composition of the first n residues (20 values; n = median SP length of the TRAINING set)
#   * hydrophobicity (Kyte-Doolittle, window 5, first 40 residues): max, mean, position of the peak
#   * charge of the first 10 residues (K,R = +1; D,E = -1)
# The cleavage site is used ONLY to compute n from the training positives; it is never a feature.
AA = "ACDEFGHIKLMNPQRSTVWY"
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
FEATURE_NAMES = [f"aac_{a}" for a in AA] + ["hydro_max", "hydro_mean", "hydro_peak_pos", "charge_1_10"]
D = len(FEATURE_NAMES)

def aac(seq, n):
    """composition of the first n residues (20-D, sums to 1)"""
    pre = seq[:n]
    c = np.array([pre.count(a) for a in AA], dtype=float)
    return c / c.sum() if c.sum() > 0 else c

def hydro_profile(seq, w=5, L=40):
    """Kyte-Doolittle profile with a sliding window w on the first L residues"""
    v = np.array([KD.get(a, 0.0) for a in seq[:L]])
    if len(v) < w:
        return np.array([v.mean() if len(v) else 0.0])
    return np.convolve(v, np.ones(w) / w, mode="valid")

def static_features(seq):
    p = hydro_profile(seq)
    head = seq[:10]
    charge = sum(a in "KR" for a in head) - sum(a in "DE" for a in head)
    return [p.max(), p.mean(), float(p.argmax()), charge]

def build_X(df, n):
    A = np.array([aac(s, n) for s in df.seq])
    S = np.array([static_features(s) for s in df.seq])
    return np.hstack([A, S])

def get_n(df):
    """n = median SP length of the positives of THIS training set"""
    return int(round(df.loc[df.y == 1, "sp_len"].median()))

n_all = get_n(train)
print(f"\nD = {D} features | n (median SP length of the training set) = {n_all}")

# Look at the range of a few features: they are measured on completely different scales.
X_all = build_X(train, n_all)
for name in ["aac_L", "hydro_max", "hydro_peak_pos", "charge_1_10"]:
    j = FEATURE_NAMES.index(name)
    print(f"{name:15s} min = {X_all[:, j].min():7.3f}   max = {X_all[:, j].max():7.3f}")

# Positives vs negatives on the extra features (also a figure for the report)
fig, axes = plt.subplots(1, 4, figsize=(15, 3.2))
for ax, f in zip(axes, ["hydro_max", "hydro_mean", "hydro_peak_pos", "charge_1_10"]):
    j = FEATURE_NAMES.index(f)
    ax.boxplot([X_all[train.y.values == 0, j], X_all[train.y.values == 1, j]], showfliers=False)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["no SP", "SP"])
    ax.set_title(f)
finish_plot("svm_feature_boxplots")

# For sections 2-6 we keep the benchmark set UNTOUCHED and, like in the class notebook, we split
# the training data: CV fold 0 is the test set, folds 1-4 are the training set
# (the folds are already stratified: same positive / negative ratio everywhere).
tr, te = train[train.fold != 0], train[train.fold == 0]
n = get_n(tr)                                   # n from the training part only
X_train, y_train = build_X(tr, n), tr.y.values
X_test, y_test = build_X(te, n), te.y.values
print("\ntraining:", X_train.shape, np.bincount(y_train), "| test:", X_test.shape, np.bincount(y_test))

# A small helper to print the scores (the positive class is 1 = signal peptide)
def evaluate(y_true, y_pred):
    print(f"accuracy  : {accuracy_score(y_true, y_pred):.3f}")
    print(f"MCC       : {matthews_corrcoef(y_true, y_pred):.3f}")
    print(f"precision : {precision_score(y_true, y_pred, zero_division=0):.3f}   (positive class = 1)")
    print(f"recall    : {recall_score(y_true, y_pred, zero_division=0):.3f}")
    print(f"F1        : {f1_score(y_true, y_pred, zero_division=0):.3f}")

def metrics(y, pred):
    return {"MCC": matthews_corrcoef(y, pred),
            "Precision": precision_score(y, pred, zero_division=0),
            "Recall": recall_score(y, pred, zero_division=0),
            "F1": f1_score(y, pred, zero_division=0),
            "Accuracy": accuracy_score(y, pred)}

# %% 2. A first SVM, on the raw features
# SVC uses the RBF kernel by default (kernel='rbf'); we set it explicitly to avoid confusion.
# class_weight='balanced': only ~11% of the proteins are positives (see the lecture slide on SVC parameters).
svc_raw = SVC(C=4, kernel="rbf", gamma="scale", class_weight="balanced")
svc_raw.fit(X_train, y_train)

y_pred = svc_raw.predict(X_test)
print("\n=== SVM on the raw features ===")
evaluate(y_test, y_pred)

# %% 3. Standardizing the features
# The RBF kernel compares two samples through their Euclidean distance: ||x - x'||^2.
# If one feature ranges over [0, 20] (position of the hydrophobic peak) and another over [0, 0.3]
# (a composition), the large one decides everything.  Standardization rescales every feature to
# mean 0 and std 1:  z = (x - mean) / std.
# IMPORTANT: mean and std are learned on the TRAINING set only, then applied to the test set.
scaler = StandardScaler()
scaler.fit(X_train)                                  # learns mean and std on the training data only

X_train_s = scaler.transform(X_train)
X_test_s = scaler.transform(X_test)                  # same mean and std as the training set

print("\ntraining, before scaling: mean", X_train.mean(axis=0)[:3].round(2), "std", X_train.std(axis=0)[:3].round(2))
print("training, after scaling:  mean", X_train_s.mean(axis=0)[:3].round(2), "std", X_train_s.std(axis=0)[:3].round(2))
print("test,     after scaling:  mean", X_test_s.mean(axis=0)[:3].round(2), "std", X_test_s.std(axis=0)[:3].round(2))
print("(the test set is close to 0 / 1, but not exactly: it was scaled with the training statistics)")

svc = SVC(C=4, kernel="rbf", gamma="scale", class_weight="balanced")
svc.fit(X_train_s, y_train)
y_pred_s = svc.predict(X_test_s)

print("\n=== SVM on the standardized features ===")
evaluate(y_test, y_pred_s)
print(f"\nMCC without scaling: {matthews_corrcoef(y_test, y_pred):.3f}")
print(f"MCC with scaling:    {matthews_corrcoef(y_test, y_pred_s):.3f}")
print("With RBF (and polynomial) kernels, always standardize the features: fit the scaler on the training data only.")

# %% 4. Looking at the trained model
# 4.1 Support vectors
print("\nsupport vectors:", svc.support_vectors_.shape)              # (number of support vectors, number of features)
print("per class:", svc.n_support_, "  (class 0 = no SP, class 1 = SP)")
print("fraction of the training set:", round(len(svc.support_vectors_) / len(X_train_s), 2))

# 4.2 The value of gamma
# We did not set gamma, so scikit-learn used gamma='scale':
#       gamma = 1 / (n_features * var(X_train))
n_features = X_train_s.shape[1]
gamma_by_hand = 1.0 / (n_features * X_train_s.var())
print("\ngamma used by the SVC:  ", svc._gamma)              # _gamma stores the numerical value actually used
print("gamma computed by hand: ", gamma_by_hand)
print("(after standardization the variance is 1, so gamma = 1/n_features =", round(1 / n_features, 4), ")")

# 4.3 The decision function, rebuilt by hand
# predict() returns a class; decision_function(x) returns the value f(x) whose SIGN gives the class:
#       f(x) = sum_{s in SV}  alpha_s y_s  K(x_s, x)  +  b
# After training, scikit-learn stores the ingredients:
#       support_vectors_  the support vectors x_s
#       dual_coef_        the products alpha_s * y_s
#       intercept_        the constant b
def rbf(A, B, gamma):
    """K(a, b) = exp(-gamma ||a - b||^2) for all pairs (rows of A) x (rows of B)"""
    d2 = (A ** 2).sum(1)[:, None] + (B ** 2).sum(1)[None, :] - 2 * A @ B.T
    return np.exp(-gamma * d2)

x_few = X_test_s[:200]
f_hand = rbf(x_few, svc.support_vectors_, svc._gamma) @ svc.dual_coef_[0] + svc.intercept_[0]
f_sklearn = svc.decision_function(x_few)

print("\nfirst 5 values   by hand:", f_hand[:5].round(4))
print("first 5 values sklearn:  ", f_sklearn[:5].round(4))
print("identical for all 200 points:", np.allclose(f_hand, f_sklearn))
print("sign of f(x) gives the prediction:", np.array_equal((f_sklearn > 0).astype(int), svc.predict(x_few)))
print("Note: f(x) is proportional to the distance from the separating surface in the feature space,")
print("      but it is not the distance itself, and it is NOT a probability.")

# %% 5. Probabilities with predict_proba
# An SVM does not produce probabilities. With probability=True scikit-learn fits an extra calibration
# step (Platt scaling) that maps f(x) to a probability, using an internal 5-fold CV:
#   - training is slower;  - the internal CV is random -> set random_state;
#   - the probabilities come from a separate fit, so predict_proba > 0.5 can DISAGREE with predict().
svc_prob = SVC(C=4, kernel="rbf", gamma="scale", class_weight="balanced", probability=True, random_state=0)
svc_prob.fit(X_train_s, y_train)

proba = svc_prob.predict_proba(X_test_s)                 # one column per class: P(class 0), P(class 1)
print("\nfirst 5 rows of predict_proba:\n", proba[:5].round(3))
print("rows sum to 1:", np.allclose(proba.sum(axis=1), 1))

pred_predict = svc_prob.predict(X_test_s)
pred_proba = (proba[:, 1] > 0.5).astype(int)
disagree = np.where(pred_predict != pred_proba)[0]
print(f"test points where predict() and predict_proba() > 0.5 disagree: {len(disagree)} of {len(y_test)}")
print("(they are borderline points; for the project use predict() / decision_function(), not the probabilities)")

# %% 6. What happens with a badly chosen gamma
# gamma sets the width of the Gaussian: each support vector influences only a tiny region when gamma is large.
svc_bad = SVC(C=4, kernel="rbf", gamma=1, class_weight="balanced")
svc_bad.fit(X_train_s, y_train)
y_pred_bad = svc_bad.predict(X_test_s)
print("\n=== gamma = 1 (on 24 standardized features) ===")
print("predicted classes:", dict(zip(*np.unique(y_pred_bad, return_counts=True))))
print("support vectors:", len(svc_bad.support_vectors_), "of", len(X_train_s), "training points")
evaluate(y_test, y_pred_bad)

# Why? Look at the typical squared distance between two training points.
sample = X_train_s[np.random.RandomState(0).choice(len(X_train_s), 300, replace=False)]
d2 = ((sample[:, None, :] - sample[None, :, :]) ** 2).sum(-1)
d2 = d2[np.triu_indices(len(sample), k=1)]
print(f"\nmedian squared distance between two training points: {np.median(d2):.1f}")
for g in [1, svc._gamma]:
    print(f"gamma = {g:.4f} -> median kernel value exp(-gamma d^2) = {np.exp(-g * np.median(d2)):.2e}")
print("With gamma = 1 almost every pair has K ~ 0: every point is 'far' from every other, f(x) ~ b for any")
print("new point and the SVM predicts a single class (MCC ~ 0).  gamma must be chosen with the scale")
print("and the number of features: 'scale' (= 1/(n_features*var)) is a sensible start; then tune by CV.")

# %% 7. Choosing the hyper-parameters: feature selection + nested cross-validation
# Pipeline of the project (lectures 04-SVM and 05-Benchmarking):
#   sequences -> features (n from the training set) -> Random-Forest ranking -> top-k features
#             -> StandardScaler -> RBF SVM (C, gamma)
# k, C, gamma are chosen by cross-validation on the TRAINING folds only (metric: MCC).
def rf_ranking(X, y):
    """feature indices sorted by RF importance (parameters fixed). One ranking serves every k."""
    rf = RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=N_JOBS).fit(X, y)
    return np.argsort(rf.feature_importances_)[::-1]

def fit_model(df, k, C, gamma):
    """Fit the WHOLE pipeline on df (training data only): n -> features -> RF top-k -> scaler -> SVM."""
    n_ = get_n(df)
    X = build_X(df, n_)
    top = rf_ranking(X, df.y.values)[:k]
    scaler_ = StandardScaler().fit(X[:, top])
    svm_ = SVC(kernel="rbf", C=C, gamma=gamma, class_weight="balanced").fit(scaler_.transform(X[:, top]), df.y.values)
    return {"n": n_, "top": top, "scaler": scaler_, "svm": svm_, "k": k, "C": C, "gamma": gamma}

def predict_model(model, df):
    X = build_X(df, model["n"])[:, model["top"]]
    Z = model["scaler"].transform(X)
    return model["svm"].predict(Z), model["svm"].decision_function(Z)

def make_grid(D):
    ks = sorted({max(1, D // 4), D // 2, D})                                # k in {D/4, D/2, D}
    return [(k, C, m / k) for k in ks for C in (1, 4, 16) for m in (0.5, 1, 2)]   # gamma in {0.5/k, 1/k, 2/k}

GRID = make_grid(D)

def _score_combo(Xf, yf, Xv, yv, ranking, k, C, gamma):
    top = ranking[:k]
    sc = StandardScaler().fit(Xf[:, top])
    svm_ = SVC(kernel="rbf", C=C, gamma=gamma, class_weight="balanced").fit(sc.transform(Xf[:, top]), yf)
    return matthews_corrcoef(yv, svm_.predict(sc.transform(Xv[:, top])))

def cv_grid_search(df, grid=GRID):
    """CV over the predefined folds (column `fold`): mean validation MCC for every (k, C, gamma).
    n and the RF ranking are recomputed on the training part of EVERY fold."""
    scores = {g: [] for g in grid}
    for f in sorted(df.fold.unique()):
        fit, val = df[df.fold != f], df[df.fold == f]
        n_ = get_n(fit)
        Xf, Xv = build_X(fit, n_), build_X(val, n_)
        ranking = rf_ranking(Xf, fit.y.values)                         # one RF per fold, reused for all k
        res = Parallel(n_jobs=N_JOBS)(
            delayed(_score_combo)(Xf, fit.y.values, Xv, val.y.values, ranking, *g) for g in grid)
        for g, r in zip(grid, res):
            scores[g].append(r)
    out = pd.DataFrame([{"k": k, "C": C, "gamma": g, "mean_MCC": np.mean(v)} for (k, C, g), v in scores.items()])
    return out.sort_values("mean_MCC", ascending=False).reset_index(drop=True)

def mean_se(rows):
    """per-fold metrics -> mean +- standard error (SE = std / sqrt(n))"""
    d = pd.DataFrame(rows)
    mean, se = d.mean(), d.std(ddof=1) / np.sqrt(len(d))
    return d, pd.Series({m: f"{mean[m]:.3f} ± {se[m]:.3f}" for m in d.columns})

if RUN_NESTED_CV:
    # Outer loop: each predefined fold is the test set once (the other 4 = training).
    # Inner loop: 4-fold CV on those 4 folds chooses (k, C, gamma). The test fold is never used for a choice.
    print("\n=== Nested cross-validation (outer 5 x inner 4) ===")
    outer_rows, chosen = [], []
    for i in sorted(train.fold.unique()):
        te_df, tr_df = train[train.fold == i], train[train.fold != i]
        inner = cv_grid_search(tr_df)
        best = inner.iloc[0]
        model = fit_model(tr_df, int(best.k), float(best.C), float(best.gamma))   # refit on the 4 training folds
        pred, _ = predict_model(model, te_df)                                      # the ONLY use of fold i
        outer_rows.append(metrics(te_df.y, pred))
        chosen.append({"outer_fold": int(i), "n": model["n"], "k": int(best.k), "C": float(best.C),
                       "gamma": round(float(best.gamma), 4), "inner_MCC": round(float(best.mean_MCC), 3),
                       "test_MCC": round(outer_rows[-1]["MCC"], 3)})
        print(chosen[-1])
    per_fold, cv_summary = mean_se(outer_rows)
    print("\nPer-fold test metrics:\n", per_fold.round(3).to_string())
    print("\nNested CV, mean ± standard error over the 5 outer folds:")
    print(cv_summary.to_string())
    pd.DataFrame(chosen).to_csv(f"{FIG_DIR}/svm_nested_cv_chosen_params.csv", index=False)
    per_fold.round(4).to_csv(f"{FIG_DIR}/svm_nested_cv_per_fold.csv", index=False)

# %% 8. Final model and benchmark
# Same data partition used for the von Heijne model. The benchmark set is NEVER used for training or for choices.
#   1. 5-fold CV on the WHOLE training set over the grid -> best (k, C, gamma)
#   2. retrain on the whole training set (n and the RF ranking are recomputed on it)
#   3. test on the benchmark set, once
print("\n=== Final model ===")
final_search = cv_grid_search(train)
print("Top 5 combinations (5-fold CV on the training set):\n", final_search.head(5).round(4).to_string())
print("If the best C / gamma sits on the EDGE of the grid, extend the grid (see the Exercises).")
final_search.round(5).to_csv(f"{FIG_DIR}/svm_final_grid_search.csv", index=False)

best = final_search.iloc[0]
final_model = fit_model(train, int(best.k), float(best.C), float(best.gamma))
print(f"\nfinal model: n={final_model['n']}, k={final_model['k']}, C={final_model['C']}, gamma={final_model['gamma']:.4f}")
print("selected features:", [FEATURE_NAMES[j] for j in final_model["top"]])
print("support vectors per class:", final_model["svm"].n_support_)

bench_pred, bench_score = predict_model(final_model, bench)
print("\n=== BENCHMARK (tested once) ===")
evaluate(bench.y.values, bench_pred)

# %% 9. Saving and reloading the model
# The model alone is not enough: new data must be scaled with the SAME scaler (and use the same n and the
# same selected features) before reaching the SVC. A Pipeline keeps scaler + SVC together:
# 'fit' fits the scaler on the training data, 'predict' applies it automatically.
pipe = Pipeline([("scaler", final_model["scaler"]), ("svc", final_model["svm"])])   # already fitted parts
to_save = {"n": final_model["n"], "top": final_model["top"], "pipeline": pipe,
           "feature_names": FEATURE_NAMES, "params": {k_: final_model[k_] for k_ in ("k", "C", "gamma")}}

with gzip.open(MODEL_FILE, "wb") as fh:            # 'with' closes the file when the block ends
    pickle.dump(to_save, fh)

with gzip.open(MODEL_FILE, "rb") as fh:
    loaded = pickle.load(fh)

X_b = build_X(bench, loaded["n"])[:, loaded["top"]]
pred_reloaded = loaded["pipeline"].predict(X_b)
print("\nsame predictions after reloading:", np.array_equal(pred_reloaded, bench_pred))
print("same support vectors per class:  ", loaded["pipeline"].named_steps["svc"].n_support_)
print(f"model saved in {MODEL_FILE}")
print("(to use it on a new protein: features with the saved n -> keep the saved columns 'top' -> pipeline.predict)")

# %% 10. Error analysis on the benchmark set
# Goal: explain WHY the SVM fails, not just count the errors.
bench = bench.copy()
bench["pred"], bench["score"] = bench_pred, bench_score
bench["outcome"] = np.select([(bench.y == 1) & (bench.pred == 1), (bench.y == 0) & (bench.pred == 1),
                              (bench.y == 1) & (bench.pred == 0)], ["TP", "FP", "FN"], default="TN")
print("\n=== Error analysis ===")
print(bench.outcome.value_counts().to_dict())

# 10.1 False positives: an N-terminal transmembrane helix looks like the hydrophobic h-region of an SP
tm = bench[TM_COL].astype(str).str.lower().isin(["true", "1", "yes"])
fp, neg = bench.outcome == "FP", bench.y == 0
print(f"\nFP with a TM helix in the first 90 residues: {tm[fp].sum()}/{fp.sum()} = {tm[fp].mean():.1%}")
print(f"TM-helix proteins are {tm[neg].mean():.1%} of the negatives; "
      f"FP rate with TM helix = {(fp & tm).sum() / max(tm[neg].sum(), 1):.1%}, "
      f"without = {(fp & ~tm).sum() / max((~tm[neg]).sum(), 1):.1%}")

# 10.2 False negatives: compare training positives, benchmark true positives and benchmark false negatives
n_fin = final_model["n"]
groups = {"train positives": build_X(train[train.y == 1], n_fin),
          "bench TP": build_X(bench[bench.outcome == "TP"], n_fin),
          "bench FN": build_X(bench[bench.outcome == "FN"], n_fin)}
sp_len = {"train positives": train.loc[train.y == 1, "sp_len"],
          "bench TP": bench.loc[bench.outcome == "TP", "sp_len"],
          "bench FN": bench.loc[bench.outcome == "FN", "sp_len"]}

fig, ax = plt.subplots(figsize=(4.5, 3.2))
ax.boxplot([sp_len[g].dropna() for g in groups], showfliers=False)
ax.set_xticks(range(1, len(groups) + 1))
ax.set_xticklabels(list(groups), rotation=15)
ax.set_ylabel("SP length")
finish_plot("svm_fn_sp_length")

show = ["hydro_max", "hydro_mean", "hydro_peak_pos", "charge_1_10"]
fig, axes = plt.subplots(1, len(show), figsize=(15, 3.2))
for ax, f in zip(axes, show):
    j = FEATURE_NAMES.index(f)
    ax.boxplot([groups[g][:, j] for g in groups], showfliers=False)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels(list(groups), rotation=15)
    ax.set_title(f)
finish_plot("svm_fn_features")

comp = pd.DataFrame({g: X[:, :20].mean(axis=0) for g, X in groups.items()}, index=list(AA))
print("\nmean N-terminal composition of the three groups:\n", comp.round(3).T.to_string())
comp.round(4).to_csv(f"{FIG_DIR}/svm_fn_composition.csv")

# 10.3 Errors per taxon
if KINGDOM_COL in bench.columns:
    rows = {}
    for kg, g in bench.groupby(KINGDOM_COL):
        rows[kg] = {"n": len(g), "positives": int((g.y == 1).sum()),
                    "FPR": (g.outcome == "FP").sum() / max((g.y == 0).sum(), 1),
                    "FNR": (g.outcome == "FN").sum() / max((g.y == 1).sum(), 1)}
    print("\nerrors per kingdom (careful: few positives in some kingdoms):\n", pd.DataFrame(rows).T.round(3).to_string())

bench.drop(columns=["seq"]).to_csv(f"{FIG_DIR}/svm_benchmark_predictions.tsv", sep="\t", index=False)
