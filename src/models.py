"""Fit and evaluate the von Heijne cleavage-site matrix.

Run from the repository root: python -m src.models
"""

from pathlib import Path

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
from sklearn.metrics import (
    accuracy_score, f1_score, matthews_corrcoef, precision_score, recall_score,
)

from .sequence_data import AA, SWISS_PROT_PERCENT, cleavage_window, load_split

WINDOW = 15
AA_TO_IDX = {aa: i for i, aa in enumerate(AA)}
BACKGROUND = np.array([SWISS_PROT_PERCENT[aa] for aa in AA])
BACKGROUND /= BACKGROUND.sum()


def build_pswm(rows, sequences):
    counts = np.ones((WINDOW, len(AA)), dtype=float)
    n = 0
    for row in rows.loc[rows.label == 1].itertuples(index=False):
        window = cleavage_window(sequences[row.accession], row.cleavage_site)
        if any(aa not in AA_TO_IDX for aa in window):
            raise ValueError(f"Nonstandard residue in training cleavage motif: {row.accession}")
        counts[np.arange(WINDOW), [AA_TO_IDX[aa] for aa in window]] += 1
        n += 1
    if n == 0:
        raise ValueError("PSWM training requires at least one positive cleavage motif")
    assert (counts.sum(axis=1) == len(AA) + n).all()
    return np.log((counts / counts.sum(axis=1, keepdims=True)) / BACKGROUND)


def scan_sequence(seq, pswm):
    encoded = np.array([AA_TO_IDX.get(aa, -1) for aa in seq[:90]], dtype=int)
    if len(encoded) < WINDOW:
        raise ValueError("Sequence is shorter than the scoring window")
    windows = sliding_window_view(encoded, WINDOW)
    windows = windows[(windows >= 0).all(axis=1)]
    if len(windows) == 0:
        raise ValueError("Sequence has no window made up of standard residues")
    return float(pswm[np.arange(WINDOW), windows].sum(axis=1).max())


def score_rows(rows, sequences, pswm):
    return np.array([scan_sequence(sequences[acc], pswm) for acc in rows.accession])


def choose_threshold(y, scores):
    """Maximize MCC over score breakpoints; ties favor the higher cutoff."""
    y = np.asarray(y, dtype=int)
    scores = np.asarray(scores)
    order = np.argsort(-scores, kind="stable")
    ranked_y, ranked_scores = y[order], scores[order]
    idx = np.r_[np.flatnonzero(ranked_scores[:-1] != ranked_scores[1:]), len(y) - 1]
    tp = ranked_y.cumsum()[idx]
    fp = idx + 1 - tp
    fn = y.sum() - tp
    tn = len(y) - y.sum() - fp
    denom = np.sqrt((tp + fp).astype(float) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = np.divide(tp * tn - fp * fn, denom, out=np.zeros(len(idx)), where=denom != 0)
    if mcc.max() <= 0:
        return float(np.nextafter(ranked_scores[0], np.inf))
    return float(ranked_scores[idx[mcc.argmax()]])


def classification_metrics(y, pred):
    return {
        "MCC": matthews_corrcoef(y, pred),
        "Precision": precision_score(y, pred, zero_division=0),
        "Recall": recall_score(y, pred, zero_division=0),
        "F1": f1_score(y, pred, zero_division=0),
        "Accuracy": accuracy_score(y, pred),
    }


def run(output_dir=Path("results")):
    train, seqs = load_split("train")
    if set(train.cv_fold) != set(range(5)):
        raise ValueError("Expected five predefined folds numbered 0 through 4")
    outer_metrics = []
    oof_scores = pd.Series(index=train.index, dtype=float)
    for test_fold in range(5):
        val_fold = (test_fold + 1) % 5
        fit = train.loc[~train.cv_fold.isin((test_fold, val_fold))]
        val = train.loc[train.cv_fold == val_fold]
        test = train.loc[train.cv_fold == test_fold]
        matrix = build_pswm(fit, seqs)
        threshold = choose_threshold(val.label, score_rows(val, seqs, matrix))
        test_pred = score_rows(test, seqs, matrix) >= threshold
        outer_metrics.append({
            "test_fold": test_fold, "validation_fold": val_fold,
            "threshold": threshold, "training_positives": int(fit.label.sum()),
            **classification_metrics(test.label, test_pred),
        })

        # Tune one final cutoff from training-only scores: each was obtained
        # using a matrix fitted on the other four folds.
        oof_matrix = build_pswm(train.loc[train.cv_fold != test_fold], seqs)
        oof_scores.loc[test.index] = score_rows(test, seqs, oof_matrix)

    threshold = choose_threshold(train.label, oof_scores)
    final_pswm = build_pswm(train, seqs)
    benchmark, bench_seqs = load_split("benchmark")
    if set(benchmark.accession) & set(train.accession):
        raise ValueError("Benchmark accessions must be independent of training")
    benchmark["score"] = score_rows(benchmark, bench_seqs, final_pswm)
    benchmark["pred"] = (benchmark.score >= threshold).astype(int)
    benchmark["threshold"] = threshold

    output_dir.mkdir(parents=True, exist_ok=True)
    cv = pd.DataFrame(outer_metrics)
    cv.to_csv(output_dir / "von_heijne_cv.csv", index=False)
    benchmark.to_csv(output_dir / "von_heijne_benchmark_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(final_pswm, index=[*range(-13, 0), 1, 2], columns=list(AA)).to_csv(
        output_dir / "von_heijne_pswm.tsv", sep="\t", index_label="position"
    )
    print(cv.to_string(index=False))
    print(f"Final threshold: {threshold:.6f}")
    print("Benchmark:", classification_metrics(benchmark.label, benchmark.pred))


if __name__ == "__main__":
    run()
