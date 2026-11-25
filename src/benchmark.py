"""Compare the two fitted methods on identical benchmark accessions."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import (
    average_precision_score, confusion_matrix, PrecisionRecallDisplay,
    RocCurveDisplay, roc_auc_score,
)

from .models import classification_metrics

RESULTS = Path("results")


def run():
    expected = pd.read_csv("data/splits/benchmark.tsv", sep="\t").set_index("accession")
    methods = {
        "von Heijne": RESULTS / "von_heijne_benchmark_predictions.tsv",
        "SVM": RESULTS / "svm" / "svm_benchmark_predictions.tsv",
    }
    scores, groups = [], []
    fig, (roc, pr) = plt.subplots(1, 2, figsize=(10, 4))
    for method, path in methods.items():
        pred = pd.read_csv(path, sep="\t").set_index("accession")
        if not pred.index.is_unique or set(pred.index) != set(expected.index):
            raise ValueError(f"{method}: predictions must cover every benchmark accession once")
        pred = pred.reindex(expected.index)
        if not pred.label.equals(expected.label) or not pred.pred.isin((0, 1)).all():
            raise ValueError(f"{method}: benchmark labels or predictions are invalid")
        y = expected.label
        tn, fp, fn, tp = confusion_matrix(y, pred.pred, labels=[0, 1]).ravel()
        scores.append({
            "method": method,
            **classification_metrics(y, pred.pred),
            "ROC_AUC": roc_auc_score(y, pred.score),
            "average_precision": average_precision_score(y, pred.score),
            "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        })
        tm = expected.tm_helix_first_90.astype("string").str.lower()
        if not tm[y == 0].isin(("true", "false")).all():
            raise ValueError("Expected a true/false TM annotation for every benchmark negative")
        for has_tm in ("true", "false"):
            subset = pred.loc[(y == 0) & tm.eq(has_tm).fillna(False)]
            positives = int(subset.pred.sum())
            groups.append({"method": method, "early_TM_helix": has_tm,
                           "negatives": len(subset), "false_positives": positives,
                           "false_positive_rate": positives / len(subset) if len(subset) else float("nan")})
        RocCurveDisplay.from_predictions(y, pred.score, ax=roc, name=method)
        PrecisionRecallDisplay.from_predictions(y, pred.score, ax=pr, name=method)
    roc.plot([0, 1], [0, 1], "--", color="gray", alpha=.6)
    pr.axhline(expected.label.mean(), linestyle="--", color="gray", alpha=.6)
    roc.set_title("ROC on the independent benchmark")
    pr.set_title("Precision–recall on the independent benchmark")
    fig.tight_layout()
    fig.savefig(RESULTS / "benchmark_curves.png", dpi=150)
    plt.close(fig)
    for name, rows in (("benchmark_metrics", scores), ("benchmark_tm_false_positives", groups)):
        table = pd.DataFrame(rows)
        table.to_csv(RESULTS / f"{name}.tsv", sep="\t", index=False, float_format="%.6f")
        print(table.round(4).to_string(index=False))


if __name__ == "__main__":
    run()
