"""Generate descriptive figures for both data splits and input for WebLogo."""

import argparse
import subprocess
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .sequence_data import AA, SWISS_PROT_PERCENT, cleavage_window, load_split

FIGURES = Path("figures")
LABELS = {0: "Negative", 1: "Positive"}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIGURES / name, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_split(name):
    frame, seqs = load_split(name)
    pos = frame.loc[frame.label == 1]
    title = "Training" if name == "train" else "Benchmark"
    prefix = "training" if name == "train" else "benchmarking"

    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.geomspace(frame.protein_length.min(), frame.protein_length.max(), 31)
    for label, desc in LABELS.items():
        ax.hist(frame.loc[frame.label == label, "protein_length"], bins=bins,
                density=True, histtype="step", linewidth=1.8, label=desc)
    ax.set(xscale="log", xlabel="Protein length (residues)", ylabel="Density",
           title=f"{title}: protein length by class")
    ax.legend()
    save(fig, f"{prefix}_protein_length_comparison.png")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(pos.cleavage_site, bins=np.arange(pos.cleavage_site.min() - .5,
            pos.cleavage_site.max() + 1.5), color="#277983")
    ax.set(xlabel="Signal peptide length (residues)", ylabel="Proteins",
           title=f"{title}: signal peptide lengths")
    save(fig, f"{prefix}_sp_length_distribution.png")

    # Remove strain/common-name suffixes from the UniProt organism field.
    frame["species"] = frame.organism.str.split(" (", n=1, regex=False).str[0]
    for column in ("kingdom", "species"):
        counts = frame.groupby([column, "label"]).size().unstack(fill_value=0)
        if column == "species":
            top = counts.sum(axis=1).nlargest(12).index
            rest = counts.loc[~counts.index.isin(top)].sum()
            counts = counts.loc[top].copy()
            counts.loc["All remaining species"] = rest
        ax = counts.rename(columns=LABELS).iloc[::-1].plot.barh(
            figsize=(8, 5), color=["#5a7eaa", "#cd663d"]
        )
        ax.set(xlabel="Proteins", ylabel=column.capitalize(),
               title=f"{title}: {column} counts by class")
        save(ax.figure, f"{prefix}_{column}_distribution.png")

    counts = Counter(aa for r in pos.itertuples(index=False)
                     for aa in seqs[r.accession][:int(r.cleavage_site)] if aa in AA)
    total = sum(counts.values())
    fig, ax = plt.subplots(figsize=(9, 4))
    x = np.arange(len(AA))
    ax.bar(x - .2, [100 * counts[a] / total for a in AA], .4, label="Signal peptides")
    ax.bar(x + .2, [SWISS_PROT_PERCENT[a] for a in AA], .4, label="Swiss-Prot 2026_04")
    ax.set(xticks=x, xticklabels=list(AA), xlabel="Amino acid", ylabel="Composition (%)",
           title=f"{title}: signal peptide amino acid composition")
    ax.legend()
    save(fig, f"{prefix}_sp_amino_acid_composition.png")

    return [cleavage_window(seqs[r.accession], r.cleavage_site)
            for r in pos.itertuples(index=False)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--logos", action="store_true", help="also run WebLogo (requires Ghostscript)")
    args = parser.parse_args()
    FIGURES.mkdir(exist_ok=True)
    for name, suffix in (("train", ""), ("benchmark", "_benchmark")):
        windows = plot_split(name)
        path = FIGURES / f"cleavage_windows_13_2{suffix}.txt"
        path.write_text("\n".join(windows) + "\n")
        print(f"{name}: {len(windows)} windows written to {path}")
        if args.logos:
            alignment = "".join(f">window_{i}\n{seq}\n" for i, seq in enumerate(windows, 1))
            output = FIGURES / f"cleavage_site_logo{suffix}.png"
            title = "Training" if name == "train" else "Benchmark"
            positions = ",".join(map(str, [*range(-13, 0), "+1", "+2"]))
            logo = subprocess.run([
                "weblogo", "-A", "protein", "-F", "png",
                f"--annotate={positions}", "-t", f"{title}: cleavage-site motif",
                "-c", "chemistry",
            ], input=alignment.encode(), stdout=subprocess.PIPE, check=True)
            output.write_bytes(logo.stdout)


if __name__ == "__main__":
    main()
