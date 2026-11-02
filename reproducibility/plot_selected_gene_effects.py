#!/usr/bin/env python3
"""Render the README's six-gene RE plot from archived five-fold CSV values."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np


GENES = ("SAMSN1", "TXNIP", "FTL", "HIST1H1D", "CENPE", "SRGN")
FOLDS = (1, 2, 3, 4, 5)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = {}
    for row in rows:
        if row["gene"] in GENES:
            if row["gene"] in selected:
                raise ValueError(f"Duplicate gene: {row['gene']}")
            selected[row["gene"]] = row
    if set(selected) != set(GENES):
        raise ValueError(f"Missing genes: {sorted(set(GENES) - set(selected))}")

    values = np.array(
        [[float(selected[gene][f"delta_fold{fold}"]) for fold in FOLDS]
         for gene in GENES],
        dtype=float,
    )
    if not np.isfinite(values).all():
        raise ValueError("Non-finite gene effect")
    for gene, gene_values in zip(GENES, values):
        agreement = max(np.count_nonzero(gene_values > 0),
                        np.count_nonzero(gene_values < 0))
        if agreement != int(selected[gene]["same_direction_folds"]):
            raise ValueError(f"Direction count disagrees for {gene}")

    # A linear, zero-centred scale preserves the meaning of small differences.
    limit = max(0.16, float(np.max(np.abs(values))))
    cmap = LinearSegmentedColormap.from_list(
        "clear_diverging",
        ["#123f73", "#4695c5", "#f5f6f6", "#e77e63", "#a82743"],
    )
    fig, ax = plt.subplots(figsize=(12.6, 5.4))
    image = ax.imshow(
        values, cmap=cmap, norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        aspect="auto", interpolation="nearest",
    )
    ax.set_xticks(range(len(FOLDS)), [f"Fold {fold}" for fold in FOLDS], fontsize=15)
    labels = [f"{gene} ({selected[gene]['same_direction_folds']}/5)" for gene in GENES]
    ax.set_yticks(range(len(GENES)), labels, fontsize=15)
    ax.set_xticks(np.arange(-0.5, len(FOLDS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(GENES), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(which="major", length=0, pad=9)
    for row_index, row_values in enumerate(values):
        for col_index, value in enumerate(row_values):
            ax.text(
                col_index, row_index, f"{value:+.3f}",
                ha="center", va="center", fontsize=15, fontweight="semibold",
                color="white" if abs(value) >= 0.09 else "#14283b",
            )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.038, pad=0.025)
    colorbar.set_label("AML minus control target-donor mean", fontsize=14, labelpad=12)
    colorbar.ax.tick_params(labelsize=12)
    fig.text(
        0.48, 0.025,
        "Model-generated reconstructions; the same target donors recur across folds.",
        ha="center", fontsize=12, color="#405367",
    )
    fig.subplots_adjust(left=0.19, right=0.90, top=0.97, bottom=0.16)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
