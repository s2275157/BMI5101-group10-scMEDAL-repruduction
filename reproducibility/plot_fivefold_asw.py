#!/usr/bin/env python3
"""Plot five-fold AML ASW summaries from the canonical metric table.

The input is produced by ``collect_metrics.py``.  The figure keeps the two
scMEDAL objectives explicit: low batch ASW is desirable for batch-invariant
methods, whereas high batch ASW is expected for scMEDAL-RE because it models
donor/batch-specific variation.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t


METHOD_ORDER = [
    "scMEDAL-FE",
    "scMEDAL-RE",
    "Harmony",
    "Scanorama",
    "scVI",
    "scANVI",
]

METHOD_COLORS = {
    "scMEDAL-FE": "#0072B2",
    "scMEDAL-RE": "#D55E00",
    "Harmony": "#7A7A7A",
    "Scanorama": "#7A7A7A",
    "scVI": "#7A7A7A",
    "scANVI": "#CC79A7",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot AML five-fold batch and cell-type ASW with 95% CIs."
    )
    parser.add_argument(
        "--metrics",
        type=Path,
        default=Path("results/metrics/test_metrics_long.csv"),
        help="Long metric table produced by collect_metrics.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/figures"),
        help="Directory for PNG, PDF, and plotted-data CSV outputs.",
    )
    parser.add_argument(
        "--folds",
        type=int,
        default=5,
        help="Number of folds used to calculate the 95%% t confidence interval.",
    )
    return parser.parse_args()


def load_asw(metrics_path: Path, folds: int) -> pd.DataFrame:
    table = pd.read_csv(metrics_path)
    required = {"method", "label", "metric", "mean", "std", "sem"}
    missing = required.difference(table.columns)
    if missing:
        raise RuntimeError(f"{metrics_path} is missing columns: {sorted(missing)}")
    if folds < 2:
        raise ValueError("--folds must be at least 2")

    asw = table[
        (table["metric"] == "silhouette")
        & table["label"].isin(["batch", "celltype"])
    ].copy()

    expected = {
        (method, label)
        for method in METHOD_ORDER
        for label in ["batch", "celltype"]
    }
    observed = set(zip(asw["method"], asw["label"]))
    if observed != expected:
        missing_pairs = sorted(expected.difference(observed))
        extra_pairs = sorted(observed.difference(expected))
        raise RuntimeError(
            "Unexpected ASW rows. "
            f"Missing method/label pairs: {missing_pairs}; extra pairs: {extra_pairs}"
        )

    t_critical = float(t.ppf(0.975, df=folds - 1))
    asw["ci95"] = t_critical * asw["sem"]
    asw["ci95_low"] = asw["mean"] - asw["ci95"]
    asw["ci95_high"] = asw["mean"] + asw["ci95"]
    asw["folds"] = folds
    return asw


def plot_panel(
    ax: plt.Axes,
    table: pd.DataFrame,
    label: str,
    title: str,
    direction_note: str,
) -> None:
    subset = table[table["label"] == label].set_index("method").loc[METHOD_ORDER]
    y = np.arange(len(METHOD_ORDER))

    ax.axvline(0, color="#B8B8B8", linewidth=1, zorder=0)
    for position, method in zip(y, METHOD_ORDER):
        row = subset.loc[method]
        ax.errorbar(
            row["mean"],
            position,
            xerr=row["ci95"],
            fmt="o",
            markersize=7,
            capsize=4,
            elinewidth=1.8,
            color=METHOD_COLORS[method],
            markeredgecolor="white",
            markeredgewidth=0.8,
            zorder=3,
        )
        horizontal_offset = 0.018 if label == "batch" else 0.008
        ax.text(
            row["ci95_high"] + horizontal_offset,
            position,
            f"{row['mean']:+.3f}",
            va="center",
            ha="left",
            fontsize=9,
            color="#333333",
        )

    display_labels = [
        "scANVI†" if method == "scANVI" else method for method in METHOD_ORDER
    ]
    ax.set_yticks(y, labels=display_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Average silhouette width (ASW)")
    ax.set_title(f"{title}\n{direction_note}", loc="left", fontsize=11, pad=10)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.7, alpha=0.75)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    low = float(subset["ci95_low"].min())
    high = float(subset["ci95_high"].max())
    span = high - low
    ax.set_xlim(low - 0.08 * span, high + 0.22 * span)


def main() -> None:
    args = parse_args()
    asw = load_asw(args.metrics, args.folds)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.8))
    plot_panel(
        axes[0],
        asw,
        label="batch",
        title="Batch ASW",
        direction_note="Lower = stronger suppression; high RE = intended batch modeling",
    )
    plot_panel(
        axes[1],
        asw,
        label="celltype",
        title="Cell-type ASW",
        direction_note="Higher = stronger separability; † label-informed",
    )
    fig.suptitle("AML five-fold test-set ASW", fontsize=14, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.015,
        "Points are five-fold means; whiskers are 95% t confidence intervals. "
        "scMEDAL-RE is not a batch-correction representation.",
        ha="center",
        fontsize=9,
        color="#555555",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.94), w_pad=3.0)

    png_path = args.output_dir / "aml_asw_fivefold_95ci.png"
    pdf_path = args.output_dir / "aml_asw_fivefold_95ci.pdf"
    data_path = args.output_dir / "aml_asw_fivefold_95ci.csv"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    asw[
        [
            "method",
            "label",
            "mean",
            "std",
            "sem",
            "ci95",
            "ci95_low",
            "ci95_high",
            "folds",
        ]
    ].to_csv(data_path, index=False)

    print(f"ASW_PNG={png_path}")
    print(f"ASW_PDF={pdf_path}")
    print(f"ASW_DATA={data_path}")
    print("PLOT_ASW_OK")


if __name__ == "__main__":
    main()
