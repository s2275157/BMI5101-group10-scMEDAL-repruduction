#!/usr/bin/env python3
"""Render presentation figures from existing AML coordinates and result tables.

This script performs no training, UMAP refitting, or new hypothesis tests.
Original coordinate CSVs remain private and are not copied to the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.analysis_root.resolve(), args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 16,
                         "axes.spines.top": False, "axes.spines.right": False})
    coords = root / "AML_six_methods_train_split1/umap_19batches_seed_5/umap_csv"
    cf = root / "AML_counterfactual_fivefold_mono_monolike"
    sources = []
    frames = {}
    for method in ("scMEDAL-FE", "scMEDAL-RE"):
        source = coords / f"X_umap_{method}_latent_train_1_19batchsample.csv"
        sources.append(source)
        frames[method] = pd.read_csv(source)
    fe, re = frames.values()
    assert len(fe) == len(re)
    for col in ("Cell", "celltype", "batch", "Patient_group"):
        assert fe[col].equals(re[col]), f"Unaligned source cells: {col}"
    group_colors = {"AML": "#E91E63", "cellline": "#40AD55", "control": "#E3BD14"}

    for column, stem in (("celltype", "fe_re_celltype"), ("Patient_group", "fe_re_patient_group")):
        labels = sorted(fe[column].unique())
        palette = list(plt.get_cmap("tab20").colors) + ["#242424"]
        colors = group_colors if column == "Patient_group" else dict(zip(labels, palette))
        fig = plt.figure(figsize=(15, 6.0))
        bottom = .23 if column == "celltype" else .12
        for j, (method, df) in enumerate(frames.items()):
            ax = fig.add_axes([.025 + j*.495, bottom, .45, .70 if column == "celltype" else .80])
            xy = df.iloc[:, :2].to_numpy()
            for label in labels:
                mask = (df[column] == label).to_numpy()
                ax.scatter(xy[mask, 0], xy[mask, 1], c=[colors[label]], s=4, alpha=.70,
                           edgecolors="none", rasterized=True)
            ax.set_title(method, fontsize=23, weight="bold", color="#0A3767", pad=9)
            ax.set_axis_off()
        handles = [Line2D([], [], marker="o", linestyle="", color=colors[label],
                          markersize=9, label=label) for label in labels]
        fig.legend(handles=handles, loc="lower center", ncol=7 if column == "celltype" else 3,
                   frameon=False, fontsize=15 if column == "celltype" else 18,
                   columnspacing=1.0, handletextpad=.3, bbox_to_anchor=(.50,.01))
        fig.savefig(out / f"{stem}.png", dpi=160, facecolor="white")
        fig.savefig(out / f"{stem}.pdf", facecolor="white")
        plt.close(fig)

    stability_file = cf / "aml_counterfactual_fivefold_gene_stability.csv"
    donor_file = cf / "aml_counterfactual_fivefold_selected_donor_means.csv"
    loo_file = cf / "aml_counterfactual_fivefold_control_leaveoneout.csv"
    sources.extend([stability_file, donor_file, loo_file])
    stability = pd.read_csv(stability_file).set_index("gene")
    genes = ["SAMSN1", "TXNIP", "FTL", "HIST1H1D", "CENPE", "SRGN"]
    selected = stability.loc[genes]
    selected.reset_index().to_csv(out / "selected_gene_stability.csv", index=False)
    matrix = selected[[f"delta_fold{i}" for i in range(1,6)]].to_numpy()
    lim = max(abs(matrix.min()), abs(matrix.max()))
    fig, ax = plt.subplots(figsize=(10.8, 5.1))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(5), [f"Fold {i}" for i in range(1,6)], fontsize=17)
    ax.set_yticks(range(len(genes)), [f"{g} ({int(selected.loc[g,'same_direction_folds'])}/5)" for g in genes], fontsize=17)
    for row in range(len(genes)):
        for col in range(5):
            ax.text(col, row, f"{matrix[row,col]:+.3f}", ha="center", va="center",
                    color="white" if abs(matrix[row,col]) > lim*.62 else "#172B3E", fontsize=15)
    cb = fig.colorbar(image, ax=ax, fraction=.04, pad=.03)
    cb.set_label("AML minus control target-donor mean", fontsize=15)
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(out / "selected_gene_fivefold_heatmap.png", dpi=160, facecolor="white")
    fig.savefig(out / "selected_gene_fivefold_heatmap.pdf", facecolor="white")
    plt.close(fig)

    donors = pd.read_csv(donor_file)
    d = donors[donors.gene == "SAMSN1"].copy()
    assert len(d) == 19
    d.to_csv(out / "samsn1_donor_summary.csv", index=False)
    fig, ax = plt.subplots(figsize=(15, 4.6))
    for i, row in enumerate(d.itertuples()):
        ax.errorbar(i, row.mean_across_folds, yerr=row.sd_across_folds, fmt="o",
                    color=group_colors[row.target_group], capsize=4, markersize=8)
    ax.axvline(11.5, color="#AFBFCB", linewidth=1)
    ax.axvline(16.5, color="#AFBFCB", linewidth=1)
    ax.set_xticks(range(19), d.target_batch, rotation=48, ha="right", fontsize=15)
    ax.set_ylabel("SAMSN1 predicted expression\n(model-scaled value)", fontsize=16)
    ax.set_ylim(-.03,.33)
    ax.grid(axis="y", color="#E7EDF2")
    ax.set_axisbelow(True)
    handles = [Line2D([], [], marker="o", linestyle="", color=group_colors[g], markersize=8, label=g)
               for g in ("AML", "control", "cellline")]
    ax.legend(handles=handles, loc="upper left", ncol=3, frameon=False, fontsize=16)
    fig.tight_layout()
    fig.savefig(out / "samsn1_target_donors.png", dpi=160, facecolor="white")
    fig.savefig(out / "samsn1_target_donors.pdf", facecolor="white")
    plt.close(fig)
    loo = pd.read_csv(loo_file)
    loo[loo.gene.isin(["SAMSN1","TXNIP","FTL"])].to_csv(out / "selected_control_donor_sensitivity.csv", index=False)
    manifest = {
        "source_cell_partition": "fold 1 training partition for UMAP",
        "source_cell_count": int(len(fe)),
        "same_source_cells_verified": True,
        "umap_refitted": False,
        "counterfactual_partition": "held-out test Mono/Mono-like cells in five folds",
        "effect_definition": "mean across AML target donors minus mean across control target donors",
        "uncertainty": "SD across source-cell/model folds; target donors recur across folds",
        "interpretation": "model-generated as-if simulation",
        "selected_examples": genes,
        "sources": [{"relative_path":str(p.relative_to(root)), "sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
                    for p in sources],
    }
    (out / "figure_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps({"aligned_umap_cells":len(fe),"figures":4,"output_dir":str(out)},indent=2))


if __name__ == "__main__":
    main()
