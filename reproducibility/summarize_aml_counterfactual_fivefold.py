#!/usr/bin/env python3
"""Summarize five AML scMEDAL-RE as-if reconstruction analyses.

Inputs are the small CSV/JSON outputs from analyze_aml_counterfactual.py.
Gene effects are summarized descriptively across held-out cell folds. The same
target donors appear in every fold, so this script does not pool p-values.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd


GROUP_ORDER = {"AML": 0, "control": 1, "cellline": 2}
GROUP_COLORS = {"AML": "#D55E00", "control": "#0072B2", "cellline": "#777777"}
DEFAULT_GENES = ["SAMSN1", "TXNIP", "FTL", "CENPE", "SRGN"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--folds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    parser.add_argument("--source-label", default="mono_monolike")
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def fold_dir(input_root, fold, source_label):
    return input_root / "fold{}_test_{}".format(fold, source_label)


def load_fold(input_root, fold, source_label):
    directory = fold_dir(input_root, fold, source_label)
    stats_path = directory / "aml_vs_control_gene_stats.csv"
    donors_path = directory / "target_batch_gene_means.csv"
    metadata_path = directory / "analysis_metadata.json"
    for path in (stats_path, donors_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError("Required fold {} output is missing: {}".format(fold, path))

    stats = pd.read_csv(stats_path)
    required_stats = {"gene_index", "gene", "aml_minus_control", "abs_difference", "fdr_bh"}
    missing = required_stats.difference(stats.columns)
    if missing:
        raise ValueError("Fold {} statistics lack columns: {}".format(fold, sorted(missing)))
    if stats["gene_index"].duplicated().any():
        raise ValueError("Fold {} has duplicate gene indices".format(fold))
    stats = stats.sort_values("gene_index").reset_index(drop=True)
    if not np.array_equal(stats["gene_index"].to_numpy(), np.arange(len(stats))):
        raise ValueError("Fold {} gene indices are not consecutive from zero".format(fold))
    for column in ("aml_minus_control", "abs_difference", "fdr_bh"):
        values = pd.to_numeric(stats[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("Fold {} has nonfinite {}".format(fold, column))

    donors = pd.read_csv(donors_path)
    # Older per-fold outputs used "index": pandas dropped the named index
    # during batch_means.loc[batch_groups.index].
    if "target_batch" not in donors.columns and "index" in donors.columns:
        donors = donors.rename(columns={"index": "target_batch"})
    if not {"target_batch", "target_group"}.issubset(donors.columns):
        raise ValueError(
            "Fold {} donor table lacks target identifiers; columns begin {}".format(
                fold, donors.columns[:5].tolist()
            )
        )
    if donors["target_batch"].duplicated().any():
        raise ValueError("Fold {} has duplicate target batches".format(fold))
    expected_gene_columns = ["gene_{}".format(i) for i in range(len(stats))]
    if not set(expected_gene_columns).issubset(donors.columns):
        raise ValueError("Fold {} donor table lacks gene-index columns".format(fold))
    donors = donors.set_index("target_batch").sort_index()
    if not donors["target_group"].isin(GROUP_ORDER).all():
        raise ValueError("Fold {} has an unexpected target group".format(fold))

    metadata = json.loads(metadata_path.read_text())
    if metadata.get("fold") != fold or metadata.get("dataset") != "test":
        raise ValueError("Fold {} metadata do not match the requested test fold".format(fold))
    if metadata.get("gene_count") != len(stats):
        raise ValueError("Fold {} metadata gene count differs from CSV".format(fold))
    if metadata.get("target_batch_count") != len(donors):
        raise ValueError("Fold {} metadata target count differs from CSV".format(fold))
    return stats, donors, metadata


def make_gene_summary(stats_by_fold, folds):
    first = stats_by_fold[folds[0]]
    summary = first[["gene_index", "gene"]].copy()
    gene_ids = summary["gene_index"].to_numpy()
    gene_names = summary["gene"].astype(str).to_numpy()
    for fold in folds:
        current = stats_by_fold[fold]
        if not np.array_equal(current["gene_index"].to_numpy(), gene_ids):
            raise ValueError("Gene index mapping differs at fold {}".format(fold))
        if not np.array_equal(current["gene"].astype(str).to_numpy(), gene_names):
            raise ValueError("Gene symbols differ at fold {}".format(fold))
        summary["delta_fold{}".format(fold)] = current["aml_minus_control"].to_numpy()
        rank = current["abs_difference"].rank(method="min", ascending=False)
        summary["abs_rank_fold{}".format(fold)] = rank.astype(int).to_numpy()

    deltas = summary[["delta_fold{}".format(f) for f in folds]].to_numpy(dtype=float)
    ranks = summary[["abs_rank_fold{}".format(f) for f in folds]].to_numpy(dtype=int)
    positive = (deltas > 0).sum(axis=1)
    negative = (deltas < 0).sum(axis=1)
    summary["mean_delta"] = deltas.mean(axis=1)
    summary["sd_delta_across_folds"] = deltas.std(axis=1, ddof=1)
    summary["mean_abs_delta"] = np.abs(deltas).mean(axis=1)
    summary["positive_folds"] = positive
    summary["negative_folds"] = negative
    summary["same_direction_folds"] = np.maximum(positive, negative)
    summary["top30_fold_count"] = (ranks <= 30).sum(axis=1)
    summary["direction"] = np.where(
        positive == len(folds),
        "AML_higher_all_folds",
        np.where(negative == len(folds), "control_higher_all_folds", "mixed"),
    )
    return summary.sort_values(
        ["same_direction_folds", "top30_fold_count", "mean_abs_delta"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def make_donor_summary(stats_by_fold, donors_by_fold, folds, genes):
    gene_table = stats_by_fold[folds[0]][["gene_index", "gene"]]
    selected = []
    for gene in genes:
        matches = gene_table.loc[gene_table["gene"].astype(str) == gene, "gene_index"]
        if len(matches) != 1:
            raise ValueError("Expected one gene symbol {}, found {}".format(gene, len(matches)))
        selected.append((gene, int(matches.iloc[0])))

    first = donors_by_fold[folds[0]]
    baseline_groups = first["target_group"].sort_index()
    for fold in folds[1:]:
        current_groups = donors_by_fold[fold]["target_group"].sort_index()
        if not current_groups.equals(baseline_groups):
            raise ValueError("Target donor/group mapping differs at fold {}".format(fold))

    records = []
    for gene, index in selected:
        gene_column = "gene_{}".format(index)
        for donor, group in baseline_groups.items():
            values = []
            for fold in folds:
                value = float(donors_by_fold[fold].loc[donor, gene_column])
                if not np.isfinite(value):
                    raise ValueError("Nonfinite donor mean for {} at fold {}".format(gene, fold))
                values.append(value)
            record = {
                "gene": gene,
                "gene_index": index,
                "target_batch": donor,
                "target_group": group,
                "mean_across_folds": float(np.mean(values)),
                "sd_across_folds": float(np.std(values, ddof=1)),
            }
            record.update({"value_fold{}".format(f): v for f, v in zip(folds, values)})
            records.append(record)
    return pd.DataFrame(records)


def make_control_leaveoneout_summary(gene_summary, donors_by_fold, folds, genes):
    """Check each selected gene after omitting one control target donor."""
    first = donors_by_fold[folds[0]]
    controls = first.index[first["target_group"] == "control"].tolist()
    if len(controls) < 3:
        raise ValueError("Need at least three control targets for leave-one-out analysis")
    records = []
    for gene in genes:
        selected = gene_summary.loc[gene_summary["gene"] == gene]
        if len(selected) != 1:
            raise ValueError("Expected one stability row for gene {}".format(gene))
        index = int(selected["gene_index"].iloc[0])
        full_mean = float(selected["mean_delta"].iloc[0])
        gene_column = "gene_{}".format(index)
        for fold in folds:
            donor_table = donors_by_fold[fold]
            aml_mean = donor_table.loc[
                donor_table["target_group"] == "AML", gene_column
            ].mean()
            control_mean = donor_table.loc[
                donor_table["target_group"] == "control", gene_column
            ].mean()
            expected = float(selected["delta_fold{}".format(fold)].iloc[0])
            if not np.isclose(aml_mean - control_mean, expected, atol=1e-6):
                raise ValueError(
                    "Fold {} donor means disagree with the gene effect for {}".format(
                        fold, gene
                    )
                )
        for omitted_control in controls:
            fold_deltas = []
            for fold in folds:
                donor_table = donors_by_fold[fold]
                aml_values = donor_table.loc[
                    donor_table["target_group"] == "AML", gene_column
                ].to_numpy(dtype=float)
                control_values = donor_table.loc[
                    (donor_table["target_group"] == "control")
                    & (donor_table.index != omitted_control), gene_column
                ].to_numpy(dtype=float)
                if len(control_values) != len(controls) - 1:
                    raise ValueError(
                        "Fold {} does not contain expected controls".format(fold)
                    )
                fold_deltas.append(float(aml_values.mean() - control_values.mean()))
            mean_loo = float(np.mean(fold_deltas))
            record = {
                "gene": gene,
                "omitted_control": omitted_control,
                "full_mean_delta": full_mean,
                "leaveoneout_mean_delta": mean_loo,
                "shift_from_full_mean": mean_loo - full_mean,
                "same_direction_folds": int(
                    max(sum(value > 0 for value in fold_deltas),
                        sum(value < 0 for value in fold_deltas))
                ),
                "direction_preserved": bool(full_mean * mean_loo > 0),
            }
            record.update({
                "leaveoneout_delta_fold{}".format(fold): value
                for fold, value in zip(folds, fold_deltas)
            })
            records.append(record)
    return pd.DataFrame(records)


def plot_effect_heatmap(summary, folds, top_n, output_dir):
    import matplotlib.pyplot as plt

    top = summary.head(top_n).iloc[::-1]
    values = top[["delta_fold{}".format(f) for f in folds]].to_numpy(dtype=float)
    limit = max(float(np.nanmax(np.abs(values))), 0.01)
    labels = []
    for _, row in top.iterrows():
        suffix = " ({}/{})".format(int(row["same_direction_folds"]), len(folds))
        labels.append("{}{}".format(row["gene"], suffix))
    fig, ax = plt.subplots(figsize=(8.2, max(7.0, 0.33 * len(top) + 1.9)))
    image = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    ax.set_xticks(np.arange(len(folds)))
    ax.set_xticklabels(["Fold {}".format(f) for f in folds])
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Held-out cell fold")
    ax.set_ylabel("Gene (folds with same effect direction)")
    ax.set_title("AML scMEDAL-RE: counterfactual gene effects across five folds")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    colorbar.set_label("AML minus control target-donor mean")
    fig.text(
        0.5, 0.015,
        "Same target donors in each fold; colors summarize model-generated as-if reconstructions.",
        ha="center", fontsize=8, color="#555555",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    png = output_dir / "aml_counterfactual_fivefold_gene_effects.png"
    pdf = output_dir / "aml_counterfactual_fivefold_gene_effects.pdf"
    fig.savefig(png, dpi=250, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_donor_means(donor_summary, genes, output_dir):
    import matplotlib.pyplot as plt

    first = donor_summary.loc[donor_summary["gene"] == genes[0]]
    ordered = sorted(
        first["target_batch"],
        key=lambda donor: (
            GROUP_ORDER[first.loc[first["target_batch"] == donor, "target_group"].iloc[0]],
            donor,
        ),
    )
    fig, axes = plt.subplots(len(genes), 1, figsize=(14.5, 2.2 * len(genes) + 1.8), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, gene in zip(axes, genes):
        rows = donor_summary.loc[donor_summary["gene"] == gene].set_index("target_batch").loc[ordered]
        for group, color in GROUP_COLORS.items():
            mask = rows["target_group"].to_numpy() == group
            x = np.flatnonzero(mask)
            ax.errorbar(
                x, rows.iloc[x]["mean_across_folds"],
                yerr=rows.iloc[x]["sd_across_folds"],
                fmt="o", color=color, markersize=4, capsize=2, lw=1, label=group,
            )
        ax.set_ylabel(gene + "\nmodel value")
        ax.grid(axis="y", alpha=0.2)
        ax.axvline(11.5, color="#BBBBBB", linewidth=1)
        ax.axvline(16.5, color="#BBBBBB", linewidth=1)
    axes[0].legend(ncol=3, frameon=False, loc="upper right")
    axes[-1].set_xticks(np.arange(len(ordered)))
    axes[-1].set_xticklabels(ordered, rotation=60, ha="right", fontsize=8)
    axes[-1].set_xlabel("Target donor/batch condition")
    fig.suptitle("AML scMEDAL-RE: selected donor effects across five source-cell folds")
    fig.text(
        0.5, 0.008,
        "Dots: mean across folds; error bars: fold SD. Model-generated as-if reconstructions.",
        ha="center", fontsize=8, color="#555555",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    png = output_dir / "aml_counterfactual_fivefold_donor_means.png"
    pdf = output_dir / "aml_counterfactual_fivefold_donor_means.pdf"
    fig.savefig(png, dpi=250, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def main():
    args = parse_args()
    if len(args.folds) < 2 or len(set(args.folds)) != len(args.folds):
        raise ValueError("Provide at least two distinct folds")
    if args.top_n < 1 or not args.genes:
        raise ValueError("Provide positive --top-n and at least one --genes entry")
    folds = sorted(args.folds)
    stats_by_fold = {}
    donors_by_fold = {}
    metadata_by_fold = {}
    for fold in folds:
        stats_by_fold[fold], donors_by_fold[fold], metadata_by_fold[fold] = load_fold(
            args.input_root, fold, args.source_label
        )

    first_meta = metadata_by_fold[folds[0]]
    for fold in folds[1:]:
        current = metadata_by_fold[fold]
        for field in ("experiment", "scenario", "method", "gene_count", "source_celltypes"):
            if current.get(field) != first_meta.get(field):
                raise ValueError("Fold {} differs on {}".format(fold, field))

    gene_summary = make_gene_summary(stats_by_fold, folds)
    donor_summary = make_donor_summary(stats_by_fold, donors_by_fold, folds, args.genes)
    leaveoneout_summary = make_control_leaveoneout_summary(
        gene_summary, donors_by_fold, folds, args.genes
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    gene_path = args.output_dir / "aml_counterfactual_fivefold_gene_stability.csv"
    donor_path = args.output_dir / "aml_counterfactual_fivefold_selected_donor_means.csv"
    leaveoneout_path = (
        args.output_dir / "aml_counterfactual_fivefold_control_leaveoneout.csv"
    )
    gene_summary.to_csv(gene_path, index=False)
    donor_summary.to_csv(donor_path, index=False)
    leaveoneout_summary.to_csv(leaveoneout_path, index=False)
    manifest = {
        "folds": folds,
        "source_cell_counts": {
            str(fold): metadata_by_fold[fold].get("source_cell_count") for fold in folds
        },
        "target_donors_per_fold": first_meta.get("target_batch_count"),
        "target_group_counts": first_meta.get("target_group_counts"),
        "gene_count": len(gene_summary),
        "all_fold_same_direction_gene_count": int(
            (gene_summary["same_direction_folds"] == len(folds)).sum()
        ),
        "selected_donor_plot_genes": args.genes,
        "interpretation": "model-generated as-if expression reconstruction",
        "summary_scope": "Descriptive effect-direction stability across held-out source-cell folds; target donors recur across folds.",
    }
    (args.output_dir / "aml_counterfactual_fivefold_summary.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    if not args.no_plots:
        plot_effect_heatmap(gene_summary, folds, args.top_n, args.output_dir)
        plot_donor_means(donor_summary, args.genes, args.output_dir)

    print("FOLDS={}".format(",".join(str(f) for f in folds)))
    print("GENES={}".format(len(gene_summary)))
    print("TARGET_DONORS={}".format(first_meta.get("target_batch_count")))
    print("ALL_FOLD_SAME_DIRECTION_GENES={}".format(manifest["all_fold_same_direction_gene_count"]))
    print("TOP_STABLE_GENES")
    columns = ["gene", "mean_delta", "mean_abs_delta", "same_direction_folds", "top30_fold_count"]
    print(gene_summary[columns].head(20).round(5).to_string(index=False))
    print("GENE_STABILITY_CSV={}".format(gene_path))
    print("DONOR_MEANS_CSV={}".format(donor_path))
    print("CONTROL_LEAVEONEOUT_CSV={}".format(leaveoneout_path))
    print("SELECTED_GENE_CONTROL_LEAVEONEOUT")
    print(leaveoneout_summary[[
        "gene", "omitted_control", "full_mean_delta", "leaveoneout_mean_delta",
        "same_direction_folds", "direction_preserved",
    ]].round(5).to_string(index=False))
    print("FIVEFOLD_COUNTERFACTUAL_SUMMARY_OK")


if __name__ == "__main__":
    main()
