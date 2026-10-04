#!/usr/bin/env python3
"""Summarize AML scMEDAL-RE counterfactual reconstructions by target group.

For a fixed set of source cells, scMEDAL-RE reconstructs expression after the
target donor/batch condition is changed.  This script averages those
model-generated as-if reconstructions within source cells, compares AML and
control target donors at the donor level, and creates a top-gene heatmap.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


GROUP_ORDER = {"AML": 0, "control": 1, "cellline": 2}
GROUP_COLORS = {"AML": "#D55E00", "control": "#0072B2", "cellline": "#7A7A7A"}
GENE_COLUMNS = ["Gene", "gene_ids", "gene_name", "gene_symbol", "symbol", "_index"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare AML and control target-donor counterfactual expression "
            "for fixed AML source-cell states."
        )
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--method", default="scMEDAL-RE")
    parser.add_argument("--fold", type=int, default=1)
    parser.add_argument(
        "--dataset", choices=["train", "val", "test"], default="test"
    )
    parser.add_argument(
        "--source-celltypes",
        nargs="+",
        default=["Mono", "Mono-like"],
        help="Source cell types held fixed across target-donor reconstructions.",
    )
    parser.add_argument("--top-genes", type=int, default=30)
    parser.add_argument("--chunk-size", type=int, default=256)
    return parser.parse_args()


def expand_path(value: str | Path) -> Path:
    return Path(os.path.expandvars(str(value))).expanduser()


def is_plain_position_index(values: pd.Series) -> bool:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.isna().any():
        return False
    return np.array_equal(numeric.to_numpy(dtype=int), np.arange(len(values)))


def load_gene_names(data_root: Path, split_dir: Path, expected: int) -> tuple[pd.Series, Path, str]:
    candidates = [data_root / "geneids.csv", split_dir / "geneids.csv"]
    diagnostics: list[str] = []
    for path in candidates:
        if not path.is_file():
            diagnostics.append(f"missing {path}")
            continue
        table = pd.read_csv(path, low_memory=False)
        if len(table) != expected:
            diagnostics.append(f"{path} has {len(table)} rows")
            continue
        for column in GENE_COLUMNS:
            if column not in table.columns:
                continue
            values = table[column].astype(str)
            if values.notna().all() and not is_plain_position_index(values):
                return values.reset_index(drop=True), path, column
        for column in table.columns:
            values = table[column].astype(str)
            if values.notna().all() and not is_plain_position_index(values):
                return values.reset_index(drop=True), path, column
        diagnostics.append(
            f"{path} columns {table.columns.tolist()} contain only positional indices"
        )
    raise RuntimeError(
        "Could not recover biological gene names. " + "; ".join(diagnostics)
    )


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def mean_selected_rows(
    array_path: Path, indices: np.ndarray, expected_genes: int, chunk_size: int
) -> np.ndarray:
    array = np.load(array_path, mmap_mode="r", allow_pickle=False)
    if array.ndim != 2 or array.shape[1] != expected_genes:
        raise RuntimeError(f"Unexpected array shape for {array_path}: {array.shape}")
    if indices.size == 0 or indices.max() >= array.shape[0]:
        raise RuntimeError(f"Source-cell indices do not align with {array_path}")
    total = np.zeros(expected_genes, dtype=np.float64)
    for start in range(0, len(indices), chunk_size):
        rows = indices[start : start + chunk_size]
        total += np.asarray(array[rows, :], dtype=np.float64).sum(axis=0)
    return total / len(indices)


def make_unique_labels(names: pd.Series) -> list[str]:
    counts: dict[str, int] = {}
    labels: list[str] = []
    for raw_name in names.astype(str):
        counts[raw_name] = counts.get(raw_name, 0) + 1
        suffix = f" [{counts[raw_name]}]" if counts[raw_name] > 1 else ""
        labels.append(raw_name + suffix)
    return labels


def plot_heatmap(
    batch_means: pd.DataFrame,
    batch_groups: pd.Series,
    top_stats: pd.DataFrame,
    output_dir: Path,
    fold: int,
    dataset: str,
    source_celltypes: list[str],
) -> tuple[Path, Path]:
    ordered_batches = sorted(
        batch_means.index,
        key=lambda batch: (GROUP_ORDER.get(batch_groups.loc[batch], 99), batch),
    )
    top_indices = top_stats["gene_index"].astype(int).tolist()
    matrix = batch_means.loc[ordered_batches, top_indices].to_numpy().T
    row_mean = matrix.mean(axis=1, keepdims=True)
    row_sd = matrix.std(axis=1, ddof=0, keepdims=True)
    row_sd[row_sd == 0] = 1.0
    z_matrix = (matrix - row_mean) / row_sd

    gene_labels = make_unique_labels(top_stats["gene"])
    group_codes = np.array(
        [[GROUP_ORDER.get(batch_groups.loc[batch], 3) for batch in ordered_batches]]
    )
    group_cmap = plt.matplotlib.colors.ListedColormap(
        [GROUP_COLORS["AML"], GROUP_COLORS["control"], GROUP_COLORS["cellline"], "#FFFFFF"]
    )

    figure_height = max(8.0, 0.30 * len(gene_labels) + 2.8)
    fig = plt.figure(figsize=(14.5, figure_height))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.18, 6.0], hspace=0.03)
    group_ax = fig.add_subplot(grid[0])
    heat_ax = fig.add_subplot(grid[1])

    group_ax.imshow(group_codes, aspect="auto", cmap=group_cmap, vmin=0, vmax=3)
    group_ax.set_xticks([])
    group_ax.set_yticks([])
    group_ax.set_ylabel("Target\ngroup", rotation=0, ha="right", va="center")
    for spine in group_ax.spines.values():
        spine.set_visible(False)

    image = heat_ax.imshow(z_matrix, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    heat_ax.set_xticks(np.arange(len(ordered_batches)), labels=ordered_batches)
    heat_ax.tick_params(axis="x", rotation=60, labelsize=9)
    heat_ax.set_yticks(np.arange(len(gene_labels)), labels=gene_labels)
    heat_ax.tick_params(axis="y", labelsize=8)
    heat_ax.set_xlabel("Target donor/batch condition")
    heat_ax.set_ylabel("Genes ranked by absolute AML-control mean difference")
    fig.suptitle(
        "scMEDAL-RE counterfactual expression\n"
        f"Fold {fold} {dataset}; fixed source cells: {' + '.join(source_celltypes)}",
        y=0.98,
        fontweight="bold",
    )
    colorbar = fig.colorbar(image, ax=heat_ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Gene-wise z-score across target donors")
    legend = [
        Patch(facecolor=color, label=group)
        for group, color in GROUP_COLORS.items()
    ]
    heat_ax.legend(
        handles=legend,
        title="Target group",
        bbox_to_anchor=(1.04, 1.0),
        loc="upper left",
        frameon=False,
    )
    fig.text(
        0.5,
        0.01,
        "Model-generated as-if reconstructions; each column is averaged over the same source cells.",
        ha="center",
        fontsize=9,
        color="#555555",
    )
    fig.subplots_adjust(left=0.12, right=0.84, bottom=0.14, top=0.84)

    png_path = output_dir / "aml_counterfactual_top_genes_heatmap.png"
    pdf_path = output_dir / "aml_counterfactual_top_genes_heatmap.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def main() -> None:
    args = parse_args()
    if args.fold < 1:
        raise ValueError("--fold must be at least 1")
    if args.top_genes < 1:
        raise ValueError("--top-genes must be at least 1")
    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be at least 1")

    manifest = json.loads(args.manifest.read_text())
    method_config = manifest["methods"].get(args.method)
    if method_config is None:
        raise RuntimeError(f"Method {args.method!r} is absent from the manifest")
    expected_folds = int(method_config["expected_splits"])
    if args.fold > expected_folds:
        raise ValueError(f"--fold exceeds the manifest fold count ({expected_folds})")

    run_path = expand_path(method_config["run_path"])
    data_root = expand_path(args.data_root)
    data_split_dir = data_root / "splits" / f"split_{args.fold}" / args.dataset
    reconstruction_dir = run_path / f"splits_{args.fold}"
    meta_path = data_split_dir / "meta.csv"

    if not run_path.is_dir():
        raise RuntimeError(f"Canonical run directory does not exist: {run_path}")
    if not meta_path.is_file():
        raise RuntimeError(f"Metadata file does not exist: {meta_path}")

    metadata = pd.read_csv(meta_path, low_memory=False)
    required_columns = {"Cell", "Patient_group", "celltype", "batch"}
    missing = required_columns.difference(metadata.columns)
    if missing:
        raise RuntimeError(f"Metadata is missing columns: {sorted(missing)}")

    source_mask = metadata["celltype"].isin(args.source_celltypes).to_numpy()
    source_indices = np.flatnonzero(source_mask)
    if source_indices.size == 0:
        available = sorted(metadata["celltype"].dropna().astype(str).unique())
        raise RuntimeError(
            f"No source cells matched {args.source_celltypes}; available={available}"
        )

    reconstruction_files = sorted(
        reconstruction_dir.glob(f"recon_batch_{args.dataset}_*.npy")
    )
    if not reconstruction_files:
        raise RuntimeError(
            f"No counterfactual arrays found in {reconstruction_dir} for {args.dataset}"
        )

    batch_group_table = metadata[["batch", "Patient_group"]].drop_duplicates()
    conflicting = batch_group_table.groupby("batch")["Patient_group"].nunique()
    if (conflicting > 1).any():
        raise RuntimeError(
            "Some batches map to multiple Patient_group values: "
            f"{conflicting[conflicting > 1].index.tolist()}"
        )
    batch_to_group = batch_group_table.set_index("batch")["Patient_group"].to_dict()

    first_array = np.load(reconstruction_files[0], mmap_mode="r", allow_pickle=False)
    if first_array.ndim != 2 or first_array.shape[0] != len(metadata):
        raise RuntimeError(
            f"Reconstruction/metadata alignment failed: {first_array.shape} vs {metadata.shape}"
        )
    n_genes = int(first_array.shape[1])
    gene_names, gene_path, gene_column = load_gene_names(
        data_root, data_split_dir, n_genes
    )

    target_means: dict[str, np.ndarray] = {}
    target_groups: dict[str, str] = {}
    prefix = f"recon_batch_{args.dataset}_"
    for path in reconstruction_files:
        if not path.stem.startswith(prefix):
            raise RuntimeError(f"Unexpected counterfactual filename: {path.name}")
        target_batch = path.stem[len(prefix) :]
        if target_batch not in batch_to_group:
            raise RuntimeError(f"Target batch {target_batch} has no Patient_group mapping")
        target_means[target_batch] = mean_selected_rows(
            path, source_indices, n_genes, args.chunk_size
        )
        target_groups[target_batch] = str(batch_to_group[target_batch])

    batch_means = pd.DataFrame.from_dict(target_means, orient="index")
    batch_means.index.name = "target_batch"
    batch_groups = pd.Series(target_groups, name="target_group")
    batch_means = batch_means.loc[batch_groups.index]
    batch_means.index.name = "target_batch"

    aml_batches = batch_groups[batch_groups.str.lower() == "aml"].index
    control_batches = batch_groups[batch_groups.str.lower() == "control"].index
    if len(aml_batches) < 2 or len(control_batches) < 2:
        raise RuntimeError(
            f"Need at least two AML and control targets; got {len(aml_batches)} and {len(control_batches)}"
        )

    aml_values = batch_means.loc[aml_batches].to_numpy()
    control_values = batch_means.loc[control_batches].to_numpy()
    statistic, p_value = mannwhitneyu(
        aml_values, control_values, axis=0, alternative="two-sided"
    )
    aml_mean = aml_values.mean(axis=0)
    control_mean = control_values.mean(axis=0)
    delta = aml_mean - control_mean

    stats = pd.DataFrame(
        {
            "gene_index": np.arange(n_genes),
            "gene": gene_names,
            "aml_target_mean": aml_mean,
            "control_target_mean": control_mean,
            "aml_minus_control": delta,
            "abs_difference": np.abs(delta),
            "mannwhitney_u": statistic,
            "p_value": p_value,
            "fdr_bh": benjamini_hochberg(p_value),
        }
    ).sort_values(["abs_difference", "p_value"], ascending=[False, True])
    top_stats = stats.head(min(args.top_genes, n_genes)).copy()

    output_dir = expand_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    batch_output = batch_means.copy()
    batch_output.insert(0, "target_group", batch_groups.loc[batch_output.index])
    batch_output.columns = ["target_group"] + [f"gene_{i}" for i in range(n_genes)]
    batch_output.reset_index().to_csv(
        output_dir / "target_batch_gene_means.csv", index=False
    )
    pd.DataFrame(
        {"gene_index": np.arange(n_genes), "gene": gene_names}
    ).to_csv(output_dir / "gene_index.csv", index=False)
    stats.to_csv(output_dir / "aml_vs_control_gene_stats.csv", index=False)
    top_stats.to_csv(output_dir / "top_genes.csv", index=False)

    source_summary = (
        metadata.loc[source_mask, ["Patient_group", "batch", "celltype"]]
        .value_counts(sort=False)
        .rename("n_source_cells")
        .reset_index()
        .sort_values(["Patient_group", "batch", "celltype"])
    )
    source_summary.to_csv(output_dir / "source_cell_summary.csv", index=False)

    png_path, pdf_path = plot_heatmap(
        batch_means,
        batch_groups,
        top_stats,
        output_dir,
        args.fold,
        args.dataset,
        args.source_celltypes,
    )

    analysis_metadata = {
        "experiment": manifest.get("experiment"),
        "scenario": manifest.get("scenario"),
        "method": args.method,
        "run_name": run_path.name,
        "fold": args.fold,
        "dataset": args.dataset,
        "source_celltypes": args.source_celltypes,
        "source_cell_count": int(source_indices.size),
        "target_batch_count": int(len(batch_means)),
        "target_group_counts": {
            str(key): int(value) for key, value in batch_groups.value_counts().items()
        },
        "gene_count": n_genes,
        "gene_table": gene_path.name,
        "gene_column": gene_column,
        "top_gene_ranking": "absolute AML-control difference across target-batch means",
        "statistical_unit": "target donor/batch mean",
        "interpretation": "model-generated as-if expression reconstruction",
    }
    metadata_path = output_dir / "analysis_metadata.json"
    metadata_path.write_text(json.dumps(analysis_metadata, indent=2) + "\n")

    print(f"SOURCE_CELLS={source_indices.size}")
    print(f"TARGET_BATCHES={len(batch_means)}")
    print(f"AML_TARGETS={len(aml_batches)}")
    print(f"CONTROL_TARGETS={len(control_batches)}")
    print(f"GENE_COLUMN={gene_column}")
    print(f"TOP_GENE={top_stats.iloc[0]['gene']}")
    print(f"HEATMAP_PNG={png_path}")
    print(f"HEATMAP_PDF={pdf_path}")
    print(f"ANALYSIS_METADATA={metadata_path}")
    print("COUNTERFACTUAL_ANALYSIS_OK")


if __name__ == "__main__":
    main()
