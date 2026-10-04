#!/usr/bin/env python3
"""Audit scMEDAL-RE counterfactual arrays without loading them into memory.

The script checks the canonical run from ``run_manifest.json`` and verifies
that every fold and dataset split has one reconstruction per target batch.
It also confirms that array rows match the corresponding input metadata and
that array columns match the supplied HVG table.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


COUNTERFACTUAL_PATTERN = re.compile(
    r"^recon_batch_(train|val|test)_(.+)\.npy$"
)
DATASET_ORDER = {"train": 0, "val": 1, "test": 2}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit canonical AML scMEDAL-RE counterfactual outputs."
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--data-root",
        type=Path,
        required=True,
        help="Author-preprocessed AML scenario directory containing splits/.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--method", default="scMEDAL-RE")
    parser.add_argument("--expected-targets", type=int, default=19)
    parser.add_argument("--expected-genes", type=int, default=2916)
    return parser.parse_args()


def compact_values(series: pd.Series) -> str:
    return "|".join(sorted(series.dropna().astype(str).unique()))


def main() -> None:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text())
    if args.method not in manifest["methods"]:
        raise RuntimeError(f"Method {args.method!r} is absent from the manifest")

    method_config = manifest["methods"][args.method]
    run_path = Path(os.path.expandvars(method_config["run_path"])).expanduser()
    data_root = Path(os.path.expandvars(str(args.data_root))).expanduser()
    expected_folds = int(method_config["expected_splits"])
    expected_files = int(method_config["expected_counterfactual_files"])

    errors: list[str] = []
    inventory_rows: list[dict[str, object]] = []
    group_rows: list[dict[str, object]] = []
    target_sets: dict[tuple[int, str], set[str]] = {}

    if not run_path.is_dir():
        raise RuntimeError(f"Canonical run directory does not exist: {run_path}")
    if not data_root.is_dir():
        raise RuntimeError(f"AML data root does not exist: {data_root}")

    for fold in range(1, expected_folds + 1):
        split_dir = run_path / f"splits_{fold}"
        if not split_dir.is_dir():
            errors.append(f"Missing output directory: splits_{fold}")
            continue

        files = sorted(split_dir.glob("recon_batch_*.npy"))
        parsed_by_dataset: dict[str, list[tuple[Path, str]]] = {
            "train": [],
            "val": [],
            "test": [],
        }
        for path in files:
            match = COUNTERFACTUAL_PATTERN.match(path.name)
            if match is None:
                errors.append(f"Unrecognized counterfactual filename: {path.name}")
                continue
            dataset, target_batch = match.groups()
            parsed_by_dataset[dataset].append((path, target_batch))

        for dataset in ["train", "val", "test"]:
            key = (fold, dataset)
            meta_path = (
                data_root
                / "splits"
                / f"split_{fold}"
                / dataset
                / "meta.csv"
            )
            gene_path = (
                data_root
                / "splits"
                / f"split_{fold}"
                / dataset
                / "geneids.csv"
            )
            if not meta_path.is_file():
                errors.append(f"Missing metadata: split_{fold}/{dataset}/meta.csv")
                continue
            if not gene_path.is_file():
                errors.append(f"Missing genes: split_{fold}/{dataset}/geneids.csv")
                continue

            metadata = pd.read_csv(meta_path, low_memory=False)
            gene_count = len(pd.read_csv(gene_path, low_memory=False))

            required_meta = {
                "Cell",
                "Patient_group",
                "celltype",
                "batch",
                "original_index",
            }
            missing_meta = required_meta.difference(metadata.columns)
            if missing_meta:
                errors.append(
                    f"split_{fold}/{dataset} metadata is missing: "
                    f"{sorted(missing_meta)}"
                )

            parsed_files = parsed_by_dataset[dataset]
            targets = {target for _, target in parsed_files}
            target_sets[key] = targets
            if len(parsed_files) != args.expected_targets:
                errors.append(
                    f"split_{fold}/{dataset} has {len(parsed_files)} files; "
                    f"expected {args.expected_targets}"
                )
            if len(targets) != len(parsed_files):
                errors.append(f"split_{fold}/{dataset} has duplicate target batches")
            if gene_count != args.expected_genes:
                errors.append(
                    f"split_{fold}/{dataset} gene table has {gene_count} rows; "
                    f"expected {args.expected_genes}"
                )

            shapes: set[tuple[int, ...]] = set()
            for path, target_batch in parsed_files:
                if path.stat().st_size == 0:
                    errors.append(f"Empty counterfactual file: {path.name}")
                    continue
                array = np.load(path, mmap_mode="r", allow_pickle=False)
                shape = tuple(int(value) for value in array.shape)
                shapes.add(shape)
                if array.ndim != 2:
                    errors.append(f"{path.name} is {array.ndim}D; expected 2D")
                    n_cells = shape[0] if shape else 0
                    n_genes = shape[1] if len(shape) > 1 else 0
                else:
                    n_cells, n_genes = shape
                if n_cells != len(metadata):
                    errors.append(
                        f"{path.name} has {n_cells} rows but metadata has "
                        f"{len(metadata)}"
                    )
                if n_genes != gene_count:
                    errors.append(
                        f"{path.name} has {n_genes} genes but geneids.csv has "
                        f"{gene_count}"
                    )

                inventory_rows.append(
                    {
                        "fold": fold,
                        "dataset": dataset,
                        "target_batch": target_batch,
                        "filename": path.name,
                        "n_cells": n_cells,
                        "n_genes": n_genes,
                        "dtype": str(array.dtype),
                        "file_size_bytes": path.stat().st_size,
                    }
                )

            if len(shapes) > 1:
                errors.append(
                    f"split_{fold}/{dataset} counterfactual shapes disagree: "
                    f"{sorted(shapes)}"
                )

            group_rows.append(
                {
                    "fold": fold,
                    "dataset": dataset,
                    "target_batch_count": len(targets),
                    "n_cells": len(metadata),
                    "n_genes": gene_count,
                    "source_batch_count": metadata["batch"].nunique()
                    if "batch" in metadata
                    else 0,
                    "source_batches": compact_values(metadata["batch"])
                    if "batch" in metadata
                    else "",
                    "celltype_count": metadata["celltype"].nunique()
                    if "celltype" in metadata
                    else 0,
                    "celltypes": compact_values(metadata["celltype"])
                    if "celltype" in metadata
                    else "",
                    "patient_groups": compact_values(metadata["Patient_group"])
                    if "Patient_group" in metadata
                    else "",
                }
            )

    if len(inventory_rows) != expected_files:
        errors.append(
            f"Found {len(inventory_rows)} counterfactual arrays; "
            f"manifest expects {expected_files}"
        )

    if target_sets:
        reference_key = sorted(target_sets)[0]
        reference_targets = target_sets[reference_key]
        for key, targets in sorted(target_sets.items()):
            if targets != reference_targets:
                errors.append(
                    f"Target set mismatch for split_{key[0]}/{key[1]}: "
                    f"missing={sorted(reference_targets - targets)}, "
                    f"extra={sorted(targets - reference_targets)}"
                )
    else:
        reference_targets = set()

    if errors:
        raise RuntimeError("Counterfactual audit failed:\n- " + "\n- ".join(errors))

    inventory = pd.DataFrame(inventory_rows).sort_values(
        ["fold", "dataset", "target_batch"],
        key=lambda column: column.map(DATASET_ORDER)
        if column.name == "dataset"
        else column,
    )
    groups = pd.DataFrame(group_rows).sort_values(
        ["fold", "dataset"],
        key=lambda column: column.map(DATASET_ORDER)
        if column.name == "dataset"
        else column,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    inventory_path = args.output_dir / "counterfactual_inventory.csv"
    groups_path = args.output_dir / "counterfactual_groups.csv"
    summary_path = args.output_dir / "counterfactual_summary.json"
    inventory.to_csv(inventory_path, index=False)
    groups.to_csv(groups_path, index=False)

    summary = {
        "experiment": manifest.get("experiment"),
        "scenario": manifest.get("scenario"),
        "method": args.method,
        "run_name": run_path.name,
        "fold_count": expected_folds,
        "dataset_splits": ["train", "val", "test"],
        "target_batch_count": len(reference_targets),
        "target_batches": sorted(reference_targets),
        "counterfactual_file_count": len(inventory),
        "expected_gene_count": args.expected_genes,
        "audit_status": "passed",
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"COUNTERFACTUAL_FILES={len(inventory)}")
    print(f"TARGET_BATCHES={len(reference_targets)}")
    print(f"INVENTORY={inventory_path}")
    print(f"GROUPS={groups_path}")
    print(f"SUMMARY={summary_path}")
    print("COUNTERFACTUAL_AUDIT_OK")


if __name__ == "__main__":
    main()
