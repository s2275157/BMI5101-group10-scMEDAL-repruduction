#!/usr/bin/env python3
"""Compare FE weight variants on matched folds and the same data subset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase-root", type=Path, required=True)
    parser.add_argument("--variants", nargs="+", required=True)
    parser.add_argument("--subset", choices=("val", "test"), default="val")
    parser.add_argument("--reference-run", type=Path)
    parser.add_argument("--output-csv", type=Path)
    return parser.parse_args()


def read_fold_scores(run_path: Path, folds: list[int], subset: str, label: str):
    import pandas as pd

    rows = []
    for fold in folds:
        score_path = run_path / f"splits_{fold}" / f"scores_{subset}_samplesize-10000.csv"
        if not score_path.is_file():
            raise FileNotFoundError(f"{label}: missing {score_path}")
        frame = pd.read_csv(score_path, header=[0, 1], index_col=0)
        if len(frame) != 1:
            raise ValueError(f"Expected one FE latent row in {score_path}; got {len(frame)}")
        row = {"variant": label, "fold": fold}
        for group in ("batch", "celltype"):
            for metric in ("silhouette", "ch", "1/db"):
                row[f"{group}_{metric}"] = float(frame[(group, metric)].iloc[0])
        row["author_objective"] = row["celltype_silhouette"] - row["batch_silhouette"]
        rows.append(row)
    return rows


def main() -> None:
    import pandas as pd

    args = parse_args()
    all_rows = []
    expected_folds = None
    for variant in args.variants:
        info_path = args.phase_root / variant / "run_info.json"
        if not info_path.is_file():
            raise FileNotFoundError(f"Missing completed-run record: {info_path}")
        info = json.loads(info_path.read_text(encoding="utf-8"))
        folds = info["folds"]
        if expected_folds is None:
            expected_folds = folds
        elif folds != expected_folds:
            raise ValueError("Variants were trained on different folds")
        all_rows.extend(read_fold_scores(Path(info["run_path"]), folds, args.subset, variant))

    if args.reference_run is not None:
        all_rows.extend(read_fold_scores(args.reference_run, expected_folds, args.subset, "original_FE"))

    detail = pd.DataFrame(all_rows)
    if args.output_csv is not None:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        detail.to_csv(args.output_csv, index=False)
        print(f"DETAIL_CSV={args.output_csv}")

    columns = [
        "batch_silhouette", "celltype_silhouette", "author_objective",
        "batch_ch", "celltype_ch", "batch_1/db", "celltype_1/db",
    ]
    summary = detail.groupby("variant", sort=False)[columns].agg(["mean", "std"])
    print("AUTHOR_OBJECTIVE = celltype_silhouette - batch_silhouette")
    print("For batch silhouette, inspect mixing alongside batch CH and biological scores.")
    print(summary.to_string(float_format=lambda value: f"{value:.4f}"))


if __name__ == "__main__":
    main()
