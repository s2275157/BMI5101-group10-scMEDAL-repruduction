#!/usr/bin/env python3
"""Collect five-fold test summaries from explicit canonical run directories."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def locate_summary(run_path: Path) -> Path:
    matches = sorted(run_path.glob("mean_scores_test_samplesize-*.csv"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one test summary in {run_path}; found {len(matches)}")
    return matches[0]


def read_summary(method: str, summary_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(summary_path, index_col=[0, 1])
    frame.index = frame.index.set_names(["label", "metric"])
    frame = frame.reset_index()

    required = {"mean", "std", "sem"}
    missing = required.difference(frame.columns)
    if missing:
        raise RuntimeError(f"{summary_path} is missing columns: {sorted(missing)}")

    frame = frame[frame["label"].isin(["batch", "celltype"])].copy()
    frame.insert(0, "method", method)
    frame.insert(1, "source_file", str(summary_path))
    return frame[["method", "source_file", "label", "metric", "mean", "std", "sem"]]


def main() -> None:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text())
    frames = []

    for method, config in manifest["methods"].items():
        run_path = Path(os.path.expandvars(config["run_path"]))
        summary_path = locate_summary(run_path)
        frames.append(read_summary(method, summary_path))

    long_table = pd.concat(frames, ignore_index=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    long_path = args.output_dir / "test_metrics_long.csv"
    long_table.to_csv(long_path, index=False)

    mean_table = long_table.pivot(index="method", columns=["label", "metric"], values="mean")
    mean_table.columns = [f"{label}_{metric}" for label, metric in mean_table.columns]
    mean_table.to_csv(args.output_dir / "test_metrics_mean.csv")

    print(f"METRICS_LONG={long_path}")
    print(f"METRICS_MEAN={args.output_dir / 'test_metrics_mean.csv'}")
    print(f"METHODS={len(frames)}")


if __name__ == "__main__":
    main()
