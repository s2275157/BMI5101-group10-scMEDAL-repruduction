#!/usr/bin/env python3
"""Validate canonical scMEDAL AML reproduction outputs without loading arrays."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


PARTITIONS = ("train", "val", "test")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def validate_method(name: str, config: dict) -> tuple[dict, list[str]]:
    run_path = Path(os.path.expandvars(config["run_path"]))
    errors: list[str] = []

    if not run_path.is_dir():
        return {"method": name, "run_path": str(run_path), "exists": False}, [
            f"{name}: run directory does not exist: {run_path}"
        ]

    split_dirs = sorted(p for p in run_path.glob("splits_*") if p.is_dir())
    latent_files = sorted(
        p for split in split_dirs for p in split.rglob("*.npy") if "latent" in p.name.lower()
    )
    summary_files = sorted(run_path.glob("*_scores_*_samplesize-*.csv"))
    zero_byte_files = sorted(p for p in run_path.rglob("*") if p.is_file() and p.stat().st_size == 0)
    counterfactual_files = sorted(
        p for split in split_dirs for p in split.rglob("recon_batch_*.npy")
    )

    expected_splits = int(config.get("expected_splits", 5))
    expected_latents = int(config.get("expected_latent_files", expected_splits * 3))

    if len(split_dirs) != expected_splits:
        errors.append(f"{name}: expected {expected_splits} split directories, found {len(split_dirs)}")
    if len(latent_files) != expected_latents:
        errors.append(f"{name}: expected {expected_latents} latent arrays, found {len(latent_files)}")
    if len(summary_files) != 6:
        errors.append(f"{name}: expected 6 score summary CSVs, found {len(summary_files)}")
    if zero_byte_files:
        errors.append(f"{name}: found {len(zero_byte_files)} zero-byte files")

    for split in split_dirs:
        names = [p.name.lower() for p in latent_files if p.parent == split]
        for partition in PARTITIONS:
            if not any(partition in filename for filename in names):
                errors.append(f"{name}: {split.name} has no latent array for {partition}")

    expected_counterfactuals = config.get("expected_counterfactual_files")
    if expected_counterfactuals is not None and len(counterfactual_files) != int(expected_counterfactuals):
        errors.append(
            f"{name}: expected {expected_counterfactuals} counterfactual arrays, "
            f"found {len(counterfactual_files)}"
        )

    result = {
        "method": name,
        "run_path": str(run_path),
        "exists": True,
        "split_directories": len(split_dirs),
        "latent_arrays": len(latent_files),
        "score_csvs": len(summary_files),
        "counterfactual_arrays": len(counterfactual_files),
        "zero_byte_files": len(zero_byte_files),
        "status": "PASS" if not errors else "FAIL",
    }
    return result, errors


def main() -> int:
    args = parse_args()
    manifest = json.loads(args.manifest.read_text())
    all_errors: list[str] = []

    print("method\tstatus\tsplits\tlatents\tscore_csvs\tcounterfactuals\tempty")
    for name, config in manifest["methods"].items():
        result, errors = validate_method(name, config)
        print(
            f"{name}\t{result.get('status', 'FAIL')}\t"
            f"{result.get('split_directories', 0)}\t{result.get('latent_arrays', 0)}\t"
            f"{result.get('score_csvs', 0)}\t{result.get('counterfactual_arrays', 0)}\t"
            f"{result.get('zero_byte_files', 0)}"
        )
        all_errors.extend(errors)

    if all_errors:
        print("\nValidation errors:", file=sys.stderr)
        for error in all_errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print("\nALL_CANONICAL_RUNS_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
