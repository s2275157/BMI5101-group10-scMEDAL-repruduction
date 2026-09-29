#!/usr/bin/env python3
"""Generate like-for-like AML UMAPs for Input PCA and six canonical runs."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--author-repo", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--analysis-name", default="AML_six_methods")
    parser.add_argument("--split", type=int, default=1, choices=range(1, 6))
    parser.add_argument("--dataset-type", default="train", choices=("train", "val"))
    parser.add_argument("--seed", type=int, default=5)
    parser.add_argument("--neighbors", type=int, default=15)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())

    author_repo = args.author_repo.expanduser().resolve()
    if not (author_repo / "analysis" / "compare_results_umap.py").is_file():
        raise FileNotFoundError(f"Author analysis code not found under: {author_repo}")

    results_path_dict = {}
    run_names_dict = {}
    for method, config in manifest["methods"].items():
        expanded = os.path.expandvars(config["run_path"])
        if "$" in expanded:
            raise RuntimeError(
                "An environment variable in run_manifest.json was not expanded. "
                "Set SCMEDAL_FORMAL_ROOT before running."
            )
        run_path = Path(expanded)
        if not run_path.is_dir():
            raise FileNotFoundError(f"Canonical run not found for {method}: {run_path}")
        results_path_dict[method] = str(run_path)
        run_names_dict[method] = run_path.name

    input_base_path = (
        author_repo
        / "data"
        / "AML_data"
        / manifest["scenario"]
        / "splits"
    )
    if not input_base_path.is_dir():
        raise FileNotFoundError(f"AML split data not found: {input_base_path}")

    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
    sys.path.insert(0, str(author_repo))
    os.chdir(author_repo)

    from analysis.compare_results_umap import get_umap

    output_root = args.output_root.expanduser().resolve()
    print(f"METHODS={len(results_path_dict)}")
    print(f"DATASET_TYPE={args.dataset_type}")
    print(f"SPLIT={args.split}")
    print(f"OUTPUT_ROOT={output_root}")

    get_umap(
        run_names_dict=run_names_dict,
        results_path_dict=results_path_dict,
        compare_models_path=str(output_root),
        input_base_path=str(input_base_path),
        analysis_name=args.analysis_name,
        n_pca_components=50,
        n_batches=19,
        n_neighbors=args.neighbors,
        rng_seed=args.seed,
        scaling="min_max",
        models=list(results_path_dict),
        types=[args.dataset_type],
        splits=[args.split],
        batch_col="batch",
        shape_col="celltype",
        color_col="celltype",
        use_rep="X_umap",
        issparse=False,
        extra_color_cols=["Patient_group"],
    )

    umap_dir = output_root / args.analysis_name / f"umap_19batches_seed_{args.seed}"
    png_count = len(list(umap_dir.rglob("*.png")))
    csv_count = len(list(umap_dir.rglob("*.csv")))
    print(f"UMAP_DIR={umap_dir}")
    print(f"PNG_COUNT={png_count}")
    print(f"CSV_COUNT={csv_count}")
    print("SIX_METHOD_UMAP_DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
