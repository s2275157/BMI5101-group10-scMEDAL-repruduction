#!/usr/bin/env python3
"""Run one isolated AML scMEDAL-FE weight experiment on the authors' code."""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
from pathlib import Path


VARIANTS = {
    "baseline": (1000, 1),
    "adv05": (1000, 0.5),
    "adv2": (1000, 2),
    "recon2000": (2000, 1),
}
SCENARIO = "log_transformed_2916hvggenes"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--author-repo", type=Path, required=True)
    parser.add_argument("--output-base", type=Path, required=True)
    parser.add_argument("--phase", choices=("pilot", "formal"), required=True)
    parser.add_argument("--variant", choices=tuple(VARIANTS), required=True)
    parser.add_argument("--seed", type=int, default=5)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    author_repo = args.author_repo.expanduser().resolve()
    output_root = args.output_base.expanduser().resolve() / args.phase / args.variant
    recon_weight, adv_weight = VARIANTS[args.variant]
    folds = [1] if args.phase == "pilot" else [1, 2, 3, 4, 5]
    evaluate_test = args.phase == "formal"

    if not (author_repo / "models" / "models.py").is_file():
        raise FileNotFoundError(f"Not the authors' repository: {author_repo}")
    if not (author_repo / "data" / "AML_data" / SCENARIO / "splits" / "split_1").is_dir():
        raise FileNotFoundError("Preprocessed AML split_1 is missing under the authors' repository")

    print(f"PHASE={args.phase}", flush=True)
    print(f"VARIANT={args.variant}", flush=True)
    print(f"FOLDS={folds}", flush=True)
    print(f"WEIGHTS=reconstruction:{recon_weight},adversarial:{adv_weight}", flush=True)
    print(f"TEST_EVALUATION={evaluate_test}", flush=True)
    print(f"OUTPUT_ROOT={output_root}", flush=True)
    if args.dry_run:
        return

    manifest_path = output_root / "run_info.json"
    if manifest_path.exists():
        raise FileExistsError(f"A completed run is already recorded at {manifest_path}")
    output_root.mkdir(parents=True, exist_ok=True)

    # A fresh PBS process is used for every variant. This aligns Python, NumPy
    # and TensorFlow seeds, although GPU kernels may still be nondeterministic.
    os.chdir(author_repo)
    sys.path.insert(0, str(author_repo))
    random.seed(args.seed)
    import numpy as np
    import tensorflow as tf

    np.random.seed(args.seed)
    tf.random.set_seed(args.seed)
    from models.models import train_model_on_named_experiment

    model = train_model_on_named_experiment(
        "scMEDAL-FE",
        "AML",
        model_kwargs={
            "n_latent_dims": 50,
            "fold_list": folds,
            "eval_test": evaluate_test,
            "loss_recon_weight": recon_weight,
            "loss_gen_weight": adv_weight,
        },
        train_kwargs={
            "quick": False,
            "save_model": True,
            "outputs_path": str(output_root),
        },
    )

    run_name = model.model_params["run_name"]
    run_path = output_root / "AML" / "latent_space" / SCENARIO / "scmedalfe" / run_name
    subset = "test" if evaluate_test else "val"
    expected = run_path / f"mean_scores_{subset}_samplesize-10000.csv"
    if not expected.is_file() or expected.stat().st_size == 0:
        raise RuntimeError(f"Training returned without the expected score file: {expected}")

    info = {
        "phase": args.phase,
        "variant": args.variant,
        "seed": args.seed,
        "folds": folds,
        "eval_test": evaluate_test,
        "loss_recon_weight": recon_weight,
        "loss_gen_weight": adv_weight,
        "run_name": run_name,
        "run_path": str(run_path),
    }
    manifest_path.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(f"RUN_PATH={run_path}", flush=True)
    print(f"RUN_INFO={manifest_path}", flush=True)
    print("FE_WEIGHT_RUN_OK", flush=True)


if __name__ == "__main__":
    main()
