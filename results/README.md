# Results policy

Commit only small, reviewable results such as aggregated CSV/JSON tables and selected figures.

Do not commit:

- AML expression matrices or annotations redistributed from the authors.
- `.npy`, `.h5ad`, checkpoints, or trained weights.
- PBS output logs containing unnecessary environment or host information.
- credentials, SSH keys, tokens, or private paths beyond the documented reproducibility manifest.

Canonical large outputs remain on NSCC scratch and must be backed up separately according to NSCC retention policy.

## Generated small results

After validating the six canonical runs, generate the versionable metric tables
and primary ASW figure from the repository root:

```bash
python reproducibility/collect_metrics.py \
  --manifest reproducibility/run_manifest.json \
  --output-dir results/metrics

python reproducibility/plot_fivefold_asw.py \
  --metrics results/metrics/test_metrics_long.csv \
  --output-dir results/figures \
  --folds 5
```

The ASW plot reports means on the five supplied test subsets and nominal 95%
Student-t intervals across overlapping folds. Harmony/Scanorama integration
includes test expression. More negative raw batch ASW does not by itself imply
better integration. High batch ASW is expected for the source-modeling RE
representation. See [protocol limitations](../docs/AUDIT_PROTOCOL.md).

Before producing counterfactual expression figures, run
`reproducibility/audit_counterfactual_outputs.py`. Commit only its small CSV/JSON
inventory and summaries; do not commit the 285 reconstructed `.npy` arrays.

The curated AML counterfactual results are in
[`counterfactual/fivefold_mono_monolike/`](counterfactual/fivefold_mono_monolike/).
They contain a five-fold gene-effect table, a selected-gene donor table, a
leave-one-control-donor-out table, two figures, and a JSON summary. The fold-1
directory contains one diagnostic 19-target heatmap and its top-gene table.
Each figure summarizes model-generated as-if reconstructions from fixed
Mono/Mono-like source expression. Each target label enters the complete RE
encoder and decoder; the latent vector is recomputed for each target.

## Reconciled public evidence

The complete fold-1 target-mean and all-gene test tables are now included. Run
`python reproducibility/verify_archived_results.py` from a configured analysis
environment to recalculate them and the FE/RE per-fold, direction-stability and
donor-sensitivity summaries. Source ZIP hashes are in `provenance/`.
Curated plots from the canonical seven coordinate tables are in `umap/`.
Per-cell coordinate/annotation tables and the ZIP itself remain local.
