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

The ASW plot reports five-fold test means and 95% Student-t confidence
intervals. Lower batch ASW denotes stronger suppression for batch-invariant
methods; high batch ASW is expected for scMEDAL-RE because it explicitly models
donor/batch-specific variation.

Before producing counterfactual expression figures, run
`reproducibility/audit_counterfactual_outputs.py`. Commit only its small CSV/JSON
inventory and summaries; do not commit the 285 reconstructed `.npy` arrays.
