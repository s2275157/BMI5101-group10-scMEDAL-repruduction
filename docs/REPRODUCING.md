# Reproducing analysis and recorded runs

Run all commands from the repository root. The public-summary route below is
independent of the historical GPU training environments.

## Verify the public result bundle

Use a separate Python 3.12 environment:

```bash
python3.12 -m venv .venv-analysis
.venv-analysis/bin/python -m pip install -r environment/analysis-review-core.txt
.venv-analysis/bin/python reproducibility/verify_archived_results.py
```

This recalculates the FE/RE fold summaries, all fold-1 gene tests and BH-FDR,
five-fold gene directions and selected donor sensitivity. It requires no private
arrays, source ZIP or GPU. It does not retrain the models.

## Reconcile the supplied result ZIP

If you possess `scMEDAL_FE_RE_analysis.zip`, import its small result tables into
a staging directory and validate them:

```bash
.venv-analysis/bin/python reproducibility/reconcile_result_archive.py \
  --archive /path/to/scMEDAL_FE_RE_analysis.zip \
  --output-root /path/to/staging

.venv-analysis/bin/python reproducibility/plot_archived_umap.py \
  --archive /path/to/scMEDAL_FE_RE_analysis.zip \
  --output-dir /path/to/staging/results/umap
```

The first script ignores Finder metadata/resource forks and copies only selected
summary files. It normalizes the legacy target-table header `index` to
`target_batch`, retaining the numeric records, and reshapes the old FE/RE
multi-header CSV into 60 explicit fold measurements. It records source and
published hashes. The second script reads stored coordinates; it does not fit a
new UMAP. Per-cell source annotations remain local. Canonical coordinate family:
`AML_six_methods_train_split1`.

## Full training is a separate workflow

Training requires the upstream implementation, author-preprocessed input splits
and the historical compatible environments. The current source version reviewed
for protocol interpretation is `bbf9a0b801123d63218cf05d93ec8f103453881d`.
The exact historical training commit and local patches are not recorded in the
public run manifest. Review [the protocol limits](AUDIT_PROTOCOL.md), including
transductive integration and scVI/scANVI scaled input, before interpreting or
repeating the original pipelines. See also the
[upstream reproduction instructions](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/blob/bbf9a0b801123d63218cf05d93ec8f103453881d/docs/ExperimentsReproducibility.md).

## Compute environment

Formal experiments were run on **NSCC ASPIRE 2A** using PBS jobs.

### scMEDAL-FE / scMEDAL-RE

- Python 3.8.20
- TensorFlow 2.13.1
- TensorFlow Probability 0.21.0
- TensorFlow Addons 0.21.0
- NumPy 1.24.3
- Scanpy 1.9.8
- AnnData 0.9.2
- CUDA 11.8 and cuDNN 8.9 modules on an NVIDIA A100 node

### Harmony / Scanorama

- Separate Python 3.11 environment: `scmedal-comparables`
- CPU formal jobs

### scVI / scANVI

- Separate Python 3.11 environment: `scmedal-scvi`
- PyTorch 2.2.0 CUDA 11.8 build
- scvi-tools 1.3.0
- GPU formal jobs

See [`environment/`](../environment/) for the minimal recorded package sets. These files describe the validated core stack; they are not a byte-for-byte export of every transitive package.

## Reproducing figures from the public summaries

This route needs no private matrices or GPU training. From this repository's root,
use a Python 3.11 environment with pandas, NumPy, SciPy and Matplotlib:

```bash
python reproducibility/plot_fivefold_asw.py \
  --metrics results/metrics/test_metrics_long.csv \
  --output-dir results/figures \
  --folds 5
```

The tables retain the reported values. Changes to plot annotations do not
constitute a model rerun.

## Reproducing training and private-output aggregation

This repository supplies downstream checks and analysis wrappers. Training also
requires the upstream implementation, its data and compatible environments.
The reviewed upstream commit is `bbf9a0b801123d63218cf05d93ec8f103453881d`.
The exact upstream commit and local patches used for the historical jobs were
not recorded in the public manifest, so the reviewed commit is an audit
reference, not a verified historical training pin. See
[the protocol review](AUDIT_PROTOCOL.md) before attempting exact reproduction.

1. Clone the authors' repository and obtain the authors' preprocessed AML splits.
2. Train or integrate each method using the five supplied scoring subsets, 50-dimensional latent spaces, and the recorded metric sample-size setting of 10,000. This is an upper bound when a partition contains fewer cells. Preserve and report each method's integration and validation procedure.
3. Set the private formal-output root, then use the exact canonical run names in [`reproducibility/run_manifest.json`](../reproducibility/run_manifest.json):

```bash
export SCMEDAL_FORMAL_ROOT=/scratch/users/<institution>/<username>/scMEDAL_formal
```

4. Validate the canonical runs:

```bash
python reproducibility/verify_outputs.py \
  --manifest reproducibility/run_manifest.json
```

5. Aggregate the five-fold test summaries:

```bash
python reproducibility/collect_metrics.py \
  --manifest reproducibility/run_manifest.json \
  --output-dir results/metrics
```

Expected small outputs:

- `results/metrics/test_metrics_long.csv`
- `results/metrics/test_metrics_mean.csv`

6. Plot the primary five-fold ASW comparison with 95% confidence intervals:

```bash
python reproducibility/plot_fivefold_asw.py \
  --metrics results/metrics/test_metrics_long.csv \
  --output-dir results/figures \
  --folds 5
```

Expected small outputs:

- `results/figures/aml_asw_fivefold_95ci.png`
- `results/figures/aml_asw_fivefold_95ci.pdf`
- `results/figures/aml_asw_fivefold_95ci.csv`

The error bars are nominal 95% Student-t intervals calculated from the
five-fold SEM. Folds overlap, so these do not estimate independent donor-level
uncertainty or demonstrate a statistically significant method difference.

7. Generate consistently configured UMAPs for Input PCA and all six canonical
   representations on the same AML cells (fold 1 training partition):

```bash
qsub -P <your_project> -q normal \
  -v AUTHOR_REPO_ROOT=/path/to/scMEDAL_for_scRNAseq,SCMEDAL_FORMAL_ROOT=/path/to/scMEDAL_formal,UMAP_OUTPUT_ROOT=/path/to/scMEDAL_visualizations \
  reproducibility/aml_six_method_umap.pbs
```

The job uses the same 19 batches, seed 5, 15 neighbours, and min-max scaling
for every representation. It writes batch-, cell-type-, and patient-group-coloured
PNG files plus the UMAP-coordinate CSVs. These UMAPs are qualitative.

8. Audit the scMEDAL-RE counterfactual arrays before biological interpretation:

```bash
python reproducibility/audit_counterfactual_outputs.py \
  --manifest reproducibility/run_manifest.json \
  --data-root /path/to/scMEDAL_for_scRNAseq/data/AML_data/log_transformed_2916hvggenes \
  --output-dir results/counterfactual \
  --expected-targets 19 \
  --expected-genes 2916
```

This opens NumPy arrays in memory-mapped mode (`mmap_mode="r"`) to inspect headers and shapes and writes a public-safe file
inventory, a fold/split summary, and a JSON audit record. It verifies that each
counterfactual array has one row per source cell and one column per HVG.

9. Analyze the held-out test cells with `Mono` and `Mono-like` labels. The same
   source cells within each fold are reconstructed under all 19 target
   donor/batch conditions. On ASPIRE 2A, the supplied PBS scripts run fold 1
   and the fold 2–5 array job:

```bash
qsub -P <your_project> -q normal \
  -v AUTHOR_REPO_ROOT=/path/to/scMEDAL_for_scRNAseq,SCMEDAL_FORMAL_ROOT=/path/to/scMEDAL_formal \
  reproducibility/aml_counterfactual_analysis.pbs

qsub -P <your_project> -q normal \
  -v AUTHOR_REPO_ROOT=/path/to/scMEDAL_for_scRNAseq,SCMEDAL_FORMAL_ROOT=/path/to/scMEDAL_formal \
  reproducibility/aml_counterfactual_folds2to5.pbs
```

Each job calls [`analyze_aml_counterfactual.py`](../reproducibility/analyze_aml_counterfactual.py).
It writes a gene-level AML minus control target-donor comparison and a
19-target heatmap. The statistics use target donors as units: 12 AML and 5
controls. The two cell lines appear descriptively in the heatmap.

10. After all five fold directories exist, summarize effect direction and
    selected donor responses:

```bash
python reproducibility/summarize_aml_counterfactual_fivefold.py \
  --input-root results/counterfactual \
  --output-dir results/counterfactual/fivefold_mono_monolike
```

The summary checks gene-index and donor mappings across folds. It reports
per-fold effect sizes, the number of folds with the same direction, how often
each gene enters the fold's top 30 by absolute difference, and the effect of
omitting each control donor for five selected genes. The same target donors
recur across folds, so these are descriptive stability summaries.
