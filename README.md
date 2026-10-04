# BMI5101 Group 10 - scMEDAL AML Reproduction

Reproduction and evaluation of **scMEDAL (single-cell Mixed Effects Deep Autoencoder Learning)** on the authors' preprocessed Acute Myeloid Leukemia (AML) scRNA-seq dataset.

This repository contains our reproducibility layer: environment notes, canonical run records, output validation, metric aggregation, and interpretation guidance. It does **not** redistribute the authors' dataset, trained weights, or large latent/reconstruction arrays.

## Project scope

We use one public dataset, AML, to reproduce and compare:

- scMEDAL-FE: batch-invariant fixed-effects representation.
- scMEDAL-RE: batch/donor-specific random-effects representation.
- Harmony.
- Scanorama.
- scVI.
- scANVI.
- Input PCA is retained as an unintegrated reference, not a trained model.

All six trained representations use the authors' five supplied folds. The AML input contains 2,916 highly variable genes, 19 donor/batches, and 21 configured cell-type categories.

## Current status (5 October 2026)

| Method | Formal five-fold run | Validation |
|---|---:|---|
| scMEDAL-FE | Complete | Log and artefacts checked |
| scMEDAL-RE | Complete | Log, five splits, and 285 counterfactual arrays checked |
| Harmony | Complete | 5 splits, 15 latent arrays, 6 score CSVs, 0 empty files |
| Scanorama | Complete | 5 splits, 15 latent arrays, 6 score CSVs, 0 empty files |
| scVI | Complete | 5 splits, 15 latent arrays, 6 score CSVs, 0 empty files |
| scANVI | Complete | 5 splits, 15 latent arrays, 6 score CSVs, 0 empty files |

Training is complete. The six-method five-fold test summaries, primary ASW
figure, and consistently configured fold-1 UMAPs have been generated. The next
stage is counterfactual/MEC visualization and donor/patient-group interpretation.

## What scMEDAL does

For a cell expression vector `x` and one-hot batch label `z`, scMEDAL learns two complementary representations:

1. **scMEDAL-FE** uses an autoencoder plus an adversarial batch classifier. The encoder is rewarded for reconstructing expression but penalized when batch can be predicted from its features. Its latent space is intended to be batch-invariant.
2. **scMEDAL-RE** is a Bayesian autoencoder conditioned on the batch label. A positively weighted batch classifier and variational/KL regularization encourage it to model batch-specific structure.

The two subnetworks are trained separately. Their latent representations may then be concatenated for downstream prediction. This is conceptually inspired by mixed-effects decomposition, but it is not a jointly estimated classical linear mixed-effects model.

### Counterfactual reconstruction

The RE decoder can hold the learned cell representation fixed while replacing `z` with another donor/batch label. The resulting reconstruction asks:

> What expression profile would this model generate for the same encoded cell if it were assigned to another donor/batch?

This is a model-based **as-if simulation**, not evidence of a causal intervention.

## Data and upstream code

- Paper: [Nature Communications (2026), DOI 10.1038/s41467-026-72666-4](https://doi.org/10.1038/s41467-026-72666-4)
- Authors' code: [DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq)
- Authors' preprocessed data and outputs: [Figshare 10.6084/m9.figshare.28414367](https://doi.org/10.6084/m9.figshare.28414367)
- AML source dataset: GEO GSE116256

Follow the upstream repository's license. This repository does not grant additional rights over the authors' code or data.

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

See [`environment/`](environment/) for the minimal recorded package sets. These files describe the validated core stack; they are not a byte-for-byte export of every transitive package.

## Reproducing the validation and metric aggregation

1. Clone the authors' repository and obtain the authors' preprocessed AML splits.
2. Train each method using the same five folds, 50-dimensional latent spaces, and sample size of 10,000 for clustering metrics.
3. Set the private formal-output root, then use the exact canonical run names in [`reproducibility/run_manifest.json`](reproducibility/run_manifest.json):

```bash
export SCMEDAL_FORMAL_ROOT=/scratch/users/<institution>/<username>/scMEDAL_formal
```

Do not commit the expanded private path to a public repository.
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

The generated tables retain only the source CSV filename, not the expanded private
filesystem path, so the small metric outputs can be reviewed before publication.
The scripts intentionally use exact canonical paths instead of selecting the newest
directory with an unconstrained wildcard.

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

The error bars are 95% Student-t confidence intervals calculated from the
five-fold SEM. scMEDAL-RE is displayed with the comparison but is interpreted
as a batch-modeling representation, not ranked as a batch-correction method.

7. Generate directly comparable UMAPs for Input PCA and all six canonical
   representations on the same AML cells (fold 1 training partition):

```bash
qsub -P <your_project> -q normal \
  -v AUTHOR_REPO_ROOT=/path/to/scMEDAL_for_scRNAseq,SCMEDAL_FORMAL_ROOT=/path/to/scMEDAL_formal,UMAP_OUTPUT_ROOT=/path/to/scMEDAL_visualizations \
  reproducibility/aml_six_method_umap.pbs
```

The job uses the same 19 batches, seed 5, 15 neighbours, and min-max scaling
for every representation. It writes batch-, cell-type-, and patient-group-coloured
PNG files plus the UMAP-coordinate CSVs. These UMAPs are qualitative; the formal
five-fold test summaries remain the quantitative comparison.

8. Audit the scMEDAL-RE counterfactual arrays before biological interpretation:

```bash
python reproducibility/audit_counterfactual_outputs.py \
  --manifest reproducibility/run_manifest.json \
  --data-root /path/to/scMEDAL_for_scRNAseq/data/AML_data/log_transformed_2916hvggenes \
  --output-dir results/counterfactual \
  --expected-targets 19 \
  --expected-genes 2916
```

This reads only NumPy headers (`mmap_mode="r"`) and writes a public-safe file
inventory, a fold/split summary, and a JSON audit record. It verifies that each
counterfactual array has one row per source cell and one column per HVG; it does
not redistribute reconstructed expression arrays.

## How to interpret the comparison

The evaluation has two different objectives:

- For batch-correction/batch-invariant methods, lower batch clustering together with preserved cell-type clustering is desirable.
- For scMEDAL-RE, strong batch/donor separation is expected because its purpose is to model batch-specific variation.

ASW is the primary metric in the paper. CH and reciprocal DB (`1/DB`) are supporting metrics. UMAP is qualitative; independently fitted UMAP axes, rotations, and distances are not directly comparable between panels.

Do not rank FE and RE as if they solve the same task.

## Repository map

```text
.
├── README.md
├── docs/
│   ├── MODEL_EXPLAINER_ZH.md
│   └── REPRODUCTION_LOG.md
├── environment/
│   ├── scmedal-aml-core.txt
│   ├── comparables-core.txt
│   └── scvi-core.txt
├── reproducibility/
│   ├── run_manifest.json
│   ├── verify_outputs.py
│   ├── collect_metrics.py
│   ├── plot_fivefold_asw.py
│   ├── audit_counterfactual_outputs.py
│   ├── plot_six_method_umap.py
│   └── aml_six_method_umap.pbs
└── results/
    └── README.md
```

## Remaining work

- Produce AML counterfactual/MEC visualizations and interpret donor/patient-group effects.
- Prepare the final research report and presentation.
