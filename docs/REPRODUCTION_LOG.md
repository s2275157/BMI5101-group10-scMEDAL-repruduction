# AML reproduction log

## Scope and evidence standard

This project reproduces model training and evaluation on the authors' **preprocessed AML data**. It does not reproduce FASTQ alignment, count-matrix construction, or the full upstream preprocessing pipeline.

A run is accepted only when the job was submitted, the saved log finished without a fatal traceback/OOM, all five folds ran, and the expected non-empty artefacts exist. A PBS job disappearing from `qstat` is not by itself evidence of success or failure.

## Timeline

### Environment and data validation

- Formal experiments ran on NSCC ASPIRE 2A under a personal allocation.
- Private account paths and scheduler identifiers are retained in the private project record, not this public repository.
- Public scripts resolve the formal output root through `SCMEDAL_FORMAL_ROOT`.
- Five supplied AML split directories and 16 expression-matrix files were found.
- scMEDAL environment used Python 3.8.20 and the TensorFlow/TFP/TFA stack recorded in `environment/scmedal-aml-core.txt`.

### scMEDAL-FE and scMEDAL-RE

- The first formal FE and RE smoke submissions failed before training because the PBS scripts did not load CUDA and cuDNN; `nvidia-smi` saw the A100 while TensorFlow reported no GPU.
- Adding `cuda/11.8.0` and `cudnn/11-8.9.4.25` restored TensorFlow GPU visibility.
- Both formal runs completed five folds.
- RE produced 285 counterfactual reconstructions: 5 folds x 19 batches x train/validation/test.
- The FE/RE metric and diagnostic UMAP analysis completed successfully.

### Harmony and Scanorama

- Python 3.8 could not import repository code using `list[str]`; a separate Python 3.11 environment was created instead of altering the validated scMEDAL environment.
- Fold-1 smoke tests completed before the formal runs.
- Both formal five-fold CPU jobs completed successfully.
- Canonical `13-32` runs each passed five-split, 15-latent, six-score-CSV, and zero-empty-file checks.
- Earlier `13-22` duplicates were preserved but excluded from analysis.

### scVI and scANVI

- Separate Python 3.11 environment with PyTorch 2.2.0 CUDA 11.8 and scvi-tools 1.3.0.
- Repository import gate printed `SCVI_SCANVI_REPO_IMPORT_OK`.
- TensorFlow duplicate cuFFT/cuDNN/cuBLAS registration warnings on the login node were recorded as non-fatal import warnings.
- A first `qsub` attempt was rejected by the queue before a job was created; the PBS queue configuration was corrected.
- The scANVI smoke job was accepted and checked before the formal runs.
- Formal scVI and scANVI jobs both finished on 24 September 2026.
- Each canonical run passed five-split, 15-latent, six-score-CSV, and zero-empty-file checks.

## Canonical runs

Exact canonical run names are stored in [`../reproducibility/run_manifest.json`](../reproducibility/run_manifest.json). The private NSCC root and scheduler identifiers are retained in the private project record rather than this public repository. Downstream scripts must use this manifest and must not choose a run using an unconstrained glob.

## Completed versus remaining

Completed:

- Six formal five-fold representations.
- FE/RE diagnostic metrics and UMAPs.
- RE counterfactual reconstruction generation.
- Canonical run identification and output validation.

Remaining:

- Unified five-fold test metric table across all methods.
- Consistent UMAP comparison including input PCA.
- AML counterfactual/MEC visualization and donor/patient-group interpretation.
- Final report and presentation.
