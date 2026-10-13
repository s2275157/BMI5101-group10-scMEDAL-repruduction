# AML scMEDAL-FE weight sensitivity

Updated 13 October 2026. The four-variant pilot values below were checked against `fe_weight_pilot_all_val.csv`, and the `adv05` PBS log was checked for its completion marker. The experimenter's local analysis directory holds copies of that CSV and log; raw logs, private scratch paths, and large outputs are not redistributed here.

## Design and status

The authors' preprocessed AML data (2,916 HVGs) and scMEDAL-FE implementation were used. The reference weights are reconstruction `1000` and adversarial/generator `1`. We changed one weight at a time. The pilot ran fold 1 on validation data; the later `adv2` experiment ran five folds. A formal test-set comparison for these weight variants has not yet been supplied. The tuning decision therefore uses validation results only.

| Variant | Reconstruction weight | Adversarial weight | Pilot fold-1 validation | Formal five-fold validation |
|---|---:|---:|---|---|
| `baseline` | 1000 | 1 | Complete | Original FE run used as comparison |
| `adv2` | 1000 | 2 | Complete | Complete |
| `recon2000` | 2000 | 1 | Complete | Not run |
| `adv05` | 1000 | 0.5 | Complete | Not run |

Pilot PBS job IDs: `25757559.pbs101` (`baseline`), `25757560.pbs101` (`adv2`), `25757561.pbs101` (`recon2000`), and `25769780.pbs101` (`adv05`). The `adv05` log ended with `FE_WEIGHT_RUN_OK` and `Finished` on 12 October 2026. The pilot output root is `scMEDAL_fe_weights/pilot` in the user's private NSCC scratch space. The exact formal `adv2` job ID is not recorded here.

The author-style objective below is `celltype_silhouette - batch_silhouette`. Larger is preferred only as this predefined combined score; batch mixing and biological retention should also be examined separately. Silhouette is not a direct measure of causality or biological truth.

## Pilot: fold 1 validation

| Variant | Batch silhouette | Cell-type silhouette | Author objective | Batch CH | Cell-type CH | Batch 1/DB | Cell-type 1/DB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `baseline` | -0.0382 | 0.0551 | 0.0933 | 103.5157 | 327.2832 | 0.1778 | 0.2613 |
| `adv05` | -0.0171 | 0.0523 | 0.0694 | 87.0449 | 248.0184 | 0.1892 | 0.2537 |
| `adv2` | -0.0443 | 0.0488 | 0.0931 | 91.7964 | 314.7250 | 0.1742 | 0.2436 |
| `recon2000` | -0.0098 | 0.0521 | 0.0618 | 90.6739 | 231.0746 | 0.1922 | 0.2554 |

Each variant has one fold here, so the sample standard deviation is undefined (`NaN`), not zero.

On this matched pilot fold, `baseline` has the highest author objective (`0.0933`), closely followed by `adv2` (`0.0931`). `adv05` (`0.0694`) and `recon2000` (`0.0618`) are lower, and both have weaker cell-type CH than `baseline`. Batch metrics move in different directions, so the pilot does not establish a clear batch-mixing improvement. One fold is a screening result, not a five-fold performance estimate.

## Formal: five-fold validation

Values are mean ± standard deviation across the five validation folds, rounded as in the supplied terminal output.

| Variant | Batch silhouette | Cell-type silhouette | Author objective | Batch CH | Cell-type CH | Batch 1/DB | Cell-type 1/DB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `original_FE` | -0.0379 ± 0.0079 | 0.0535 ± 0.0029 | 0.0913 ± 0.0066 | 93.6629 ± 4.4539 | 323.1557 ± 17.1032 | 0.1776 ± 0.0033 | 0.2550 ± 0.0034 |
| `adv2` | -0.0367 ± 0.0058 | 0.0544 ± 0.0045 | 0.0911 ± 0.0046 | 93.3137 ± 6.2785 | 322.6210 ± 26.9277 | 0.1814 ± 0.0075 | 0.2555 ± 0.0093 |

`adv2` is essentially tied with the original FE run on this validation summary: the author objective is lower by 0.0002, cell-type silhouette is higher by 0.0009, and supporting batch/cell-type metrics are mixed. This does not establish an improvement. The original FE run and the new tuning jobs were separate trainings; their random-seed comparability has not been audited. The table alone cannot support a statistical significance claim.

## Decision and reproducibility

Keep the original FE configuration as the primary benchmark. The `adv2` five-fold validation result is essentially tied with the original FE run; the single-fold `adv05` and `recon2000` pilots do not justify advancing either variant to five-fold training in this bounded study. This is a pragmatic screening decision, not proof that those weights are inferior on every fold or dataset. No test result was used to select weights.

The workflow remains available in [`fe_weight_train.py`](../reproducibility/fe_weight_train.py), [`aml_fe_weight.pbs`](../reproducibility/aml_fe_weight.pbs), and [`fe_weight_compare.py`](../reproducibility/fe_weight_compare.py). The recorded `adv05` run used reconstruction weight 1000, adversarial weight 0.5, seed 5, fold 1, and validation-only evaluation.

The repository's PBS template takes the project and private paths at submission time; it does not contain a personal NSCC account or scratch path. Example commands are in the main README.
