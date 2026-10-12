# AML results and presentation interpretation

Updated 12 October 2026. This guide records the evidence used in presentation version 5; the course deliverable is the PowerPoint file.

## 1. What was evaluated

The six canonical formal runs use the authors' preprocessed AML dataset and five supplied train/validation/test folds. The input has 2,916 log-transformed highly variable genes and 19 donor/batch categories. The formal test metrics come from `results/metrics/test_metrics_long.csv`; five-fold means are in `test_metrics_mean.csv`. The exact runs remain in `reproducibility/run_manifest.json`.

FE aims to learn batch-invariant structure through reconstruction and adversarial batch prediction. RE aims to learn batch-specific structure through a batch-conditioned Bayesian autoencoder. The two subnetworks were trained separately. Their names reflect a mixed-effects interpretation; this reproduction does not fit a jointly estimated classical linear mixed-effects model.

## 2. Quantitative comparison

Test ASW, rounded to three decimals, reported as five-fold mean ± fold SD:

| Representation | Batch ASW | Cell-type ASW | Interpretation |
|---|---:|---:|---|
| scMEDAL-FE | −0.038 ± 0.006 | 0.054 ± 0.003 | Batch-invariant objective |
| scMEDAL-RE | 0.685 ± 0.026 | −0.077 ± 0.007 | Batch-specific objective |
| Harmony | −0.107 ± 0.001 | 0.018 ± 0.002 | Integration |
| Scanorama | −0.170 ± 0.014 | 0.010 ± 0.011 | Integration |
| scVI | −0.236 ± 0.007 | 0.029 ± 0.012 | Probabilistic integration |
| scANVI | −0.180 ± 0.005 | 0.214 ± 0.004 | Label-informed integration |

Slide 10 displays means with 95% Student-t intervals computed from fold SEM; slide 11 displays fold SD. These are different uncertainty summaries. Overlapping source populations and recurring donors limit an independent-cohort interpretation of the intervals.

FE has higher cell-type ASW than Harmony, Scanorama, and scVI in this AML evaluation. Its raw batch ASW is closer to zero than theirs. Under the paper's raw batch-ASW readout this is consistent with greater residual batch structure, which is also visible in donor-colored UMAPs. Raw batch ASW is sensitive to cell composition and is not a complete measure of integration quality: interpretation also requires within-cell-type mixing and biological checks.

RE's strong batch separation agrees with its batch-specific objective. Its low cell-type ASW and fragmented cell-type UMAPs are expected consequences of emphasizing donor structure. RE is shown as a complementary representation, with its own objective.

scANVI has the highest cell-type ASW among these evaluated representations. It uses cell-type labels during training, so label availability must be stated when comparing it with unsupervised methods. The result supports this particular AML evaluation, with generalization to other datasets still to be tested.

Combined FE+RE downstream prediction was not evaluated in these results. The representation plots do not establish that the concatenated representation improves disease prediction.

## 3. UMAP figures and provenance

The paired FE/RE panels preserve the existing fold-1 training UMAP coordinates. The plotting script checks that all 23,049 rows have identical cell identifiers, cell types, batches, and patient-group labels in the two coordinate files. No UMAP refit or model training was performed during presentation preparation.

- `results/presentation/fe_re_celltype.png`: shared legend for 21 cell-type categories, with FE/RE panels.
- `results/presentation/fe_re_patient_group.png`: the same coordinates colored by AML, control, and cell-line groups.

Each representation has its own fitted UMAP geometry. Orientation, panel coordinates, and inter-island distances cannot be compared quantitatively across FE/RE. The interpretation focuses on label patterns and is paired with the formal five-fold test metrics. Patient group is linked to donor identity here; RE group separation alone does not prove disease-specific biology.

## 4. Donor-conditioned reconstruction

The analysis reads the 285 existing RE reconstruction arrays: five folds × three partitions × 19 target labels. The interpreted subset uses held-out test Mono/Mono-like source cells (996–999 per fold). Within each fold, the same source cells are reconstructed under all 19 target donor/batch labels.

For each gene, the effect is the equally weighted mean across 12 AML target donors minus the equally weighted mean across five control target donors. The two cell lines are displayed descriptively. Values are on the model-scaled expression scale.

Across 2,916 genes, 216 have the same effect direction in all five folds. Selected examples:

| Gene | Mean AML-target minus control-target difference | Same direction |
|---|---:|---:|
| SAMSN1 | −0.07132 | 5/5 folds |
| TXNIP | −0.08562 | 5/5 folds |
| FTL | −0.07749 | 5/5 folds |
| HIST1H1D | +0.06470 | 5/5 folds |
| CENPE | See selected-gene table | 4/5 folds |
| SRGN | See selected-gene table | 4/5 folds |

These genes are selected descriptive examples; selection was not preregistered. Fold 1 has zero genes with BH-FDR below 0.05 among the 2,916 tests. Direction stability across folds is a description of model outputs and does not itself establish statistical significance or causal disease effects.

`samsn1_target_donors.png` displays donor-specific predictions with fold SD. Removing BM5 changes the SAMSN1 mean difference from −0.07132 to −0.05649. SAMSN1, TXNIP, and FTL retain their negative mean direction when each control donor is omitted in turn, while magnitudes vary. The same target donors recur across all folds, so fold variation is not uncertainty across independent donor cohorts.

The appropriate interpretation is a model-generated **as-if simulation**: a hypothetical reconstructed profile conditioned on a substituted donor/batch label. Independent biological validation and a design separating disease, donor, and technical factors would be needed for stronger disease claims.

## 5. Re-render the presentation figures

Use an isolated plotting environment with `environment/presentation-plotting.txt`. This environment is separate from the recorded training environments.

```bash
python reproducibility/plot_presentation_results.py \
  --analysis-root /path/to/scMEDAL_FE_RE_analysis \
  --output-dir results/presentation
```

Required existing inputs under the analysis root:

- `AML_six_methods_train_split1/umap_19batches_seed_5/umap_csv/X_umap_scMEDAL-FE_latent_train_1_19batchsample.csv`
- Corresponding `scMEDAL-RE` coordinate CSV.
- `AML_counterfactual_fivefold_mono_monolike/aml_counterfactual_fivefold_gene_stability.csv`
- `aml_counterfactual_fivefold_selected_donor_means.csv` and `aml_counterfactual_fivefold_control_leaveoneout.csv` in that same directory.

The output manifest records relative source names and SHA-256 hashes. The repository publishes derived PNG/PDF figures and small selected summary tables. Original per-cell coordinate CSVs, expression data, reconstruction arrays, and weights remain outside this presentation archive.

## 6. Discussion and future work

The central finding is that FE and RE expose distinct representation patterns. FE retains some cross-donor cell-type organization; RE emphasizes donor-associated variation. Donor-conditioned reconstruction provides a way to inspect that modeled variation and its sensitivity to particular target donors.

Current boundaries are one AML dataset, author-provided cell-level folds, recurring target donors, label-informed scANVI, and a limited metric set. Prospective work can test donor-held-out generalization, marker/state preservation and within-cell-type mixing, combined FE+RE downstream prediction, and predefined validation-based FE loss-weight tuning. No improvement experiment is presented as completed.

## References

- [Andrade et al., scMEDAL, Nature Communications (2026)](https://doi.org/10.1038/s41467-026-72666-4).
- [Korsunsky et al., Harmony (2019)](https://doi.org/10.1038/s41592-019-0619-0).
- [Hie et al., Scanorama (2019)](https://doi.org/10.1038/s41587-019-0113-3).
- [Lopez et al., scVI (2018)](https://doi.org/10.1038/s41592-018-0229-2).
- [Xu et al., scANVI (2021)](https://doi.org/10.15252/msb.20209620).

All numerical findings above are derived from this repository's AML result tables and the existing local analysis outputs; the cited papers describe the underlying methods.
