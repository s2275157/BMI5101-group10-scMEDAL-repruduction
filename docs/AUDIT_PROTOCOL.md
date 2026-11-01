# AML protocol review and interpretation boundaries

Reviewed on 1 November 2026. Public reproduction snapshot:
`3cad8227decdf0ab9365730a1822cc55032e5c84`. Upstream source reviewed:
[`bbf9a0b`](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/tree/bbf9a0b801123d63218cf05d93ec8f103453881d).

## What this project establishes

The archived summary tables show that FE and RE encode different source-related
structures on the authors' preprocessed AML data. FE has higher mean cell-type
ASW than Harmony, Scanorama and scVI in these implemented runs. RE has high batch
ASW. These are descriptive results for the implemented pipelines.

The study covers one preprocessed AML dataset, six representations, five supplied
cell-level folds and an exploratory source-condition reconstruction analysis.
Preprocessing from FASTQ, donor-held-out generalization and a combined FE+RE
downstream predictor are outside the completed scope. Exact numerical agreement
with every experiment in the paper has not been demonstrated.

## Source versions

The supplied `scMEDAL.pdf` is the 19 March 2025 preprint,
[doi:10.21203/rs.3.rs-6081478/v1](https://doi.org/10.21203/rs.3.rs-6081478/v1).
It describes a two-neuron bottleneck in Methods 4.5 and AML batch-ASW values of
about -0.30 (FE) and +0.32 (RE) in section 2.4. This reproduction's run names
record 50 latent dimensions and its archived means are -0.0377 and +0.6850.
Those differences must be acknowledged when discussing reproduction fidelity.

The later journal article is
[Andrade et al., Nature Communications 17, 7075 (2026)](https://www.nature.com/articles/s41467-026-72666-4).
Publication metadata were verified. Complete numerical reconciliation with that
later article was not completed in this audit because direct full-text access
was inconsistent. Do not attribute the preprint's particular settings or figure
numbers to the later article without checking that version.

The current upstream source matches the corresponding local reference files
after trailing-newline normalization. The public training manifest does not
record the historical upstream commit, patches, input checksums or full package
locks. Inspection of current source is therefore evidence of that pipeline's
behavior, rather than proof that an unrecorded historical patch never changed it.

## The comparison protocols differ

| Representation | Construction and data access in reviewed source | Consequence |
|---|---|---|
| FE / RE | Fit using the training partition and validation procedure in the TensorFlow pipeline | Test-cell encoding is separate from model fitting, subject to upstream preprocessing choices |
| Harmony | Concatenate train, validation and test cells, run PCA and Harmony, then split back for scoring | Transductive integration includes test expression |
| Scanorama | Concatenate all three partitions, run PCA and Scanorama, then split back | Transductive integration includes test expression |
| scVI / scANVI | Instantiate on the supplied training partition and reserve `validation_size=0.1` internally | The supplied validation subset does not control early stopping in the same way as FE/RE |
| scANVI | Uses cell-type labels during training | Its cell-type ASW reflects a different supervision condition |

Evidence:

- [Harmony concatenation and integration](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/blob/bbf9a0b801123d63218cf05d93ec8f103453881d/comparables/harmony.py)
- [Scanorama concatenation and integration](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/blob/bbf9a0b801123d63218cf05d93ec8f103453881d/comparables/scanorama.py)
- [scVI training and input registration](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/blob/bbf9a0b801123d63218cf05d93ec8f103453881d/comparables/scVI.py)
- [scANVI training and input registration](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/blob/bbf9a0b801123d63218cf05d93ec8f103453881d/comparables/scANVI.py)

## scVI/scANVI input limitation

The run names specify `log_transformed_2916hvggenes`, `scaling-min_max` and a
ZINB likelihood. `load_data()` scales `X` and constructs an AnnData object from
that matrix. The scVI/scANVI wrappers call `setup_anndata()` without selecting
a raw-count layer. Standard count-based scvi-tools models ordinarily require
raw counts. See the
[official scvi-tools data guidance](https://docs.scvi-tools.org/en/latest/tutorials/notebooks/quick_start/api_overview.html).

The current scores remain a record of the inherited pipeline. They do not
establish how standard raw-count scVI/scANVI would perform. A corrected-input
benchmark would require new training and a separately identified result table.
Retrospectively relabeling the archived scores as corrected would be inaccurate.

## What the RE reconstruction actually holds fixed

The generator calls `trained_model.predict((x, z_target))` for each target.
The RE model's `call()` uses that target label in both stages:

```text
r_target = re_encoder((x, z_target))
xhat_target = re_decoder((r_target, z_target))
```

The source expression and selected source-cell rows remain fixed. The latent
representation changes with the target label. A fixed-latent decoder-only swap
is a distinct proposed analysis and is not the generator of the archived arrays.

Evidence: [generation loop](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/blob/bbf9a0b801123d63218cf05d93ec8f103453881d/utils/model_train_utils.py),
[RE model call](https://github.com/DeepLearningForPrecisionHealthLab/scMEDAL_for_scRNAseq/blob/bbf9a0b801123d63218cf05d93ec8f103453881d/models/scMEDAL/scMEDAL.py).

The 19 target labels include 12 AML donors, 5 controls and 2 cell lines. For
each target, the analysis first averages across the selected Mono/Mono-like
source cells; it then compares the 12 AML target means with the 5 control target
means. The cell lines are descriptive targets and are excluded from this group
difference. These are learned source-condition associations. Disease, individual
biology and technical effects remain confounded.

## Metrics, intervals and inferential scope

- Raw ASW measures label clustering. Higher cell-type ASW indicates clearer
  annotated grouping. Weak batch-label separation can be useful, but a more
  negative batch ASW is not automatically better correction.
- `test_metrics_long.csv` contains 36 aggregated rows: 6 methods, 2 labels and
  3 metrics. The supplied archive additionally provides the FE/RE fold measurements,
  now in `results/metrics/fe_re_test_per_fold.csv`. The four other methods have
  aggregate summaries only in this archive.
- The plotted interval is `mean ± t(0.975, df=4) × SEM`. Folds reuse sources and
  overlap. Treat this as a nominal descriptive interval, without a claim of
  independent donor-level coverage or statistical significance between methods.
- The selected six gene rows on PPT slide 13 agree with the archived five-fold
  gene-stability table. Exactly 216 of 2,916 genes have a common sign in all
  five folds. Sign agreement alone does not establish biological replication.
- The supplied result ZIP contains complete fold-1 target means and all 2,916
  gene tests. These are now public. All U statistics, p-values and BH adjustments
  were recalculated, with zero genes below FDR 0.05 and minimum FDR 0.1142.
  This resolves the earlier evidence gap limited to the public top-30 table.
- SAMSN1's full mean is -0.0713217288 and its BM5-excluded mean is
  -0.0564944827 (rounded). The corrected display uses -0.0713 and -0.0565.
- Tuning summaries agree between the README, tuning document and PPT. Original
  tuning score CSVs/logs are absent, and training seeds were not matched in the
  historical comparison. The reported near-tie is a descriptive conclusion.

## What a stronger follow-up would need

Record the historical source commit and applied patches, dataset checksums,
the four other methods' per-fold metrics and tuning validation tables. Complete
fold-1 all-gene statistics and FE/RE fold metrics have now been added.
For a new benchmark, use raw counts for count-based models, state transductive
versus inductive access explicitly, reserve donors for testing where the model
supports it, and evaluate a defined downstream task for FE+RE. A method that
requires learned one-hot donor effects needs an explicit strategy for previously
unseen donor labels before claiming donor-held-out generalization.

## Supplied archive reconciliation

The archive `scMEDAL_FE_RE_analysis.zip` has SHA-256
`e4994533be602735abd6c7bc52c87e40d4022fd4f89cea107e54f626aec619c5`.
The complete target-mean and gene-statistic tables resolve the prior fold-1 FDR
evidence gap. Seven coordinate CSVs contain the same ordered 23,049 training
cells, 19 sources and 21 cell types. The source-only and six-method analyses
contain different coordinate fits; the six-method family is canonical for
presentation. See [the reconciliation record](RESULTS_RECONCILIATION.md).
