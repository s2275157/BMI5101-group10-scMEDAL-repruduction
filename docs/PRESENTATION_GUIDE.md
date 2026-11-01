# Presentation guide

The central narrative is: source labels mix technical and biological variation;
FE and RE emphasize different information; the reported scores describe the
implemented pipelines; source-condition reconstructions provide exploratory
model associations.

| Slides | Purpose | Evidence or qualification |
|---|---|---|
| 1-4 | Dataset problem and motivation | Sources include patients, controls and cell lines. The illustrations are conceptual |
| 5-6 | Model mechanisms | Both networks reconstruct the input. FE/RE are trained separately. RE uses probabilistic source-dependent effects |
| 7 | Scope and protocol | Author-preprocessed AML, five cell-level folds, 50-dimensional outputs. Protocols differ across methods |
| 8 | Quantitative comparison | `results/metrics/test_metrics_long.csv`. Nominal fold intervals. Transductive Harmony/Scanorama and scaled scVI/scANVI inputs limit comparison |
| 9-11 | Qualitative FE/RE UMAPs | Fold-1 training cells. Coordinates are independently fitted by representation and reused when recoloring a representation |
| 12 | Source-condition reconstruction | Hold expression `x` fixed, supply `z_target` to both encoder and decoder, recompute the latent |
| 13 | Directional stability | `aml_counterfactual_fivefold_gene_stability.csv`. 216 of 2,916 genes. Same target sources recur |
| 14 | SAMSN1 donor sensitivity | `aml_counterfactual_fivefold_control_leaveoneout.csv`. -0.0713 becomes -0.0565 after omitting BM5 |
| 15 | Bounded weight screening | Documented validation summaries. Raw tuning result CSVs are not publicly archived |
| 16-17 | Interpretation and conclusion | Descriptive pipeline comparison, exploratory source effects, untested new-donor generalization |
| 18 | Contributions | Kept blank at the user's request |
| 19-20 | Reflection and references | Distinguish the supplied 2025 preprint from the later 2026 journal publication |

## Short answers for likely questions

**Did you reproduce the whole paper?** We reproduced selected AML computations
using the authors' preprocessed data and added a scoped benchmark and stability
analysis. We did not redo preprocessing, every dataset or the FE+RE predictor.
Exact equality with the paper's numbers has not been established.

**Is the comparison fully held out?** The same test subsets are scored, but
Harmony and Scanorama jointly integrate all partitions. They are transductive.
scVI/scANVI also inherit scaled input, which limits interpretation of their scores.

**What stays fixed in the reconstruction?** The source expression matrix and
selected source cells stay fixed. The target source label enters both encoder
and decoder, so the latent changes for each target.

**Is lower batch ASW always better?** Raw batch ASW measures grouping of labels.
Interpret it with cell-type preservation, composition and source confounding.
A more negative value alone is insufficient to rank integration quality.

**Are the 216 genes validated biomarkers?** They have a consistent modeled
direction across five overlapping folds. This is exploratory agreement, with
reused target donors and no independent biological validation.

**Did tuning improve FE?** The documented validation means are near-equal for
the original setting and adversarial weight 2. We retained the original setting.
Raw tuning result tables are needed for independent verification of the run.

## Reconciliation with the supplied results

The full 2,916-gene fold-1 table and 19 target-source means now support direct
verification of the reported zero significant genes (minimum BH-FDR 0.1142).
Slide 13 notes can cite this complete calculation. The archived six-method
coordinates are canonical for slides 9-11; their 23,049 cells align exactly.
The sensitivity slide is a separately documented experiment whose source score
CSVs are absent from the supplied archive. The contribution table remains blank.
