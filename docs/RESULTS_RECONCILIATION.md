# Reconciliation with the supplied AML results

Reconciled on 1 November 2026 against `scMEDAL_FE_RE_analysis.zip`.
Source SHA-256: `e4994533be602735abd6c7bc52c87e40d4022fd4f89cea107e54f626aec619c5`.
The ZIP contains 147 meaningful files after excluding Finder metadata and macOS
resource forks. The public repository before reconciliation was commit
`3cad8227decdf0ab9365730a1822cc55032e5c84`.

## Results that agree

| Evidence | Check | Outcome |
|---|---|---|
| Six-method aggregates | The archive's 36-row metric table and six-method mean table versus the original public tables | Exact byte matches |
| FE/RE fold measurements | Recompute all six metrics' means and SDs from five folds for each of FE and RE | Matches the published summaries |
| Fold-1 full gene analysis | Recompute all 2,916 AML/control target-mean differences, U statistics, p-values and BH-FDR values | Matches to numerical precision |
| Fold-1 FDR | Count BH-FDR below 0.05 across all genes | 0; minimum BH-FDR 0.1142387 |
| Five-fold gene stability | Recompute each gene's mean, SD and sign count | 216 of 2,916 genes have a common sign in all five folds |
| Selected target responses | Recompute the 95 selected gene/source means and SDs; compare fold-1 values with the complete target table | Matches |
| Control sensitivity | Recompute all 25 selected-gene/control-omission rows | Matches |
| SAMSN1 excluding BM5 | Compare the actual mean with the original PPT label | Correct display is -0.0565, replacing -0.057 |
| Canonical coordinate tables | Check all seven methods' ordered cell identifiers and labels | Same 23,049 cells, 19 sources, 21 cell types |
| Original PPT source-color figures | Compare both slide-9 PNGs with archive members | Exact byte matches to the six-method no-axis PNGs |

Run `python reproducibility/verify_archived_results.py` to repeat the public
numeric checks. The stored report is
[`results/provenance/reconciliation_report.json`](../results/provenance/reconciliation_report.json).

## Resolved evidence gap

The earlier public repository contained only the fold-1 top-30 gene list. The
supplied ZIP includes the complete fold-1 target-mean and all-gene test tables;
both are now added. The zero-significant-gene statement can therefore be
independently recalculated from public files. This strengthens traceability of
the existing exploratory result and does not make it independent biological
validation.

The old target-mean CSV calls its source column `index`. The curated copy calls
it `target_batch`; its numerical records are preserved. Both source and
published hashes are recorded. The older FE/RE multi-header score table is
reshaped into an explicit long table with 60 measurements. The value 10,000 is
the requested metric sample-size setting, rather than a count of training cells.

## UMAP versions and presentation

`AML_FE_RE_only_2026-09-23` and `AML_six_methods_train_split1` contain the same
cells but different coordinate fits. Their plotted locations cannot be spliced
into a single coordinate system. The six-method family is canonical for the
presentation and public curated figures.

The `X_umap_*` plots show a fitted UMAP. The separate `X_modelrawlatent_*` plots
show a different projection and are kept out of the curated UMAP set. Axis and
no-axis copies of a given plot are display variants, rather than new results.
`plot_archived_umap.py` replots the stored coordinates without fitting UMAP or
rerunning a model, and uses identical colors across representations.

Presentation slide 13 and its notes now cite the complete recalculated FDR
table. Slides 10 and 11 use reproducible plots from the same canonical stored
coordinates as slide 9. The supplementary sensitivity slide explicitly states
that its original score files are outside the supplied archive. The contribution
table remains blank, and the presentation itself stays local under repository
policy.

## Evidence that remains outside this ZIP

- Individual fold measurements for Harmony, Scanorama, scVI and scANVI. Their
  six-method aggregate summaries match exactly, but the underlying fold values
  cannot be recalculated from this ZIP.
- FE tuning score CSVs and completion logs. The retained sensitivity numbers
  remain separately documented summaries.
- Historical training source commit/patches, original model weights, expression
  input matrices and the 285 reconstruction arrays. Training completion and
  full-array validation were not repeated in this reconciliation.

The protocol limitations in [AUDIT_PROTOCOL.md](AUDIT_PROTOCOL.md) continue to
apply: transductive Harmony/Scanorama, scaled scVI/scANVI input, supervision and
validation differences, overlapping folds, and source/disease confounding.
