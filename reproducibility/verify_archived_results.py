#!/usr/bin/env python3
"""Recalculate the public AML summaries without private model arrays or a GPU."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

METHODS = ('scMEDAL-FE', 'scMEDAL-RE', 'Harmony', 'Scanorama', 'scVI', 'scANVI')
LABELS = ('batch', 'celltype')
METRICS = ('ch', '1/db', 'silhouette')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, message):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(a.shape == b.shape, message + ': shape mismatch')
    require(np.isfinite(a).all() and np.isfinite(b).all(), message + ': nonfinite value')
    require(np.allclose(a, b, rtol=1e-8, atol=1e-10), message)


def bh_adjust(p_values):
    values = np.asarray(p_values, dtype=float)
    require(np.isfinite(values).all() and ((0 <= values) & (values <= 1)).all(),
            'P-values must be finite and between zero and one')
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.empty_like(values)
    adjusted[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0, 1)
    return adjusted


def verify(repo_root):
    root = Path(repo_root)
    metrics = pd.read_csv(root / 'results/metrics/test_metrics_long.csv')
    expected = {(m, l, s) for m in METHODS for l in LABELS for s in METRICS}
    keys = list(zip(metrics['method'], metrics['label'], metrics['metric']))
    require(len(keys) == len(expected) and set(keys) == expected,
            'Expected one row for each of the 36 method/label/metric combinations')
    close(metrics['sem'], metrics['std'] / np.sqrt(5), 'Metric SEM differs from SD / sqrt(5)')
    require((metrics['std'] >= 0).all(), 'Metric SD must be nonnegative')

    raw = pd.read_csv(root / 'results/metrics/fe_re_test_per_fold.csv')
    raw_expected = {(m, f, l, s) for m in METHODS[:2] for f in range(1, 6)
                    for l in LABELS for s in METRICS}
    raw_keys = list(zip(raw['method'], raw['fold'], raw['label'], raw['metric']))
    require(len(raw_keys) == len(raw_expected) and set(raw_keys) == raw_expected,
            'Expected 60 unique FE/RE fold/label/metric measurements')
    require((raw['sample_size_setting'] == 10000).all(), 'Unexpected metric sample-size setting')
    for (method, label, metric), rows in raw.groupby(['method', 'label', 'metric']):
        summary = metrics[(metrics['method'] == method) & (metrics['label'] == label)
                          & (metrics['metric'] == metric)].iloc[0]
        values = rows['value'].to_numpy()
        close([values.mean(), values.std(ddof=1)], [summary['mean'], summary['std']],
              'Per-fold mean/SD mismatch: {}/{}/{}'.format(method, label, metric))

    directory = root / 'results/counterfactual/fold1_test_mono_monolike'
    targets = pd.read_csv(directory / 'target_batch_gene_means.csv')
    if 'target_batch' not in targets and 'index' in targets:
        targets = targets.rename(columns={'index': 'target_batch'})
    require(not targets['target_batch'].duplicated().any(), 'Duplicate target source')
    require(targets['target_group'].value_counts().to_dict() == {'AML': 12, 'control': 5, 'cellline': 2},
            'Expected 12 AML, five control and two cell-line targets')
    gene_columns = ['gene_{}'.format(i) for i in range(2916)]
    require(set(targets.columns) == {'target_batch', 'target_group', *gene_columns},
            'Unexpected target-gene table schema')
    aml = targets.loc[targets['target_group'] == 'AML', gene_columns].to_numpy()
    control = targets.loc[targets['target_group'] == 'control', gene_columns].to_numpy()
    require(np.isfinite(targets[gene_columns].to_numpy()).all(), 'Nonfinite target-gene mean')
    u, p_values = mannwhitneyu(aml, control, axis=0, alternative='two-sided')
    adjusted = bh_adjust(p_values)
    stats = pd.read_csv(directory / 'aml_vs_control_gene_stats.csv').sort_values('gene_index')
    require(np.array_equal(stats['gene_index'].to_numpy(), np.arange(2916)), 'Gene index mismatch')
    for column, calculated in [('aml_target_mean', aml.mean(0)), ('control_target_mean', control.mean(0)),
                               ('aml_minus_control', aml.mean(0) - control.mean(0)),
                               ('abs_difference', np.abs(aml.mean(0) - control.mean(0))),
                               ('mannwhitney_u', u), ('p_value', p_values), ('fdr_bh', adjusted)]:
        close(stats[column], calculated, 'Recalculated fold-1 ' + column + ' differs from archive')
    genes = pd.read_csv(directory / 'gene_index.csv').sort_values('gene_index')
    require(genes[['gene_index', 'gene']].reset_index(drop=True).equals(
        stats[['gene_index', 'gene']].reset_index(drop=True)), 'Gene-name mapping mismatch')
    top = pd.read_csv(directory / 'top_genes.csv')
    wanted_top = stats.sort_values(['abs_difference', 'p_value'], ascending=[False, True]).head(30)
    require(top['gene_index'].tolist() == wanted_top['gene_index'].tolist(), 'Top-30 ranking mismatch')
    close(top['fdr_bh'], wanted_top['fdr_bh'], 'Top-30 FDR mismatch')
    source_counts = pd.read_csv(directory / 'source_cell_summary.csv')
    require(set(source_counts['celltype']) == {'Mono', 'Mono-like'}, 'Unexpected source cell type')
    require(int(source_counts['n_source_cells'].sum()) == 999, 'Fold-1 source count differs from 999')

    five = root / 'results/counterfactual/fivefold_mono_monolike'
    stable = pd.read_csv(five / 'aml_counterfactual_fivefold_gene_stability.csv').sort_values('gene_index')
    require(stable[['gene_index', 'gene']].reset_index(drop=True).equals(
        stats[['gene_index', 'gene']].reset_index(drop=True)), 'Five-fold gene mapping mismatch')
    delta_columns = ['delta_fold{}'.format(f) for f in range(1, 6)]
    deltas = stable[delta_columns].to_numpy()
    close(stable['delta_fold1'], stats['aml_minus_control'], 'Fold-1 effects differ across tables')
    close(stable['mean_delta'], deltas.mean(1), 'Five-fold mean effect mismatch')
    close(stable['sd_delta_across_folds'], deltas.std(1, ddof=1), 'Five-fold effect SD mismatch')
    same_direction = np.maximum((deltas > 0).sum(1), (deltas < 0).sum(1))
    require(np.array_equal(stable['same_direction_folds'].to_numpy(), same_direction),
            'Per-gene direction counts differ')
    consistent = int((same_direction == 5).sum())
    require(consistent == 216, 'Common-direction gene count differs from 216')
    ranks = []
    for fold in range(1, 6):
        calculated = stable['delta_fold{}'.format(fold)].abs().rank(method='min', ascending=False)
        close(stable['abs_rank_fold{}'.format(fold)], calculated, 'Per-fold absolute-effect rank differs')
        ranks.append(calculated.to_numpy())
    close(stable['top30_fold_count'], (np.array(ranks) <= 30).sum(0), 'Top-30 fold count differs')

    donors = pd.read_csv(five / 'aml_counterfactual_fivefold_selected_donor_means.csv')
    require(len(donors) == 95 and not donors.duplicated(['gene', 'target_batch']).any(),
            'Expected 95 unique selected-gene/target-source rows')
    donor_values = donors[['value_fold{}'.format(f) for f in range(1, 6)]].to_numpy()
    close(donors['mean_across_folds'], donor_values.mean(1), 'Donor mean mismatch')
    close(donors['sd_across_folds'], donor_values.std(1, ddof=1), 'Donor SD mismatch')
    target_lookup = targets.set_index('target_batch')
    for row in donors.itertuples():
        require(row.target_group == target_lookup.loc[row.target_batch, 'target_group'],
                'Selected donor group differs from canonical target mapping')
        close([row.value_fold1], [target_lookup.loc[row.target_batch, 'gene_{}'.format(row.gene_index)]],
              'Selected fold-1 target mean mismatch')

    loo = pd.read_csv(five / 'aml_counterfactual_fivefold_control_leaveoneout.csv')
    require(len(loo) == 25 and not loo.duplicated(['gene', 'omitted_control']).any(),
            'Expected 25 unique selected-gene/control-omission rows')
    stable_lookup = stable.set_index('gene')
    for row in loo.itertuples():
        gene_rows = donors[donors['gene'] == row.gene]
        a = gene_rows.loc[gene_rows['target_group'] == 'AML', ['value_fold{}'.format(f) for f in range(1, 6)]].mean().to_numpy()
        c = gene_rows.loc[(gene_rows['target_group'] == 'control')
                         & (gene_rows['target_batch'] != row.omitted_control),
                         ['value_fold{}'.format(f) for f in range(1, 6)]].mean().to_numpy()
        predicted = a - c
        stored = [getattr(row, 'leaveoneout_delta_fold{}'.format(f)) for f in range(1, 6)]
        close(stored, predicted, 'Leave-one-control-out fold values differ')
        close([row.leaveoneout_mean_delta, row.full_mean_delta],
              [predicted.mean(), stable_lookup.loc[row.gene, 'mean_delta']],
              'Leave-one-control-out summary differs')
        close([row.shift_from_full_mean], [predicted.mean() - row.full_mean_delta],
              'Leave-one-control-out shift differs')
        require(int(row.same_direction_folds) == max(int((predicted > 0).sum()), int((predicted < 0).sum())),
                'Leave-one-control-out direction count differs')
        require(bool(row.direction_preserved) == bool(row.full_mean_delta * predicted.mean() > 0),
                'Leave-one-control-out direction flag differs')
    samsn = loo[(loo['gene'] == 'SAMSN1') & (loo['omitted_control'] == 'BM5')].iloc[0]
    return {
        'status': 'PASS', 'scope': 'Public archived summaries recalculated; no model training or private array revalidation',
        'six_method_aggregate_rows': 36, 'fe_re_per_fold_measurements': 60,
        'fold1_genes_retested': 2916, 'fold1_target_group_counts': {'AML': 12, 'control': 5, 'cellline': 2},
        'fold1_source_cells': 999, 'fold1_bh_fdr_min': float(adjusted.min()),
        'fold1_bh_fdr_below_0_05': int((adjusted < 0.05).sum()),
        'fivefold_common_direction_genes': consistent,
        'samsn1_full_mean_delta': float(samsn['full_mean_delta']),
        'samsn1_without_bm5_mean_delta': float(samsn['leaveoneout_mean_delta']),
        'other_four_methods_per_fold_measurements_supplied': False,
        'tuning_score_csvs_supplied': False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = verify(args.repo_root)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
