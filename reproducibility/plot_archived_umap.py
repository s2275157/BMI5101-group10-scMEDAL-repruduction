#!/usr/bin/env python3
"""Plot stored UMAP coordinates, keeping source-cell annotations local."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

PREFIX = ('scMEDAL_FE_RE_analysis/AML_six_methods_train_split1/'
          'umap_19batches_seed_5/umap_csv/')
METHOD_KEYS = [('Input PCA', 'input_pca'), ('scMEDAL-FE', 'scMEDAL-FE'),
               ('scMEDAL-RE', 'scMEDAL-RE'), ('Harmony', 'Harmony'),
               ('Scanorama', 'Scanorama'), ('scVI', 'scVI'), ('scANVI', 'scANVI')]


def load_coordinates(archive_path):
    frames, provenance = {}, []
    reference = None
    with ZipFile(archive_path) as archive:
        for name, key in METHOD_KEYS:
            suffix = '{}_latent_train_1_19batchsample.csv'.format(key) if key != 'input_pca' else 'input_pca_train_1_19batchsample.csv'
            member = PREFIX + 'X_umap_' + suffix
            data = archive.read(member)
            frame = pd.read_csv(io.BytesIO(data), usecols=lambda c: c.startswith('UMAP') or c in
                                {'Cell', 'batch', 'celltype', 'Patient_group'})
            xy = [c for c in frame if c.startswith('UMAP')]
            if len(xy) != 2 or not np.isfinite(frame[xy].to_numpy()).all():
                raise ValueError('Expected two finite UMAP coordinates for ' + name)
            frame = frame.rename(columns={xy[0]: 'x', xy[1]: 'y'})
            labels = frame[['Cell', 'batch', 'celltype', 'Patient_group']]
            if labels.isna().any().any() or frame['Cell'].duplicated().any():
                raise ValueError('Missing labels or duplicate cells for ' + name)
            if reference is None:
                reference = labels
            elif not reference.equals(labels):
                raise ValueError('Cell order or annotation mismatch for ' + name)
            if len(frame) != 23049 or frame['batch'].nunique() != 19 or frame['celltype'].nunique() != 21:
                raise ValueError('Unexpected fold-1 training dimensions for ' + name)
            frames[name] = frame
            provenance.append({'method': name, 'archive_path': member,
                               'source_sha256': hashlib.sha256(data).hexdigest(),
                               'rows': len(frame), 'sources': 19, 'celltypes': 21})
    return frames, provenance


def palette(column, labels):
    if column == 'Patient_group':
        return ['AML', 'control', 'cellline'], {'AML': '#ef1760', 'control': '#dfb800', 'cellline': '#3aa657'}
    categories = sorted(labels.astype(str).unique())
    colors = list(plt.get_cmap('tab20').colors) + ['#111111']
    return categories, dict(zip(categories, colors))


def render(frames, methods, column, output_dir, stem):
    reference = frames[methods[0]]
    categories, colors = palette(column, reference[column])
    paired = len(methods) == 2
    fig, axes = plt.subplots(1 if paired else 2, 2 if paired else 4,
                             figsize=(12, 4.8) if paired else (14, 8.0), squeeze=False)
    for ax, name in zip(axes.flat, methods):
        frame = frames[name]
        for category in categories:
            rows = frame[column] == category
            ax.scatter(frame.loc[rows, 'x'], frame.loc[rows, 'y'], s=0.55 if paired else 0.40,
                       c=[colors[category]], alpha=0.8, linewidths=0, rasterized=True)
        ax.set_title(name, fontweight='bold', fontsize=14 if paired else 12, color='#123b69', pad=9)
        ax.set_axis_off()
        ax.margins(.025)
    for ax in list(axes.flat)[len(methods):]:
        ax.set_visible(False)
    ncol = 7 if column == 'celltype' else (5 if column == 'batch' else 3)
    handles = [Line2D([], [], marker='o', linestyle='', markersize=4,
                      color=colors[c], label=c) for c in categories]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .03),
               ncol=ncol, frameon=False, fontsize=8, columnspacing=1.5, handletextpad=.3)
    fig.text(.5, .005, 'Fold-1 training cells (n = 23,049). Independently fitted UMAP by representation.',
             ha='center', fontsize=7.5, color='#667788')
    fig.subplots_adjust(left=.015, right=.985, top=.91, bottom=.23 if paired else .16,
                        wspace=.08, hspace=.20)
    destination = output_dir / (stem + '.png')
    fig.savefig(destination, dpi=200, facecolor='white')
    plt.close(fig)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    frames, provenance = load_coordinates(args.archive)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    assets = []
    for column, suffix in [('batch', 'batch'), ('celltype', 'celltype'), ('Patient_group', 'patient_group')]:
        for methods, stem in [(['scMEDAL-FE', 'scMEDAL-RE'], 'fe_re'),
                              ([m for m, _ in METHOD_KEYS], 'all_methods')]:
            image = render(frames, methods, column, args.output_dir,
                           '{}_{}_fold1_train'.format(stem, suffix))
            assets.append({'filename': image.name, 'sha256': hashlib.sha256(image.read_bytes()).hexdigest(),
                           'color_label': column, 'methods': methods})
    ref = frames['scMEDAL-FE']
    (ref.groupby(['Patient_group', 'batch', 'celltype'], observed=True).size()
        .rename('n_cells').reset_index().to_csv(args.output_dir / 'fold1_train_cell_counts.csv', index=False))
    source_groups = ref[['batch', 'Patient_group']].drop_duplicates()
    if source_groups['batch'].duplicated().any():
        raise ValueError('A source maps to multiple patient groups')
    source_groups.sort_values('batch').to_csv(args.output_dir / 'source_groups.csv', index=False)
    record = {'coordinate_family': 'AML_six_methods_train_split1', 'dataset': 'train', 'fold': 1,
              'cell_alignment_verified_across_all_methods': True, 'source_coordinates': provenance,
              'figures': assets, 'per_cell_coordinates_and_annotations_published': False,
              'scope': 'Replot stored coordinates without fitting UMAP or rerunning models'}
    (args.output_dir / 'figure_manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    print('ARCHIVED_UMAP_PLOTS_OK; figures={}; cells={}'.format(len(assets), len(ref)))


if __name__ == '__main__':
    main()
