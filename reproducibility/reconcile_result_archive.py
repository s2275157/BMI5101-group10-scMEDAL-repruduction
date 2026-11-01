#!/usr/bin/env python3
"""Import small AML summaries from the supplied result ZIP and recalculate them.

Raw per-cell annotations, source ZIP files and model arrays are never copied to
the public results tree. macOS resource forks and Finder metadata are ignored.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import stat
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

import pandas as pd

from verify_archived_results import verify

ROOT = 'scMEDAL_FE_RE_analysis/'
CF1 = 'AML_counterfactual_fold1_test_mono_monolike/'
CF5 = 'AML_counterfactual_fivefold_mono_monolike/'
BENCHMARK = 'AML_six_methods_fivefold_metrics/metrics/'


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def import_archive(archive_path, output_root):
    output_root = Path(output_root)
    copied = []
    with ZipFile(archive_path) as archive:
        inventory = []
        for member in archive.infolist():
            path = PurePosixPath(member.filename)
            if member.is_dir() or path.parts[0] == '__MACOSX' or path.name in {'.DS_Store'} or path.name.startswith('._'):
                continue
            if path.is_absolute() or '..' in path.parts or stat.S_ISLNK(member.external_attr >> 16):
                raise ValueError('Unsafe ZIP member: ' + member.filename)
            if not member.filename.startswith(ROOT):
                raise ValueError('Unexpected archive root: ' + member.filename)
            inventory.append({'archive_path': member.filename[len(ROOT):], 'bytes': member.file_size,
                              'sha256': hashlib.sha256(archive.read(member)).hexdigest()})

        def write(source, destination, rename_index=False):
            data = archive.read(ROOT + source)
            source_hash = hashlib.sha256(data).hexdigest()
            renamed = rename_index and data.startswith(b'index,target_group,')
            if renamed:
                data = b'target_batch' + data[len(b'index'):]
            target = output_root / destination
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            copied.append({'archive_path': source, 'repository_path': destination,
                           'source_sha256': source_hash, 'published_sha256': hashlib.sha256(data).hexdigest(),
                           'transformation': 'rename legacy index header to target_batch' if renamed else 'byte-preserving copy'})

        for name in ['test_metrics_long.csv', 'test_metrics_mean.csv']:
            write(BENCHMARK + name, 'results/metrics/' + name)
        for name in ['analysis_metadata.json', 'gene_index.csv', 'source_cell_summary.csv',
                     'aml_vs_control_gene_stats.csv', 'target_batch_gene_means.csv', 'top_genes.csv',
                     'aml_counterfactual_top_genes_heatmap.png', 'aml_counterfactual_top_genes_heatmap.pdf']:
            write(CF1 + name, 'results/counterfactual/fold1_test_mono_monolike/' + name,
                  rename_index=name == 'target_batch_gene_means.csv')
        for member in sorted(archive.namelist()):
            prefix = ROOT + CF5
            if not member.startswith(prefix) or member.endswith('/'):
                continue
            name = member[len(prefix):]
            if '/' in name or not name.endswith(('.csv', '.json', '.png', '.pdf')):
                continue
            write(CF5 + name, 'results/counterfactual/fivefold_mono_monolike/' + name)

        source = 'AML_FE_RE_only_2026-09-23/test_allscores.csv'
        raw = pd.read_csv(io.BytesIO(archive.read(ROOT + source)), header=[0, 1])
        def column(name):
            matches = [c for c in raw.columns if c[0] == name]
            if len(matches) != 1:
                raise ValueError('Ambiguous archived column: ' + name)
            return matches[0]
        models = {'scmedalfe': 'scMEDAL-FE', 'scmedalre': 'scMEDAL-RE'}
        rows = []
        for _, row in raw.iterrows():
            method = models[str(row[column('latent_name')])]
            for label in ['batch', 'celltype']:
                for metric in ['ch', '1/db', 'silhouette']:
                    rows.append({'method': method, 'fold': int(row[column('fold')]), 'label': label,
                                 'metric': metric, 'value': float(row[(label, metric)]),
                                 'sample_size_setting': int(row[column('sample_size')])})
        fold_table = pd.DataFrame(rows).sort_values(['method', 'fold', 'label', 'metric'])
        fold_path = output_root / 'results/metrics/fe_re_test_per_fold.csv'
        fold_table.to_csv(fold_path, index=False)
        copied.append({'archive_path': source, 'repository_path': str(fold_path.relative_to(output_root)),
                       'source_sha256': hashlib.sha256(archive.read(ROOT + source)).hexdigest(),
                       'published_sha256': sha256_file(fold_path),
                       'transformation': 'reshape multi-header fold table into 60 explicit metric rows; sample size is the requested metric setting'})

        original_manifest = json.loads(archive.read(ROOT + 'AML_FE_RE_only_2026-09-23/analysis_manifest.json'))
        sanitized = {key: original_manifest[key] for key in ['created_at', 'python', 'analysis_name',
                                                            'model_folders', 'metrics_scope', 'umap_scope']}
        sanitized['private_output_paths_removed'] = True
        provenance = output_root / 'results/provenance'
        provenance.mkdir(parents=True, exist_ok=True)
        (provenance / 'fe_re_analysis_manifest.json').write_text(json.dumps(sanitized, indent=2) + '\n')
        pd.DataFrame(inventory).to_csv(provenance / 'archive_inventory.csv', index=False)
        source_record = {'source_archive': Path(archive_path).name, 'source_sha256': sha256_file(archive_path),
                         'meaningful_archive_files': len(inventory), 'copied_or_derived': copied,
                         'excluded_from_public_copy': ['per-cell UMAP CSVs with source annotations',
                                                      'macOS resource forks and Finder metadata',
                                                      'legacy/duplicate figure variants', 'source ZIP'],
                         'scope': 'Archived summary evidence; no fresh model training'}
        (provenance / 'source_archive.json').write_text(json.dumps(source_record, indent=2) + '\n')

    report = verify(output_root)
    report['source_archive_sha256'] = source_record['source_sha256']
    (provenance / 'reconciliation_report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True,
                        help='Repository root or a new staging directory for the selected summaries')
    args = parser.parse_args()
    print(json.dumps(import_archive(args.archive, args.output_root), indent=2))


if __name__ == '__main__':
    main()
