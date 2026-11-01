"""Regression checks for corrupt evidence tables and the archived CSV schema."""
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'reproducibility'))
from verify_archived_results import verify


class ArchiveVerificationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        shutil.copytree(ROOT / 'results', self.root / 'results',
                        ignore=shutil.ignore_patterns('*.png', '*.pdf'))

    def tearDown(self):
        self.temp.cleanup()

    def test_archived_evidence_recalculates(self):
        result = verify(self.root)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['fold1_bh_fdr_below_0_05'], 0)
        self.assertEqual(result['fivefold_common_direction_genes'], 216)

    def test_tampered_fdr_is_rejected(self):
        path = self.root / 'results/counterfactual/fold1_test_mono_monolike/aml_vs_control_gene_stats.csv'
        table = pd.read_csv(path)
        table.loc[0, 'fdr_bh'] = 0.001
        table.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, 'fdr_bh'):
            verify(self.root)

    def test_duplicate_fold_measurement_is_rejected(self):
        path = self.root / 'results/metrics/fe_re_test_per_fold.csv'
        table = pd.read_csv(path)
        pd.concat([table, table.iloc[[0]]], ignore_index=True).to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, '60 unique'):
            verify(self.root)

    def test_legacy_target_header_remains_readable(self):
        path = self.root / 'results/counterfactual/fold1_test_mono_monolike/target_batch_gene_means.csv'
        table = pd.read_csv(path).rename(columns={'target_batch': 'index'})
        table.to_csv(path, index=False)
        self.assertEqual(verify(self.root)['status'], 'PASS')


if __name__ == '__main__':
    unittest.main()
