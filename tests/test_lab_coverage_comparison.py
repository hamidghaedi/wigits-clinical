import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('c',Path(__file__).parents[1]/'workflow/compare_lab_coverage.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)


class CoverageComparisonTests(unittest.TestCase):
    def test_numeric_formatting_does_not_create_a_difference(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a';b=Path(tmp)/'b'
            a.write_text('GENOME_TERRITORY\tMEAN_COVERAGE\n100\t20.0\n\ncoverage\thigh_quality_coverage_count\n0\t10\n20\t90\n')
            b.write_text(a.read_text().replace('20.0','20.000000'))
            r=c.compare(a,b);self.assertTrue(r['metrics_equal']);self.assertTrue(r['histogram_equal'])

    def test_metric_and_histogram_changes_are_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a';b=Path(tmp)/'b'
            a.write_text('GENOME_TERRITORY\tMEAN_COVERAGE\n100\t20\n\ncoverage\thigh_quality_coverage_count\n0\t10\n20\t90\n')
            b.write_text(a.read_text().replace('100\t20','100\t19').replace('20\t90','20\t89'))
            r=c.compare(a,b);self.assertFalse(r['metrics_equal']);self.assertFalse(r['histogram_equal'])
            self.assertEqual(r['metric_differences'][0]['field'],'MEAN_COVERAGE')

    def test_missing_measurements_are_not_silently_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a';a.write_text('GENOME_TERRITORY\tMEAN_COVERAGE\n100\t20\n')
            with self.assertRaises(ValueError):c.sections(a)
