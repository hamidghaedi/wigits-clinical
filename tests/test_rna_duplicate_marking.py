import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('rna_duplicate_test',ROOT/'workflow/compare_rna_duplicate_marking.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class RnaDuplicateMarkingTests(unittest.TestCase):
    def test_synthetic_fixture_is_coordinate_sorted_and_paired(self):
        rows=[line.split('\t') for line in m.synthetic_sam().splitlines() if not line.startswith('@')]
        self.assertEqual(len(rows),30)
        positions=[int(row[3]) for row in rows]
        self.assertEqual(positions,sorted(positions))
        names={r[0] for r in rows}
        self.assertEqual(len(names),15)
        for name in names:
            self.assertEqual({int(r[1]) for r in rows if r[0]==name},{99,147})

    def test_policy_changes_only_optical_detection_in_rna_process(self):
        config=(ROOT/'workflow/conf/rna_duplicate_marking.config').read_text()
        self.assertIn("withName: '.*:READ_ALIGNMENT_RNA:GATK4_MARKDUPLICATES'",config)
        self.assertIn("ext.args = '--READ_NAME_REGEX null'",config)
        for setting in ['REMOVE_DUPLICATES','REMOVE_SEQUENCING_DUPLICATES','SAGE','REDUX','STAR_ALIGN']:
            self.assertNotIn(setting,config)

    def test_full_native_entry_loads_the_shared_policy_before_explicit_overrides(self):
        launcher=(ROOT/'workflow/run_full_native.sbatch').read_text()
        self.assertIn('-c "$PACKAGE/workflow/conf/rna_duplicate_marking.config" "$@"',launcher)
