import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('rna_resources', ROOT/'workflow/build_input_resources.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
GIB = 1024**3


class RnaResourceTests(unittest.TestCase):
    def test_large_sorted_bam_has_time_margin_and_no_64_thread_picard_claim(self):
        allocation = m.rna_allocation('GATK4_MARKDUPLICATES', 167*GIB)
        self.assertEqual((allocation['cpus'], allocation['memory_gb'], allocation['time_hours']), (8,128,22))
        self.assertGreater(m.rna_allocation('GATK4_MARKDUPLICATES', 300*GIB)['time_hours'],22)

    def test_rna_steps_use_different_actual_workloads(self):
        self.assertEqual(m.rna_allocation('STAR_ALIGN', 63*GIB)['time_hours'],6)
        self.assertEqual(m.rna_allocation('SAMTOOLS_SORT', 1633*GIB)['time_hours'],6)
        self.assertEqual(m.rna_allocation('SAMTOOLS_SORT', 1633*GIB)['memory_gb'],80)
        self.assertEqual(m.rna_allocation('ISOFOX', 299*GIB)['memory_gb'],176)

    def test_invalid_or_excessive_size_is_not_silently_clamped(self):
        for size in (0, -1, 100_000*GIB):
            with self.assertRaises(ValueError):m.rna_allocation('GATK4_MARKDUPLICATES',size)

    def test_dynamic_rules_do_not_change_dna_or_claim_intermediate_sizes_at_launch(self):
        rows=[dict(role=r,sample_id=r,read_pairs=None,input_bytes=80*GIB) for r in ['tumor','normal','rna']]
        full=m.plan(rows);dna=m.plan(rows[:2])
        for name in dna:
            # SAGE_APPEND previously depended on max(tumour,RNA); identical here.
            self.assertEqual(full[name],dna[name])
        for name in m.RNA_RESOURCES:
            self.assertTrue(full[name]['dynamic'])
            self.assertNotIn('memory_gb',full[name])
        config=m.render_config(full,m.redux_plan(rows))
        self.assertIn('reads_fwds, reads_revs',config)
        self.assertIn('bam instanceof Collection',config)
        self.assertIn('aln instanceof Collection',config)
        self.assertNotIn("withName: 'GATK4_MARKDUPLICATES'",config)

    def test_wgs_config_has_no_rna_variables(self):
        rows=[dict(role=r,sample_id=r,read_pairs=None,input_bytes=80*GIB) for r in ['tumor','normal']]
        config=m.render_config(m.plan(rows),m.redux_plan(rows))
        self.assertNotIn('reads_fwds',config)
        self.assertNotIn('RNA input-v5',config)

    def test_completed_measurements_fit_the_recorded_envelope(self):
        # Rounded, de-identified observed envelopes, not patient fixtures.
        for stage,gib,hours,rss in [
            ('STAR_ALIGN',79,3.26,96.3),('SAMTOOLS_SORT',2374,3.20,53.3),
            ('GATK4_MARKDUPLICATES',138,8.15,93.8),
            ('GATK4_MARKDUPLICATES',240,8.17,93.8),('ISOFOX',299,.52,103.1)]:
            allocation=m.rna_allocation(stage,gib*GIB)
            self.assertGreater(allocation['time_hours'],hours*1.5)
            self.assertGreater(allocation['memory_gb'],rss*1.25)
