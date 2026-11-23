import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('reference_dispatch',
    Path(__file__).resolve().parents[1]/'workflow/submit_reference_deepsomatic.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ReferenceComparisonTests(unittest.TestCase):
    def test_allocation_grows_from_measured_workload(self):
        small = m.resources(m.BASELINE_BYTES)['stages']
        large = m.resources(2*m.BASELINE_BYTES)['stages']
        self.assertEqual(small['examples']['cpus'], 64)
        self.assertEqual(small['inference']['cpus'], 16)
        for stage in small:
            self.assertGreater(large[stage]['hours'], small[stage]['hours'])
            self.assertGreaterEqual(large[stage]['mem_gb'], small[stage]['mem_gb'])
        with self.assertRaises(ValueError):
            m.resources(m.BASELINE_BYTES*20)

    def test_cpu_gpu_separation_and_full_caller(self):
        cmds = m.commands(Path('/frozen'), Path('/case.json'), Path('/run/results'),
                          m.resources(m.BASELINE_BYTES))
        self.assertEqual(list(cmds), ['examples', 'inference', 'postprocess'])
        for stage, cmd in cmds.items():
            self.assertEqual('--gres=gpu:1' in cmd, stage == 'inference')
            self.assertIn('/frozen/deepsomatic_full.py', cmd)
            self.assertNotIn('--regions', ' '.join(cmd))
