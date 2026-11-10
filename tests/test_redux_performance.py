import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('benchmark',ROOT/'workflow/benchmark_redux_performance.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ReduxExperimentTests(unittest.TestCase):
    def test_trials_keep_scientific_arguments_and64_threads(self):
        source = 'redux -Xmx999 -threads 12 -input_bam source.bam -form_consensus -unmap_regions regions.tsv'
        for mode,flags in module.MODES.items():
            result = module.experiment_script(source,mode)
            self.assertEqual(result,source.replace('-Xmx999','-Xmx128849018880').replace('-threads 12','-threads 64 '+flags))
        self.assertNotIn('-partition_ratio',module.experiment_script(source,'baseline'))

    def test_non_native_or_preconfigured_command_is_refused(self):
        for source in ['redux -Xmx999','redux -threads 12',
                       'redux -Xmx999 -threads 12 -parallel_concat']:
            with self.assertRaises(ValueError):
                module.experiment_script(source,'baseline')
