import copy
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow"))
import optional_extensions as planner
import extension_tasks as tasks
import neo_interface as ni


class OptionalExtensionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.case_root = self.root / "native"
        shutil.copytree(ROOT / "tests/fixtures/neo_interface/subject_a", self.case_root)
        layout = ni.CaseLayout.from_results_dir(str(self.case_root))
        self.case = dict(case_id="SYNTHETIC", tumour_sample=layout.tumour_dna, normal_sample=layout.normal_dna,
                         rna_sample=layout.tumour_rna, results_case_dir=str(self.case_root), native_genome="GRCh38_hmf",
                         options=dict(blood_ch=True, rna_neoantigens=True),
                         laboratory=dict(source="synthetic_laboratory.csv", normal_specimen_type="Blood", tumour_rna_available=True))
        for sample in (layout.tumour_dna, layout.normal_dna):
            bam = self.case_root / "alignments" / sample / (sample + ".redux.bam")
            bam.parent.mkdir(parents=True)
            bam.write_bytes(b"synthetic, never parsed as BAM")
            Path(str(bam) + ".bai").touch()
        resources = dict(bind_roots=[str(self.root)])
        for key in ("reference_fasta", "samtools_image", "sage_image", "pave_image", "hmf_resources", "vep_image", "vep_cache", "vep_plugins", "pvactools_image", "genes_of_interest"):
            path = self.root / key
            path.touch()
            resources[key] = str(path)
        self.config = dict(schema_version=planner.SCHEMA, resources=resources, cases=[self.case])

    def plan(self):
        return planner.build_plan(self.config)

    def test_both_extensions_ready_independently(self):
        plan = self.plan()
        self.assertEqual([plan["cases"][0][k]["status"] for k in ("blood_ch", "rna_neoantigens")], ["READY", "READY"])
        self.case["options"]["blood_ch"] = False
        self.assertEqual(self.plan()["cases"][0]["blood_ch"]["status"], "DISABLED")
        self.assertEqual(self.plan()["cases"][0]["rna_neoantigens"]["status"], "READY")
        self.case["options"]["rna_neoantigens"] = False
        self.assertEqual(self.plan()["cases"][0]["rna_neoantigens"]["status"], "DISABLED")

    def test_blood_is_never_inferred_from_normal_id(self):
        for specimen in ("", "unknown", "saliva", "normal"):
            self.case["laboratory"]["normal_specimen_type"] = specimen
            self.assertEqual(self.plan()["cases"][0]["blood_ch"]["status"], "INELIGIBLE")

    def test_recorded_rna_requires_actual_native_rna(self):
        self.case["laboratory"]["tumour_rna_available"] = False
        self.assertEqual(self.plan()["cases"][0]["rna_neoantigens"]["status"], "INELIGIBLE")
        self.case["laboratory"]["tumour_rna_available"] = True
        self.case["rna_sample"] = "ANOTHER_RNA"
        self.assertEqual(self.plan()["cases"][0]["rna_neoantigens"]["status"], "MISSING_INPUTS")

    def test_unknown_lab_source_is_not_eligibility(self):
        del self.case["laboratory"]["source"]
        plan = self.plan()
        self.assertEqual([plan["cases"][0][k]["status"] for k in ("blood_ch", "rna_neoantigens")], ["INELIGIBLE", "INELIGIBLE"])

    def test_missing_reference_does_not_become_negative(self):
        Path(self.config["resources"]["reference_fasta"]).unlink()
        self.assertEqual(self.plan()["cases"][0]["blood_ch"]["status"], "MISSING_INPUTS")

    def test_mismatched_dna_binding_fails(self):
        self.case["normal_sample"] = "WRONG_NORMAL"
        with self.assertRaisesRegex(ValueError, "identities disagree"):
            self.plan()

    def test_custom_reference_not_implicitly_admitted(self):
        self.case["native_genome"] = "GRCh38_d1_vd1"
        with self.assertRaisesRegex(ValueError, "native GRCh38_hmf"):
            self.plan()

    def test_resources_increase_with_actual_workload(self):
        for stage in ("blood_ch", "neo_prepare", "neo_predict"):
            low, high = planner.resource_allowance(stage, 1), planner.resource_allowance(stage, 1200)
            self.assertGreater(high["minutes"], low["minutes"])
            self.assertEqual(high["cpus"], low["cpus"])

    def test_selected_missing_task_is_incomplete(self):
        status = planner.collect_status(self.plan(), self.root / "outputs")
        self.assertEqual(status["cases"][0]["blood_ch"]["status"], "INCOMPLETE")

    def test_failed_and_disabled_remain_distinct(self):
        self.case["options"]["rna_neoantigens"] = False
        out = self.root / "outputs/SYNTHETIC/blood_ch/status.json"
        planner.write_json(out, dict(status="FAILED", reason="caller failed"))
        status = planner.collect_status(self.plan(), self.root / "outputs")
        self.assertEqual(status["cases"][0]["blood_ch"]["status"], "FAILED")
        self.assertEqual(status["cases"][0]["rna_neoantigens"]["status"], "DISABLED")

    def test_typing_failure_does_not_run_prediction(self):
        path = self.root / "outputs/SYNTHETIC/neo_export/status.json"
        planner.write_json(path, dict(status="NOT_ASSESSED", reason="failed HLA"))
        status = planner.collect_status(self.plan(), self.root / "outputs")
        self.assertEqual(status["cases"][0]["rna_neoantigens"]["status"], "NOT_ASSESSED")

    def test_failed_command_has_its_own_receipt(self):
        out = self.root / "chip"
        with patch.object(tasks, "chip", side_effect=RuntimeError("synthetic tool failure")):
            with self.assertRaisesRegex(RuntimeError, "synthetic tool failure"):
                tasks.execute(self.plan(), "SYNTHETIC", "blood_ch", out, 8)
        receipt = json.loads((out / "status.json").read_text())
        self.assertEqual(receipt["status"], "FAILED")
        self.assertTrue(self.case_root.is_dir())

    def test_unpublished_failure_is_recovered_from_trace(self):
        trace = self.root / "trace.tsv"
        trace.write_text("process\ttag\tstatus\texit\tworkdir\nOPTIONAL_EXTENSIONS:BLOOD_CH\tSYNTHETIC\tFAILED\t137\t" + str(self.root / "failed_task") + "\n")
        planner.write_json(self.root / "failed_task/blood_ch/status.json", dict(status="FAILED", reason="synthetic memory failure"))
        status = planner.collect_status(self.plan(), self.root / "outputs", trace)
        self.assertEqual(status["cases"][0]["blood_ch"]["status"], "FAILED")
        self.assertEqual(status["cases"][0]["blood_ch"]["reason"], "synthetic memory failure")
        self.assertEqual(status["cases"][0]["blood_ch"]["task_exit"], "137")

    def test_chip_uses_native_blood_binding_and_explicit_resources(self):
        class Capture:
            def run(self, argv, env):
                self.argv, self.env = argv, env
        calls = Capture()
        tasks.chip(self.case, self.config["resources"], self.root, 8, calls)
        self.assertEqual(calls.argv[3], self.case["normal_sample"])
        self.assertEqual(calls.env["PAVE_SIF"], self.config["resources"]["pave_image"])
        self.assertEqual(calls.env["THREADS"], "8")

    def test_chip_publication_never_copies_staged_input_bams(self):
        out = self.root / "chip_output"
        source = self.root / "original.bam"
        source.write_bytes(b"preserve original")
        stage = out / "stage"
        stage.mkdir(parents=True)
        (stage / "input.bam").symlink_to(source)
        class Capture:
            def run(self, argv, env):
                pass
        tasks.chip(self.case, self.config["resources"], out, 8, Capture())
        self.assertFalse(stage.exists())
        self.assertEqual(source.read_bytes(), b"preserve original")

    def test_unknown_file_in_input_stage_is_preserved(self):
        out = self.root / "chip_output"
        stage = out / "stage"
        stage.mkdir(parents=True)
        source = stage / "unknown.bam"
        source.write_bytes(b"do not delete")
        class Capture:
            def run(self, argv, env):
                pass
        with self.assertRaisesRegex(ValueError, "Unexpected non-link"):
            tasks.chip(self.case, self.config["resources"], out, 8, Capture())
        self.assertTrue(source.is_file())

    def test_native_head_only_runs_extensions_with_explicit_config(self):
        text = (ROOT / "workflow/run_full_native.sbatch").read_text()
        self.assertLess(text.index('bash "$PACKAGE/workflow/run_oncoanalyser_frontenac.sh"'), text.index('if [[ -f "$PACKAGE/optional_extensions.json" ]]'))
        self.assertIn('native_workflow_completed.json', text)
        self.assertIn('run_optional_extensions.sh', text)

    def test_prediction_keeps_native_annotations_and_local_predictors(self):
        class Capture:
            def __init__(self):
                self.calls = []
            def run(self, argv, env=None):
                self.calls.append([str(v) for v in argv])
            def container(self, image, *argv):
                self.calls.append([image] + [str(v) for v in argv])
        prepared = self.root / "prepared"
        prepared.mkdir()
        planner.write_json(prepared / "SYNTHETIC.interface_manifest.json", dict(hla=dict(alleles=["HLA-A*02:01"])))
        calls = Capture()
        tasks.neo_predict(self.case, self.config["resources"], self.root, 8, calls, prepared)
        vep = [c for c in calls.calls if c[0] == "vep_image"]
        self.assertEqual(len(vep), 2)
        for command in vep:
            self.assertIn("--use_given_ref", command)
            self.assertIn("Wildtype", command)
            self.assertIn("Frameshift", command)
            self.assertIn("--tsl", command)
            self.assertIn("--biotype", command)
        prediction = next(c for c in calls.calls if c[:3] == ["pvactools_image", "pvacseq", "run"])
        self.assertIn("--genes-of-interest-file", prediction)
        self.assertIn("MHCflurry", prediction)
        self.assertIn("MHCflurryEL", prediction)
        self.assertIn("MHCnuggetsI", prediction)
        self.assertNotIn("NetMHCpan", prediction)
        self.assertIn("--normal-sample-name", prediction)
        self.assertIn(self.case["normal_sample"], prediction)


if __name__ == "__main__":
    unittest.main()
