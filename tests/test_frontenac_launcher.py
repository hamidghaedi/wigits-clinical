"""Check configuration discovery before module loading or any scheduler work."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class LauncherConfigTests(unittest.TestCase):
    def test_checkout_and_slurm_spool_paths(self):
        launcher = Path(__file__).resolve().parents[1] / "workflow/run_oncoanalyser_frontenac.sh"
        prefix = launcher.read_text().split("source /cvmfs/", 1)[0]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workflow = root / "checkout/workflow"
            (workflow / "conf").mkdir(parents=True)
            (workflow / "conf/frontenac_paths.env").write_text("CONFIG_FOUND=yes\n")
            spool = root / "spool"
            spool.mkdir()
            for name, script_dir, extra, expected in [
                ("direct", workflow, {}, 0),
                ("sbatch", spool, {"SLURM_SUBMIT_DIR": str(workflow.parent)}, 0),
                ("override", spool, {"WIGITS_WORKFLOW_DIR": str(workflow)}, 0),
                ("missing", spool, {"SLURM_SUBMIT_DIR": str(root)}, 2),
            ]:
                with self.subTest(name=name):
                    script = script_dir / "slurm_script"
                    script.write_text(prefix + '\nprintf "%s" "$CONFIG_FOUND"\n')
                    env = {k: v for k, v in os.environ.items()
                           if k not in {"SLURM_SUBMIT_DIR", "WIGITS_WORKFLOW_DIR"}}
                    result = subprocess.run(["bash", str(script)], env={**env, **extra},
                                            text=True, capture_output=True)
                    self.assertEqual(result.returncode, expected, result.stderr)
                    if expected == 0:
                        self.assertEqual(result.stdout, "yes")
                    else:
                        self.assertIn("workflow configuration not found", result.stderr)


if __name__ == "__main__":
    unittest.main()
