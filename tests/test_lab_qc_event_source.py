import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

SPEC=importlib.util.spec_from_file_location('event_source',Path(__file__).parents[1]/'workflow/prepare_lab_qc_event_source.py')
tool=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(tool)


class EventSourceTests(unittest.TestCase):
    def test_native_source_is_preserved_and_hook_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'native';base.mkdir()
            subprocess.run(['git','init','--quiet',str(base)],check=True)
            for key,value in [('user.name','Synthetic Reviewer'),('user.email','synthetic@example.invalid')]:
                subprocess.run(['git','-C',str(base),'config',key,value],check=True)
            (base/'lib').mkdir();original='workflow {\n'+tool.ANCHOR+'}\n';(base/'main.nf').write_text(original)
            subprocess.run(['git','-C',str(base),'add','.'],check=True)
            subprocess.run(['git','-C',str(base),'commit','--quiet','-m','Synthetic native fixture'],check=True)
            result=tool.prepare(base,root/'sources');again=tool.prepare(base,root/'sources')
            self.assertEqual(result,again);self.assertEqual((base/'main.nf').read_text(),original)
            source=Path(result['source']);self.assertIn("channel.topic('bamtools_metrics_dir')",(source/'main.nf').read_text())
            self.assertEqual(subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True),'')
            self.assertEqual(tool.prepare(source,root/'sources'),result)
