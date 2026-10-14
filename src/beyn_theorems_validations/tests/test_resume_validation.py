"""Resume boundaries and checkpoint failure controls; no numerical sweep."""
import sys,json,tempfile,unittest
from pathlib import Path
from dataclasses import replace,fields
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'runners'))
import resume_validation as r

class ResumeTests(unittest.TestCase):
    def test_exact_manifest_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(r.f.CONFIG,output_directory=tmp)
            r.f.write_run_manifest(Path(tmp),cfg,'interrupted')
            loaded,_=r.load_config(r.f,tmp)
            self.assertEqual(loaded,cfg)
            p=Path(tmp)/'run_manifest.json';data=json.loads(p.read_text());data['configuration_fingerprint']='wrong';p.write_text(json.dumps(data))
            with self.assertRaises(ValueError):r.load_config(r.f,tmp)

    def test_atomic_paper_case_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(r.f.CONFIG,output_directory=tmp)
            row=r.f._paper_empty_row(.1,.6,cfg,'ambiguous','test control')
            with patch.object(r.f,'evaluate_paper_sweep_geometry',return_value=row) as calculate:
                first=r.paper_checkpoint_worker((.6,.1,cfg));second=r.paper_checkpoint_worker((.6,.1,cfg))
                self.assertEqual(calculate.call_count,1);self.assertEqual(first.status,second.status)

    def test_error_not_checkpointed(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(r.f.CONFIG,output_directory=tmp)
            row=r.f._paper_empty_row(.1,.6,cfg,'error','injected failure')
            with patch.object(r.f,'evaluate_paper_sweep_geometry',return_value=row):
                with self.assertRaises(ValueError):r.paper_checkpoint_worker((.6,.1,cfg))
            self.assertFalse(list(Path(tmp).rglob('completed.json')))

    def test_case_geometry_mismatch_rejected(self):
        row=r.f._paper_empty_row(.1,.6,r.f.CONFIG,'ambiguous','control')
        self.assertFalse(r.valid_paper(row,.5,.1,r.f.CONFIG))

if __name__=='__main__':unittest.main()
