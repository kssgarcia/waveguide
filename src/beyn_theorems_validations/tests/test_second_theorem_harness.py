"""Small execution/geometry regressions; no theorem sweeps."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import second_theorem_X as x
import second_theorem_Y as y
import execution_harness as h


class HarnessTests(unittest.TestCase):
    def test_unique_runs_preserve_previous_results(self):
        for module in (x,y):
            with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
                cfg=replace(module.CONFIG,output_directory=tmp)
                a=h.new_run(module,cfg);p=Path(a.output_directory)/'summary.csv';p.write_text('previous')
                b=h.new_run(module,cfg)
                self.assertNotEqual(a.output_directory,b.output_directory)
                self.assertEqual(p.read_text(),'previous')
                pointer=json.loads((Path(tmp)/'latest_run.json').read_text())
                self.assertEqual(pointer['directory'],b.output_directory)

    def test_checkpoint_inputs_config_tamper_and_failure(self):
        for module in (x,y):
            with tempfile.TemporaryDirectory() as tmp:
                cfg=replace(module.CONFIG,output_directory=tmp)
                result=(np.array([1.+2.j]),{'resolved':False,'nonexistence_supported':False,'residual':float('nan')})
                with patch.object(module,'validate_epsilon',return_value=result) as calculate:
                    h.case(module,cfg,'test',(.1,))
                    cached=h.case(module,cfg,'test',(.1,))
                    self.assertEqual(calculate.call_count,1)
                    np.testing.assert_array_equal(cached[0],result[0])
                    self.assertFalse(cached[1]['nonexistence_supported'])
                    h.case(module,cfg,'test',(.2,));self.assertEqual(calculate.call_count,2)
                    h.case(module,replace(cfg,harmonic_order=cfg.harmonic_order+1),'test',(.2,))
                    self.assertEqual(calculate.call_count,3)
                    path=Path(tmp)/'epsilon_cases/test/completed.json'
                    saved=json.loads(path.read_text());saved['payload']+=' '
                    path.write_text(json.dumps(saved))
                    h.case(module,cfg,'test',(.2,));self.assertEqual(calculate.call_count,4)
                path.unlink()
                with patch.object(module,'validate_epsilon',side_effect=RuntimeError('failed refinement')):
                    with self.assertRaises(RuntimeError):h.case(module,cfg,'test',(.2,))
                self.assertFalse(path.exists())

    def test_y_continuation_seed_is_part_of_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(y.CONFIG,output_directory=tmp)
            with patch.object(y,'validate_epsilon',return_value=(False,)) as calculate:
                h.case(y,cfg,'test',(.1,),{'previous_mode':None})
                h.case(y,cfg,'test',(.1,),{'previous_mode':'changed seed'})
                self.assertEqual(calculate.call_count,2)

    def test_failure_and_interrupt_are_not_completed(self):
        for exception in (RuntimeError('failure'),KeyboardInterrupt()):
            with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
                cfg=replace(x.CONFIG,output_directory=tmp)
                with self.assertRaises(type(exception)):
                    h.run_entrypoint(x,cfg,lambda cfg: (_ for _ in ()).throw(exception))
                run=next(p for p in Path(tmp).iterdir() if p.is_dir())
                status=json.loads((run/'run_manifest.json').read_text())['stage']
                self.assertEqual(status,'interrupted' if isinstance(exception,KeyboardInterrupt) else 'failed')

    def test_cached_geometry_preserves_full_and_odd_matrices(self):
        for module in (x,y):
            cfg=module.CONFIG
            module._kernel_geometry.cache_clear()
            def direct(*args):return module.obstacle_geometry(*args)
            for k in (2.,2.+.01j):
                args=(k,.05,8,cfg) if module is x else (k,.05,-.002,8,cfg)
                before=module.assemble_full_matrix(*args)
                with patch.object(module,'_kernel_geometry',side_effect=direct):
                    after=module.assemble_full_matrix(*args)
                np.testing.assert_array_equal(before,after)
            self.assertGreater(module._kernel_geometry.cache_info().hits,0)
            if module is x:
                odd=module.assemble_odd_matrix(2.+.01j,.05,8,cfg)
                Q=module.odd_lift_matrix(8)
                np.testing.assert_array_equal(odd,Q.conj().T@module.assemble_full_matrix(2.+.01j,.05,8,cfg)@Q)


if __name__=='__main__':unittest.main()
