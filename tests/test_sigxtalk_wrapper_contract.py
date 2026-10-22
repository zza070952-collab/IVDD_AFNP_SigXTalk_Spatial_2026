"""Structural tests only: temporary mock CSV values are not research results."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

BASE=Path(__file__).resolve().parents[1]
SCRIPT=BASE/'run_sigxtalk_pipeline.py'
if not SCRIPT.exists():SCRIPT=BASE/'scripts/run_sigxtalk_pipeline.py'
spec=importlib.util.spec_from_file_location('wrapper',SCRIPT)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class ContractTests(unittest.TestCase):
    def test_nonfinite_prediction_is_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'structural_fixture.csv'
            p.write_text('Receptor,TF,TG,pred_label\nMOCK_R,MOCK_S,MOCK_T,nan\n')
            with self.assertRaises(ValueError):module.count_csv(p)

    def test_cli_output_csv_does_not_imply_R_closure(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'structural_fixture.csv'
            p.write_text('Receptor,TF,TG,pred_label\nMOCK_R,MOCK_S,MOCK_T,0.8\nMOCK_R,MOCK_S,MOCK_U,0.5\n')
            self.assertEqual(module.count_csv(p),(2,1))
            self.assertFalse((Path(temp)/'NEW_R_PYTHON_R_COMPLETION.json').exists())

    def test_pinned_source_mismatch_is_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);raw=b'nonbiological structural fixture\n'
            (p/'own_fixture.txt').write_bytes(raw)
            blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
            lock={'repository':'fixture_not_upstream','commit':'test_only',
                  'git_blob_sha1':{'own_fixture.txt':blob}}
            lp=p/'lock.json';lp.write_text(json.dumps(lock))
            self.assertEqual(module.verify_upstream(p,lp)['source_file_identity_check'],'PASS')
            (p/'own_fixture.txt').write_text('changed')
            with self.assertRaises(ValueError):module.verify_upstream(p,lp)


if __name__=='__main__':unittest.main()
