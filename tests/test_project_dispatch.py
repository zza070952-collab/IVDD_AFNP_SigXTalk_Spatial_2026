"""Technical receipt/guard fixtures only; not biological observations."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import subprocess

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "scripts"))
import project_common as common
import run_project as driver


class DispatcherTests(unittest.TestCase):
    def test_incomplete_depth_is_not_redrawn(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "DEPTH_ANALYSIS_STATUS.json").write_text('{"status":"RUNNING"}')
            config, roi = root / "depth.json", root / "roi.json"
            config.write_text('{"seed":20261022}')
            roi.write_text('{}')
            values = {"depth_output": root, "depth_config": config, "roi_config": roi}
            with self.assertRaises(RuntimeError):
                driver.verify_completed_depth(values, root)

    def test_completed_depth_reuses_exact_recorded_draw_and_protected_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "output"
            output.mkdir()
            config, roi, old = root / "depth.json", root / "roi.json", root / "TECHNICAL_PROTECTED_SOURCE.tsv"
            config.write_text('{"seed":20261022}')
            roi.write_text('{}')
            old.write_text('structure_only\n')
            status = {"status":"DONE_SINGLE_PREDECLARED_DEPTH_CHECK", "target_total_UMI":100,
                "single_draw_only":True,"protected_old_inputs_unchanged":True,"random_seed":20261022,
                "config_sha256":common.sha(config),"ROI_sha256":common.sha(roi)}
            common.dump(status, output / "DEPTH_ANALYSIS_STATUS.json")
            for name in driver.DEPTH_FILES:
                (output / name).write_text('TECHNICAL_FIXTURE_NOT_A_RESEARCH_TABLE\n')
            (output / "PRIVATE_FROZEN_INPUT_HASH_AUDIT.tsv").write_text(
                'path\tsha256_before\tsha256_after\n'+str(old)+'\t'+common.sha(old)+'\t'+common.sha(old)+'\n')
            values={"depth_output":output,"depth_config":config,"roi_config":roi}
            actual=driver.verify_completed_depth(values,root)
            self.assertFalse(actual["single_draw_repeated"])
            self.assertTrue(actual["status"].startswith("SKIPPED_VERIFIED"))
            old.write_text('changed\n')
            with self.assertRaises(ValueError):
                driver.verify_completed_depth(values,root)

    def test_correct_private_B_configuration_is_not_rewritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory).resolve()
            path=root/'actual.json'
            path.write_text(json.dumps({'project_root':str(root)},separators=(',',':')))
            before=path.read_bytes()
            used,_=common.resolved_b_config({'_root':root,'modules':{'sigxtalk_config':str(path)}})
            self.assertEqual(used,path)
            self.assertEqual(path.read_bytes(),before)

    def test_prediction_scores_without_trained_input_hash_cannot_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ['counts.mtx.gz','genes.tsv','meta.tsv','pairing.tsv']:
                (root/name).write_bytes(b'TECHNICAL_ONLY')
            values={'input_dir':'.','pairing_tsv':'pairing.tsv'}
            provenance={'input_kind':'real_GSE230808_raw_counts','new_input_file_hashes':{}}
            with self.assertRaises(ValueError):
                driver.check_b_input_hashes({'_root':root},provenance,values)

    def test_existing_external_pairing_and_thresholds_are_read_only_inputs(self):
        import pandas as pd
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory).resolve()
            root,old=folder/'new',folder/'old'
            root.mkdir();old.mkdir()
            pairs=pd.DataFrame([{'donor_id':f'STRUCT_{i}','AF_gsm':f'AF{i}','NP_gsm':f'NP{i}'} for i in range(6)])
            pairing,thresholds=old/'pairs.tsv',old/'thresholds.tsv'
            pairs.to_csv(pairing,sep='\t',index=False)
            q=pairs[['donor_id']].assign(q25=.2,q75=.8)
            q.to_csv(thresholds,sep='\t',index=False)
            manifest=folder/'manifest.tsv'
            pd.concat([pairs[['donor_id','AF_gsm']].rename(columns={'AF_gsm':'gsm'}).assign(tissue='AF'),
                       pairs[['donor_id','NP_gsm']].rename(columns={'NP_gsm':'gsm'}).assign(tissue='NP')]).to_csv(manifest,sep='\t',index=False)
            selected=root/'input';selected.mkdir()
            q.to_csv(selected/'RECOMPUTED_DONOR_ELIGIBILITY.tsv',sep='\t',index=False)
            bcfg=folder/'b.json';bcfg.write_text(json.dumps({'project_root':str(root),'module_dir':'module',
                'pairing_tsv':str(pairing),'state_thresholds_tsv':str(thresholds)}))
            rawcfg=folder/'raw.json';rawcfg.write_text(json.dumps({'raw_selection':{'sample_manifest':str(manifest),'output_dir':'input'}}))
            before=[pairing.read_bytes(),thresholds.read_bytes()]
            driver.materialize_pairing({'_root':root,'modules':{'sigxtalk_config':str(bcfg)}},rawcfg)
            self.assertEqual([pairing.read_bytes(),thresholds.read_bytes()],before)

    def test_cached_downloader_records_no_transfer_without_unused_raw_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory).resolve()
            cfg=root/'empty.json';cfg.write_text('{}')
            completed=subprocess.run([sys.executable,str(BASE/'scripts/00_download_public_data.py'),
                '--config',str(cfg),'--project-root',str(root),'--from-cached-inputs'],capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stdout+completed.stderr)
            proof=json.loads((root/'00_audit/PUBLIC_DOWNLOAD_CACHE_SKIP.json').read_text())
            self.assertEqual(proof['raw_network_transfers'],0)
            self.assertTrue(proof['status'].startswith('SKIPPED'))


if __name__ == '__main__':
    unittest.main()
