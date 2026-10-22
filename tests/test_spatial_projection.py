"""Small synthetic technical fixture only; never a scientific result table."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
from scipy import sparse

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'

class SparseProjectionTest(unittest.TestCase):
    def test_fixed_available_member_mean_and_zero_detection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            folder = root/'prepared48/GSM_TECHNICAL_FIXTURE'
            folder.mkdir(parents=True)
            counts = np.array([[5, 0, 95], [0, 20, 80]], dtype=np.int64)
            sparse.save_npz(folder/'counts48_symbols_qc.npz', sparse.csr_matrix(counts))
            obs = pd.DataFrame({'spot_id':['s1','s2'], 'capture_id':[folder.name]*2,
                'genotype':['TECHNICAL']*2, 'x_fullres':[1.,2.], 'y_fullres':[3.,4.],
                'n_counts':[100,100], 'n_genes':[2,2]})
            obs.to_csv(folder/'BIN48_METADATA_QC.tsv.gz', sep='\t', index=False)
            pd.DataFrame({'mouse_gene_symbol':['Genea','Geneb','Other']}).to_csv(
                folder/'MOUSE_SYMBOL_FEATURES.tsv', sep='\t', index=False)
            roi = obs.assign(disc_id='technical_disc', histology_core_region=['NP_CORE','AF_CORE'],
                             anatomical_region=['NP_CORE','AF_CORE'])
            roi.to_csv(root/'roi.tsv', sep='\t', index=False)
            frozen = pd.DataFrame({'program':['P','P'], 'rank':[1,2], 'gene':['GENEA','GENEB']})
            frozen.to_csv(root/'frozen.tsv', sep='\t', index=False)
            pd.DataFrame({'capture_id':[folder.name]*2, 'human_gene':['GENEA','GENEB'],
                'matched_native_mouse_symbol':['Genea','Geneb'], 'mapping_status':['MGI_ONE_TO_ONE']*2,
                'n_detected_bins':[1,1]}).to_csv(root/'coverage.tsv', sep='\t', index=False)
            result = subprocess.run([sys.executable, str(SCRIPTS/'spatial_programs.py'),
                '--prepared48', str(root/'prepared48'), '--roi-metadata', str(root/'roi.tsv'),
                '--coverage', str(root/'coverage.tsv'), '--programs-file', str(root/'frozen.tsv'),
                '--programs','P','--candidates','GENEA','GENEB','--out',str(root/'output')],
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            values = pd.read_csv(root/'output/SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz', sep='\t')
            expected = np.log1p(counts[:,:2]*100).mean(axis=1)
            np.testing.assert_allclose(values.P_mean_log1p_CP10k, expected, atol=1e-12, rtol=0)
            genes = pd.read_csv(root/'output/CANDIDATE_GENES_BY_REGION.tsv', sep='\t')
            absent = genes.loc[genes.human_gene.eq('GENEA') & genes.region.eq('AF_CORE')].iloc[0]
            self.assertEqual(absent.n_detected_bins, 0)
            self.assertEqual(absent.mean_log1p_CP10k, 0)

if __name__ == '__main__':
    unittest.main()
