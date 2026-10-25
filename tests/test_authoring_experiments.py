"""Small synthetic regression cases; no observed experiment values are fixtures."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
from scipy import stats

BASE=Path(__file__).resolve().parents[1]/'scripts/authoring/v4/experiments'
def load(name):
    spec=importlib.util.spec_from_file_location('authoring_'+name,BASE/(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

class SequencePairing(unittest.TestCase):
    def test_strands_are_not_taken_from_next_paragraph(self):
        parse=load('sequence_parser').parse_sequence_pairs
        paragraphs=['siExampleA S ACGUACGUdTdT\nAS UGCAUGCAdTdT','siExampleB S AAAACCCCdTdT\nAS GGGGUUUUdTdT']
        rows=parse(paragraphs)
        self.assertEqual(rows[0]['reverse_5to3'],'UGCAUGCAdTdT')
        self.assertEqual(rows[1]['reverse_5to3'],'GGGGUUUUdTdT')
        self.assertEqual([r['source_paragraph'] for r in rows],[1,2])
    def test_incomplete_paragraph_is_not_completed_by_neighbour(self):
        rows=load('sequence_parser').parse_sequence_pairs(['siExampleA S ACGUACGUdTdT','AS UGCAUGCAdTdT'])
        self.assertEqual(rows,[])
    def test_primers_pair_with_same_name_and_metadata(self):
        rows=load('sequence_parser').parse_sequence_pairs(['ExampleA（123c1，8 bp）\nF ACGTACGT R TGCATGCA','ExampleB(456c2, 9 bp)\nF AAAACCCC R GGGGTTTT'])
        self.assertEqual([(r['record'],r['primerbank_id'],r['amplicon_bp']) for r in rows],[('ExampleA','123c1',8),('ExampleB','456c2',9)])
        self.assertEqual(rows[0]['forward_5to3'],'ACGTACGT')
    def test_duplicate_name_fails_instead_of_silent_overwrite(self):
        with self.assertRaises(ValueError):load('sequence_parser').parse_sequence_pairs(['siExample S ACGUdTdT AS UGCAdTdT']*2)

class BatchStatistics(unittest.TestCase):
    def test_primary_ci_exponentiation_reverses_bounds(self):
        values=np.array([-1.0,0.0,1.0,2.0])
        result=load('analyze').summarize(values)
        mean=values.mean();se=values.std(ddof=1)/2;radius=stats.t.ppf(.975,3)*se
        self.assertEqual(result['n_batches'],4)
        self.assertAlmostEqual(result['deltaCt_CI95_low'],mean-radius)
        self.assertAlmostEqual(result['fold_CI95_low'],2**(-(mean+radius)))
        self.assertAlmostEqual(result['fold_CI95_high'],2**(-(mean-radius)))
    def test_adjustments_keep_original_row_order(self):
        module=load('analyze')
        np.testing.assert_allclose(module.adjust([.03,.001,.04],'Holm'),[.06,.003,.06])
        np.testing.assert_allclose(module.adjust([.03,.001,.04],'BH'),[.04,.003,.04])

if __name__=='__main__':unittest.main()
