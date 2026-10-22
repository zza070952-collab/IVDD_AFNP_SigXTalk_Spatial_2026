"""Technical structure tests only; never substitute these for research outputs."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts/human" / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


class HumanTests(unittest.TestCase):
    def test_frozen_members_and_original_order(self):
        folder = ROOT / "resources/frozen_k5"
        lock = json.loads((folder / "PROGRAM_LOCK.json").read_text())
        self.assertEqual(hashlib.sha256((folder / "PROGRAM_LOCK.json").read_bytes()).hexdigest(),
                         "12905c6475cd4f3dd07c660efeee15e60d59f592c8f74c6cc60202f549110da2")
        table = pd.read_csv(folder / "FROZEN_K5_PROGRAMS.tsv", sep="\t")
        self.assertEqual(len(table), 125)
        for name, genes in lock["programs"].items():
            self.assertEqual(table[table.program == name].sort_values("rank").gene.tolist(), genes)

    def test_sparse_decimal_parse_and_conservative_identity(self):
        parser = module("human_parser", "parse_tu.py")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "structural_fixture.tsv.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("gene\tC1\tC2\nACAN\t1.25\t0\nCOL2A1\t1\t0\nPTPRC\t0\t2\nLST1\t0\t1\n")
            matrix, genes, cells, info = parser.load_sparse_tsv(path, "TEST_ONLY")
            self.assertEqual(matrix.shape, (4, 2))
            self.assertEqual(matrix[0, 0], 1.25)
            self.assertFalse(info["all_submitted_values_integer"])
            self.assertEqual(parser.annotate(matrix, genes).identity.tolist(),
                             ["NP_matrix_marker_supported", "immune_marker_supported"])

    def test_alias_requires_one_present_feature(self):
        parser = module("human_alias", "parse_tu.py")
        aliases = [{"original_gene": "OLD", "mapped_gene": "NEW1"},
                   {"original_gene": "OLD", "mapped_gene": "NEW2"}]
        self.assertEqual(parser.match_genes(["OLD"], ["NEW1"], aliases)[0][1], "NEW1")
        self.assertEqual(parser.match_genes(["OLD"], ["NEW1", "NEW2"], aliases)[0][1], "")
        self.assertEqual(parser.match_genes(["OLD"], ["OTHER"], aliases)[0][1], "")

    def test_summary_does_not_rename_aucell_or_fill_missing(self):
        summary = module("human_summary", "summarize.py")
        metadata = pd.DataFrame({"cell_id": ["a", "b"], "donor_id": ["d1", "d2"],
                                 "gsm": ["g1", "g2"], "identity": ["test", "test"]})
        scores = pd.DataFrame({"cell_id": ["a", "b"], **{f"GEP{i}_top25_AUCell": [0.2, np.nan if i == 3 else 0.4] for i in range(1, 6)}})
        actual, _ = summary.summarize_cells(metadata, scores, "Tu_TEST_ONLY", pd.Series([True, True]), "AUCell_decoupler_2.2.0", "AUCell")
        self.assertEqual(set(actual.score_method), {"AUCell_decoupler_2.2.0"})
        empty = actual[(actual.program == "GEP3") & (actual.donor_id == "d2")].iloc[0]
        self.assertEqual(empty.n_valid_scores, 0)
        self.assertTrue(pd.isna(empty.score_median))


if __name__ == "__main__":
    unittest.main()
