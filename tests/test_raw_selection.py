"""Tiny structural tests; not biological results or a substitute for a raw run."""
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import io, sparse

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("public_raw_preparer", ROOT / "scripts/prepare_sigxtalk_raw.py")
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)


class RawSelectionTests(unittest.TestCase):
    def test_distributed_manifest_is_eight_verified_pairs(self):
        table = prep.validate_manifest(ROOT / "resources/metadata/SIGXTALK_GSE230808_SOURCE_MANIFEST.tsv")
        self.assertEqual(len(table), 16)
        self.assertEqual(table.donor_id.nunique(), 8)

    def test_seeded_selection_is_recomputed_and_ignores_input_order(self):
        cells = ["STRUCT_NP:" + str(i) for i in range(8)]
        metadata = pd.DataFrame({"cell_id": cells, "donor_id": "STRUCT_ONLY", "gsm": "STRUCT_NP",
            "high_confidence_NP": True, "n_counts": 10, "n_features": 3}).set_index("cell_id", drop=False)
        scores = pd.DataFrame({"cell_id": cells, "GEP2_top25_UCell": np.arange(8) / 10}).set_index("cell_id", drop=False)
        af_cells = ["STRUCT_AF:" + str(i) for i in range(8)]
        af = pd.DataFrame({"cell_id": af_cells, "donor_id": "STRUCT_ONLY", "gsm": "STRUCT_AF",
            "broad_stromal": True, "n_counts": 10, "n_features": 3}).set_index("cell_id", drop=False)
        first, eligibility = prep.select_cells(metadata, scores, af, minimum=2, cap=3, expected_eligible=1)
        second, _ = prep.select_cells(metadata.iloc[::-1], scores.iloc[::-1], af.iloc[::-1], minimum=2, cap=3, expected_eligible=1)
        self.assertEqual(first.index.tolist(), second.index.tolist())
        self.assertEqual(len(first), 7)
        self.assertEqual(set(first.bridge_group), {"NP_GEP2_low", "NP_GEP2_high", "AF_stromal"})
        self.assertTrue(eligibility.eligible.all())

    def test_real_sparse_values_and_duplicate_qc_audit_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, out, audit = root / "raw", root / "input", root / "audit"
            for folder in [raw, out, audit]:
                folder.mkdir()
            manifest, rows = [], []
            for gsm, tissue, values, expected_features in [
                ("STRUCT_NP", "NP", [1, 2, 0], 1), ("STRUCT_AF", "AF", [4, 0, 1], 2)]:
                prefix = raw / (gsm + "_test_")
                with gzip.open(str(prefix) + "features.tsv.gz", "wt") as handle:
                    handle.write("ENS1\tDUP\tGene Expression\nENS2\tDUP\tGene Expression\nENS3\tMT-Z\tGene Expression\n")
                with gzip.open(str(prefix) + "barcodes.tsv.gz", "wt") as handle:
                    handle.write("BARCODE\n")
                with gzip.open(str(prefix) + "matrix.mtx.gz", "wb") as handle:
                    io.mmwrite(handle, sparse.coo_matrix(np.array(values).reshape(3, 1)), field="integer")
                manifest.append({"gsm": gsm, "donor_id": "STRUCT_ONLY", "tissue": tissue})
                rows.append({"cell_id": gsm + ":BARCODE", "gsm": gsm, "donor_id": "STRUCT_ONLY",
                    "tissue": tissue, "bridge_group": "NP_GEP2_high" if tissue == "NP" else "AF_stromal",
                    "n_counts": sum(values), "n_features": expected_features,
                    "GEP2": .4 if tissue == "NP" else np.nan})
            meta = pd.DataFrame(rows).set_index("cell_id", drop=False)
            result = prep.reconstruct_counts(meta, pd.DataFrame(manifest), raw, out, audit)
            self.assertEqual(result["shape_genes_by_cells"], [3, 2])
            self.assertEqual(result["raw_UMI_equal_to_cached"], 2)
            self.assertEqual(result["historical_symbol_policy_features_equal_to_cached"], 2)
            self.assertEqual(result["raw_features_equal_to_cached"], 1)
            with gzip.open(out / "counts.mtx.gz", "rb") as handle:
                actual = io.mmread(handle).toarray()
            np.testing.assert_array_equal(actual, [[1, 4], [2, 0], [0, 1]])
            self.assertEqual((out / "genes.tsv").read_text().splitlines(), ["DUP", "DUP-1", "MT-Z"])

    def test_duplicate_name_suffix_collision_is_explicit(self):
        self.assertEqual(prep.unique_symbols(["A", "A", "A-1", "A"]).tolist(), ["A", "A-2", "A-1", "A-3"])


if __name__ == "__main__":
    unittest.main()
