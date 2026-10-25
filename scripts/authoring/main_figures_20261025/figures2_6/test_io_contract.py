"""Synthetic-only tests of public path and byte-copy contracts."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import copy_retained_figures as copy
from _io_guard import independent_output, reject_output_links


class IOContractTests(unittest.TestCase):
    def test_disjoint_trees(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            independent_output(root/'input', root/'output')
            for src, out in [(root, root), (root, root/'child'),
                             (root/'child', root)]:
                with self.assertRaises(ValueError):
                    independent_output(src, out)

    def test_synthetic_copy_and_input_preservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'source'
            output = Path(tmp)/'output'
            expected = {}
            for number, (folder, stem) in copy.TARGETS.items():
                for suffix in copy.FORMATS:
                    src = source/'04_figures'/folder/(stem+'.'+suffix)
                    src.parent.mkdir(parents=True, exist_ok=True)
                    payload = ('synthetic artifact '+str(number)+' '+suffix).encode()
                    src.write_bytes(payload)
                    expected[src] = payload
            with patch.object(sys, 'argv', ['copy', '--source-root', str(source),
                                            '--output-root', str(output)]):
                with contextlib.redirect_stdout(io.StringIO()):
                    copy.main()
            for src, original in expected.items():
                self.assertEqual(src.read_bytes(), original)
            receipt = json.loads((output/'audit/FIGURES3_5_BYTE_IDENTITY.json').read_text())
            self.assertEqual(len(receipt), len(copy.TARGETS)*len(copy.FORMATS))
            for row in receipt:
                self.assertEqual((source/row['input']).read_bytes(),
                                 (output/row['output']).read_bytes())
                self.assertTrue(row['byte_identical'])
                self.assertFalse(row['panel_content_changed'])

    def test_missing_input_preflight_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp)/'source', Path(tmp)/'output'
            with patch.object(sys, 'argv', ['copy', '--source-root', str(source),
                                            '--output-root', str(output)]):
                with self.assertRaises(FileNotFoundError):
                    copy.main()
            self.assertFalse(output.exists())

    def test_existing_output_symlink_rejected_when_supported(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'output'
            output.mkdir()
            original = Path(tmp)/'input'
            original.write_text('synthetic')
            try:
                (output/'link').symlink_to(original)
            except (OSError, NotImplementedError):
                self.skipTest('Symlink creation unavailable for this account.')
            with self.assertRaises(ValueError):
                reject_output_links(output)


if __name__ == '__main__':
    unittest.main()
