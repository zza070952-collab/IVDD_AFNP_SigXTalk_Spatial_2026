"""Technical temporary-repository checks; no study observations or manuscripts."""
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1]/'scripts/authoring/audit_public_release.py'
SPEC = importlib.util.spec_from_file_location('publication_scope_gate', SCRIPT)
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)

@unittest.skipUnless(shutil.which('git'), 'Git is required for index-blob tests')
class PublicReleaseScopeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='ivdd-release-scope-')
        self.repo = Path(self.temp.name).resolve()
        self.assertTrue(self.repo.name.startswith('ivdd-release-scope-'))
        self.git('init', '-q')

    def tearDown(self):
        self.assertTrue(Path(self.temp.name).resolve() == self.repo)
        self.assertTrue(self.repo.name.startswith('ivdd-release-scope-'))
        self.temp.cleanup()

    def git(self, *args):
        return subprocess.check_output(['git', '-C', str(self.repo), '-c', 'core.autocrlf=false', *args],
                                       stderr=subprocess.STDOUT)

    def write(self, name, text, *, stage=True):
        path = self.repo/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding='utf-8')
        if stage:
            self.git('add', '--', name)
        return path

    def test_disallowed_document_and_data_suffixes(self):
        for name in ('draft.docx', 'source.xlsx', 'figure.pdf', 'matrix.npz'):
            self.write(name, 'not a real document or research fixture')
        result = gate.scan_repository(self.repo)
        self.assertEqual(result['files_scanned'], 4)
        self.assertEqual(result['status'], 'FAIL_PUBLICATION_SCOPE')
        self.assertTrue(all(any(x['kind']=='forbidden_data_document_binary_or_runtime_suffix'
                                for x in f['issues']) for f in result['files']))

    def test_private_payload_filename_is_rejected(self):
        self.write('build/EXPERIMENT_PAYLOAD.json', '{}')
        self.assertGreater(gate.scan_repository(self.repo)['issue_count'], 0)

    def test_untracked_analysis_and_smoke_are_not_scanned(self):
        self.write('scripts/clean.py', "field = 'Raw_Ct'\n")
        self.write('analysis/.smoke/source.xlsx', 'untracked private cache', stage=False)
        self.write('analysis/private.txt', 'ghp_'+'x'*36, stage=False)
        result = gate.scan_repository(self.repo)
        self.assertEqual(result['issue_count'], 0)
        self.assertEqual([f['path'] for f in result['files']], ['scripts/clean.py'])
        self.assertFalse(result['working_tree_and_untracked_scanned'])

    def test_code_fields_frozen_resources_and_method_configs_are_allowed(self):
        self.write('scripts/reader.py', "columns = ['Raw_Ct', 'Raw_OD450']\n")
        self.write('resources/frozen/programs.tsv', 'program\trank\tgene\nExampleProgram\t1\tEXAMPLE_GENE\n')
        self.write('config/method.example.json', '{"python":"/usr/bin/python3","fields":["Raw_Ct"]}')
        self.assertEqual(gate.scan_repository(self.repo)['issue_count'], 0)

    def test_real_format_secret_and_private_path_are_rejected(self):
        self.write('config/secret.txt', 'ghp_'+'z'*36)
        private_path = '/'+'home'+'/'+'technical_fixture_account'+'/private/input.tsv'
        self.write('config/local.txt', private_path)
        result = gate.scan_repository(self.repo)
        kinds = {x['kind'] for f in result['files'] for x in f['issues']}
        self.assertIn('credential_pattern', kinds)
        self.assertIn('private_absolute_path', kinds)

    def test_staged_blob_is_used_in_both_working_tree_directions(self):
        path = self.write('scripts/example.py', 'value = 1\n')
        path.write_text("token = '"+'ghp_'+'q'*36+"'\n", encoding='utf-8')
        self.assertEqual(gate.scan_repository(self.repo)['issue_count'], 0)
        self.git('add', 'scripts/example.py')
        path.write_text('value = 1\n', encoding='utf-8')
        self.assertGreater(gate.scan_repository(self.repo)['issue_count'], 0)

    def test_staged_only_after_clean_commit_and_HEAD_mode(self):
        self.write('scripts/base.py', 'base = True\n')
        self.git('-c', 'user.name=Technical fixture', '-c', 'user.email=fixture@example.invalid',
                 'commit', '-qm', 'Technical fixture baseline')
        self.write('scripts/change.py', 'change = True\n')
        self.assertEqual(gate.scan_repository(self.repo, staged_only=True)['files_scanned'], 1)
        head = gate.scan_repository(self.repo, source='HEAD')
        self.assertEqual([f['path'] for f in head['files']], ['scripts/base.py'])
        self.assertEqual(head['issue_count'], 0)

    def test_method_literal_review_note_is_not_failure(self):
        method = 'Portable method description. '*24
        self.write('scripts/method.py', 'DESCRIPTION = '+repr(method)+'\n')
        result = gate.scan_repository(self.repo)
        self.assertEqual(result['issue_count'], 0)
        self.assertEqual(result['manual_review_notes'], 1)

    def test_observation_payload_is_distinct_from_a_field_name(self):
        self.write('config/observations.json', '{"record":{"Raw_Ct":1.25}}')
        result = gate.scan_repository(self.repo)
        self.assertTrue(any(x['kind']=='experimental_observation_payload'
                            for f in result['files'] for x in f['issues']))

if __name__ == '__main__':
    unittest.main()
