#!/usr/bin/env python3
"""Read-only publication scope gate over Git index blobs (or committed HEAD).

Untracked files and working-tree modifications are never inputs to this gate.
Long Python literals receive review notes, not automatic rejection of methods.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess

FORBIDDEN_SUFFIXES = {
    '.xlsx', '.xls', '.docx', '.doc', '.pdf', '.png', '.jpg', '.jpeg',
    '.tif', '.tiff', '.gif', '.zip', '.gz', '.tar', '.tgz', '.7z',
    '.npz', '.npy', '.mtx', '.parquet', '.loom', '.h5', '.h5ad', '.hdf5',
    '.rds', '.rdata', '.pt', '.pth', '.pkl', '.pickle', '.bundle', '.pem', '.key',
    '.pyc', '.exe', '.dll', '.so',
}
PRIVATE_DIRECTORIES = {
    'manuscript', 'manuscripts', 'revised_manuscript', 'current_user_inputs',
    'internal_review', 'source_unresolved', 'data_cache', 'node_modules', '__pycache__',
}
PRIVATE_PAYLOAD = re.compile(
    r'^(?:raw_extract|experiment_payload|compact_workbook_payload|workbook_payload|'
    r'manuscript_content|supplied_methods_transcription|supplied_sequences_exact)\.'
    r'|^(?:raw_qpcr_ct|raw_cck8_od450|sirna_qpcr_raw|program_qpcr_raw|cck8_raw)\.', re.I)
PRIVATE_PATHS = re.compile(
    r'(?:[A-Za-z]:[/\\]Users[/\\]|/Users/|/home/|/mnt/[A-Za-z]/Users/)'
    r'[A-Za-z0-9][A-Za-z0-9_. -]*[/\\]'
    r'|/root/(?:\.ssh|\.aws|\.config|Desktop|Documents)(?:/|\b)', re.I)
CREDENTIALS = [
    re.compile(r'gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}'),
    re.compile(r'\bAKIA[0-9A-Z]{16}\b'),
    re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    re.compile(r'https?://[^\s/@:\"\x27]+:[^\s/@\"\x27]+@', re.I),
    re.compile(r'Authorization\s*:\s*Bearer\s+[A-Za-z0-9._-]{20,}', re.I),
]
SECRET_ASSIGNMENT = re.compile(
    r'''\b(?:password|passwd|api[_-]?key|access[_-]?token|client[_-]?secret)\b'''
    r'''\s*[:=]\s*["']([^"'\r\n]+)["']''', re.I)
OBSERVATION_FIELDS = {'raw_ct', 'raw_od450'}

def _git(repo, *args, input_bytes=None):
    result = subprocess.run(['git', '-C', str(repo), *args], input=input_bytes,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError('Git could not read the requested repository object')
    return result.stdout

def _entries(repo, staged_only=False, source='index'):
    if source == 'HEAD':
        if staged_only:
            raise ValueError('--staged-only requires index source')
        records = _git(repo, 'ls-tree', '-r', '-z', 'HEAD').split(b'\0')
        entries = []
        for record in filter(None, records):
            meta, name = record.split(b'\t', 1)
            mode, kind, oid = meta.decode('ascii').split()
            entries.append((name.decode('utf-8', 'surrogateescape'), mode, oid, kind))
        return entries
    selected = None
    if staged_only:
        selected = set(_git(repo, 'diff', '--cached', '--name-only', '--diff-filter=ACMRT', '-z').split(b'\0'))
    entries = []
    for record in filter(None, _git(repo, 'ls-files', '--stage', '-z').split(b'\0')):
        meta, name = record.split(b'\t', 1)
        mode, oid, stage = meta.decode('ascii').split()
        if selected is not None and name not in selected:
            continue
        if stage != '0':
            raise ValueError('Unmerged index entries cannot be cleared for publication')
        entries.append((name.decode('utf-8', 'surrogateescape'), mode, oid,
                        'commit' if mode == '160000' else 'blob'))
    return entries

def _blobs(repo, entries):
    """Read object IDs from the index/tree through one binary-safe cat-file call."""
    wanted = [e for e in entries if e[3] == 'blob']
    if not wanted:
        return {}
    raw = _git(repo, 'cat-file', '--batch', input_bytes=''.join(e[2]+'\n' for e in wanted).encode('ascii'))
    result, offset = {}, 0
    for _, _, oid, _ in wanted:
        end = raw.index(b'\n', offset)
        header = raw[offset:end].decode('ascii').split()
        if len(header) != 3 or header[0] != oid or header[1] != 'blob':
            raise RuntimeError('Git returned an unexpected object type')
        size = int(header[2]);start = end+1
        result[oid] = raw[start:start+size]
        offset = start+size+1
    return result

def _placeholder(value):
    v = value.strip()
    return (not v or v.startswith(('$', '<')) or v.upper() in
            {'REDACTED', 'PLACEHOLDER', 'EXAMPLE', 'CHANGEME', 'YOUR_TOKEN', 'YOUR_API_KEY'})

def _json_observation(value):
    if isinstance(value, dict):
        if any(str(k).lower() in OBSERVATION_FIELDS and isinstance(v, (float, int))
               and not isinstance(v, bool) for k, v in value.items()):
            return True
        return any(_json_observation(v) for v in value.values())
    return isinstance(value, list) and any(_json_observation(v) for v in value)

def _numeric_text(value):
    try:
        float(value)
        return bool(value.strip())
    except (TypeError, ValueError):
        return False

def scan_blob(name, content):
    path = PurePosixPath(name)
    issues, notes = [], []
    if path.suffix.lower() in FORBIDDEN_SUFFIXES:
        issues.append({'kind': 'forbidden_data_document_binary_or_runtime_suffix'})
    if PRIVATE_DIRECTORIES.intersection(p.lower() for p in path.parts):
        issues.append({'kind': 'private_payload_directory'})
    if PRIVATE_PAYLOAD.search(path.name):
        issues.append({'kind': 'private_payload_filename'})
    try:
        text = content.decode('utf-8-sig')
        if '\0' in text:
            raise UnicodeError('binary')
    except UnicodeError:
        issues.append({'kind': 'nontext_blob'})
        return issues, notes
    for line_number, line in enumerate(text.splitlines(), 1):
        if PRIVATE_PATHS.search(line):
            issues.append({'kind': 'private_absolute_path', 'line': line_number})
        if any(pattern.search(line) for pattern in CREDENTIALS):
            issues.append({'kind': 'credential_pattern', 'line': line_number})
        if any(not _placeholder(m.group(1)) for m in SECRET_ASSIGNMENT.finditer(line)):
            issues.append({'kind': 'hardcoded_secret_assignment', 'line': line_number})
    if path.suffix.lower() == '.json':
        try:
            if _json_observation(json.loads(text)):
                issues.append({'kind': 'experimental_observation_payload'})
        except json.JSONDecodeError:
            issues.append({'kind': 'invalid_JSON'})
    elif path.suffix.lower() in {'.tsv', '.csv'}:
        reader = csv.DictReader(io.StringIO(text), delimiter='\t' if path.suffix.lower()=='.tsv' else ',')
        observation_columns = [x for x in reader.fieldnames or [] if x.lower() in OBSERVATION_FIELDS]
        if observation_columns and any(any(_numeric_text(row.get(k)) for k in observation_columns) for row in reader):
            issues.append({'kind': 'experimental_observation_table'})
    elif path.suffix.lower() == '.py':
        try:
            tree = ast.parse(text)
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant) and isinstance(node.value, str) and len(node.value)>500:
                    notes.append({'kind': 'long_literal_manual_scope_review', 'line': node.lineno,
                                  'characters': len(node.value)})
        except SyntaxError as error:
            issues.append({'kind': 'python_syntax', 'line': error.lineno})
    return issues, notes

def scan_repository(repo, *, staged_only=False, source='index'):
    entries = _entries(Path(repo), staged_only=staged_only, source=source)
    contents = _blobs(Path(repo), entries)
    files = []
    for name, mode, oid, kind in entries:
        if kind != 'blob':
            issues, notes, content = [{'kind': 'non_blob_git_entry'}], [], b''
        else:
            content = contents[oid]
            issues, notes = scan_blob(name, content)
        files.append({'path': name, 'mode': mode, 'git_object': oid, 'bytes': len(content),
                      'sha256': hashlib.sha256(content).hexdigest(), 'issues': issues,
                      'review_notes': notes})
    issue_count = sum(len(f['issues']) for f in files)
    return {'status': 'FAIL_PUBLICATION_SCOPE' if issue_count else 'PASS_TRACKED_BLOB_SCOPE',
            'source': 'staged_index' if staged_only else source, 'files_scanned': len(files),
            'issue_count': issue_count, 'manual_review_notes': sum(len(f['review_notes']) for f in files),
            'working_tree_and_untracked_scanned': False, 'files': files,
            'limits': 'Pattern screening does not replace review of manuscript prose, permissions or reproducibility claims.'}

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', '--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--staged-only', action='store_true')
    parser.add_argument('--source', choices=['index', 'HEAD'], default='index')
    args = parser.parse_args(argv)
    report = scan_repository(args.repo_root, staged_only=args.staged_only, source=args.source)
    print(json.dumps(report, indent=2, ensure_ascii=True))
    return 1 if report['issue_count'] else 0

if __name__ == '__main__':
    raise SystemExit(main())
