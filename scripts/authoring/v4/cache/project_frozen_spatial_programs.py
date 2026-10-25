#!/usr/bin/env python3
"""Expose the existing strict sparse spatial projector for all five fixed programs.

This entry point writes a new output cache only. It does not prepare counts,
assign/move ROIs, fit programs, train a model, or repeat a depth-sampling draw.
The scoring implementation remains scripts/spatial_programs.py.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import subprocess
import sys

SCRIPTS = Path(__file__).resolve().parents[3]
REPO = SCRIPTS.parent
CANONICAL_PROGRAMS = REPO / 'resources/frozen_k5/FROZEN_K5_PROGRAMS.tsv'
IMPLEMENTATION = SCRIPTS / 'spatial_programs.py'
PROGRAMS = ['GEP1', 'GEP2', 'GEP3', 'GEP4', 'GEP5']
CANDIDATES = 'FN1 ITGB1 PRSS23 CCN2 PLOD2 THBS1 LOX TGFBI GJA1 CD44 IGFBP3'.split()

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def signatures(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not {'program', 'rank', 'gene'}.issubset(reader.fieldnames or []):
            raise ValueError('Frozen program table requires program/rank/gene')
        rows = list(reader)
    result = {}
    for program in PROGRAMS:
        selected = sorted((r for r in rows if r['program'] == program), key=lambda r: int(r['rank']))
        genes = [r['gene'] for r in selected]
        if [int(r['rank']) for r in selected] != list(range(1, 26)) or len(set(genes)) != 25:
            raise ValueError('Expected 25 unique ranked members in every frozen program')
        result[program] = genes
    if len(rows) != 125:
        raise ValueError('Expected exactly the five fixed 25-member programs')
    return result

def resolve_inputs(args):
    base = args.input_work / '03_spatial_GSE287499' if args.input_work else None
    prepared = args.prepared48 or (base / 'prepared48' if base else None)
    roi = args.roi_metadata or (base / 'roi_validation/DISC_AF_NP_ROI_METADATA.tsv.gz' if base else None)
    coverage = args.coverage or (base / 'HUMAN_MOUSE_ORTHOLOG_COVERAGE.tsv' if base else None)
    if any(p is None for p in (prepared, roi, coverage)):
        raise ValueError('Supply --input-work, or explicit prepared48, ROI metadata and coverage files')
    return prepared.resolve(), roi.resolve(), coverage.resolve(), args.programs_file.resolve()

def protect_output(output, prepared, input_files):
    output = output.resolve()
    if output == prepared or output.is_relative_to(prepared):
        raise ValueError('Output must be outside the prepared-input directory')
    for path in input_files:
        if path == output or path.is_relative_to(output):
            raise ValueError('Output must not contain an original input')
        if output == path.parent or output.is_relative_to(path.parent):
            raise ValueError('Output must be separate from the original input file directory')
    if output.exists() and not output.is_dir():
        raise ValueError('Output directory points to an existing file')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be new or empty; completed input caches are not overwritten')
    return output

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-work', '--work', type=Path, help='Optional old-layout input root, read-only')
    parser.add_argument('--prepared48', type=Path)
    parser.add_argument('--roi-metadata', type=Path)
    parser.add_argument('--coverage', '--mapping-file', dest='coverage', type=Path,
                        help='Per-capture HUMAN_MOUSE_ORTHOLOG_COVERAGE.tsv, not the global mapping-only table')
    parser.add_argument('--programs-file', type=Path, default=CANONICAL_PROGRAMS)
    parser.add_argument('--programs', nargs='+', default=PROGRAMS)
    parser.add_argument('--candidates', nargs='+', default=CANDIDATES)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if len(set(args.programs)) != len(args.programs) or not set(args.programs).issubset(PROGRAMS):
        parser.error('Programs must be a unique subset of GEP1 through GEP5')
    prepared, roi, coverage, programs = resolve_inputs(args)
    for path in (roi, coverage, programs, CANONICAL_PROGRAMS, IMPLEMENTATION):
        if not path.is_file():
            raise ValueError('A required explicit input or public source file is missing')
    if signatures(programs) != signatures(CANONICAL_PROGRAMS):
        raise ValueError('Input membership differs from the distributed fixed 125-member definition')
    paths = {'roi_metadata': roi, 'per_capture_ortholog_coverage': coverage,
             'frozen_programs': programs, 'distributed_frozen_programs': CANONICAL_PROGRAMS,
             'scoring_implementation': IMPLEMENTATION,
             'projection_entry_point': Path(__file__).resolve()}
    folders = sorted(p.parent for p in prepared.glob('*/counts48_symbols_qc.npz'))
    if not folders:
        raise ValueError('No prepared sparse capture matrices were found')
    for i, folder in enumerate(folders, 1):
        for name in ('counts48_symbols_qc.npz', 'BIN48_METADATA_QC.tsv.gz', 'MOUSE_SYMBOL_FEATURES.tsv'):
            path = folder / name
            if not path.is_file():
                raise ValueError('A prepared capture is missing counts, bin metadata or feature symbols')
            paths[f'capture_{i}/{name}'] = path
    if (prepared / 'SPATIAL_INPUT_QC.tsv').exists():
        paths['prepared_input_QC'] = prepared / 'SPATIAL_INPUT_QC.tsv'
    output = protect_output(args.output_dir, prepared, list(paths.values()))
    before = {key: sha(path) for key, path in paths.items()}
    output.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(IMPLEMENTATION), '--prepared48', str(prepared),
               '--roi-metadata', str(roi), '--coverage', str(coverage),
               '--programs-file', str(programs), '--programs', *args.programs,
               '--candidates', *args.candidates, '--out', str(output)]
    # The existing public implementation computes unweighted mean log1p(CP10k),
    # with measured zeros and fixed available strict-one-to-one members retained.
    subprocess.run(command, check=True)
    after = {key: sha(path) for key, path in paths.items()}
    if before != after:
        raise RuntimeError('A protected input changed during projection')
    private_pattern = re.compile(r'[A-Za-z]:[/\\]Users[/\\][^/\\\s]+[/\\]|/home/[^/\s]+/|/mnt/(?:data)/')
    for path in output.iterdir():
        if path.suffix in ('.tsv', '.json') and private_pattern.search(path.read_text(encoding='utf-8')):
            raise RuntimeError('Output contains a machine-specific path; review input metadata before sharing')
    products = {p.name: {'bytes': p.stat().st_size, 'sha256': sha(p)}
                for p in sorted(output.iterdir()) if p.is_file()}
    versions = {}
    for package in ('numpy', 'pandas', 'scipy'):
        try: versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: versions[package] = 'unavailable'
    receipt = {'status': 'DONE_FROZEN_SPATIAL_PROGRAM_PROJECTION',
               'finished_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
               'implementation': 'scripts/spatial_programs.py', 'programs': args.programs,
               'protected_input_sha256': before, 'protected_inputs_unchanged': True,
               'products': products, 'versions': versions,
               'counts_prepared_here': False, 'ROI_assignment_or_geometry_changed': False,
               'program_fitting_or_model_training': False, 'depth_sampling_performed': False,
               'genotype_inference_performed': False, 'absolute_machine_paths_recorded': False}
    (output / 'FROZEN_SPATIAL_PROJECTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    main()
