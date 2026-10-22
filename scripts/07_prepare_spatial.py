#!/usr/bin/env python3
"""Public raw/cached spatial entry point; frozen geometry, sparse 16-to-48 um.

Raw mode rebuilds in a new baseline cache. Cached mode validates a previously
prepared baseline and explicitly reports the skipped raw-processing stages.
The 08 entry point then writes a separate depth-results directory.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import urllib.request
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def load_config(config, project_root):
    config = Path(config).resolve()
    cfg = yaml.safe_load(config.read_text(encoding='utf8'))
    root = Path(project_root).resolve()
    for key, value in list(cfg.items()):
        if key.endswith('_config'):
            cfg[key] = (config.parent/value).resolve()
    for key in ['raw_tar', 'baseline_cache', 'cached_spatial', 'depth_output', 'figure_output']:
        cfg[key] = (root/cfg[key]).resolve()
    return cfg

def check_hash(path, expected):
    if not path.is_file() or sha(path) != expected:
        raise ValueError(f'Input absent or frozen SHA256 mismatch: {path.name}')

def checked_archive(cfg, allow_download):
    path = cfg['raw_tar']
    if not path.exists():
        if not allow_download:
            raise FileNotFoundError('Public GEO archive is absent; use --download-public or the public-data downloader')
        # The default URL is public GEO; credentials and signed URLs are unused.
        if not str(cfg['raw_url']).startswith('https://ftp.ncbi.nlm.nih.gov/'):
            raise ValueError('The built-in archive fetcher accepts public NCBI HTTPS only')
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name+'.partial')
        with urllib.request.urlopen(cfg['raw_url'], timeout=120) as response, partial.open('wb') as out:
            for block in iter(lambda: response.read(1024*1024), b''):
                out.write(block)
        if partial.stat().st_size != cfg['raw_bytes'] or sha(partial) != cfg['raw_sha256']:
            raise ValueError('Fetched archive size/hash differs from the inspected public release; partial file retained')
        partial.replace(path)
    check_hash(path, cfg['raw_sha256'])
    if path.stat().st_size != cfg['raw_bytes']:
        raise ValueError('Public GEO archive byte count mismatch')
    return path

def run(script, arguments, logs):
    command = [sys.executable, str(SCRIPT_DIR/script), *map(str, arguments)]
    print('Running '+script, flush=True)
    completed = subprocess.run(command, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    logs.parent.mkdir(parents=True, exist_ok=True)
    logs.write_text(completed.stdout, encoding='utf8')
    if completed.returncode:
        raise RuntimeError(f'{script} exited {completed.returncode}; inspect {logs.name}')
    return {'script': script, 'exit_code': completed.returncode, 'log_filename': logs.name}

def validate_cached(root):
    required = ['SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz', 'NP_PROGRAM_BY_REGION.tsv',
                'CANDIDATE_GENES_BY_REGION.tsv', 'PROGRAM_MEMBERS_USED_IN_SPATIAL_SCORE.tsv',
                'HUMAN_MOUSE_ORTHOLOG_COVERAGE.tsv', 'roi_validation/DISC_AF_NP_ROI_METADATA.tsv.gz']
    for name in required:
        if not (root/name).is_file():
            raise FileNotFoundError('Cached baseline is incomplete: '+name)
    captures = sorted(p.parent for p in (root/'prepared48').glob('*/counts48_symbols_qc.npz'))
    if len(captures) != 2:
        raise ValueError('Cached baseline must contain the two real sparse capture count matrices')
    for folder in captures:
        for name in ['BIN48_METADATA_QC.tsv.gz', 'MOUSE_SYMBOL_FEATURES.tsv']:
            if not (folder/name).is_file():
                raise FileNotFoundError('Missing count-aligned capture metadata: '+name)
    return {'captures': [p.name for p in captures],
            'source_table_sha256': sha(root/'SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz')}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True, type=Path)
    ap.add_argument('--project-root', required=True, type=Path)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--from-raw', action='store_true')
    mode.add_argument('--from-cached-inputs', action='store_true')
    ap.add_argument('--cached-spatial', type=Path, help='Read-only baseline override; never overwritten')
    ap.add_argument('--download-public', action='store_true')
    a = ap.parse_args()
    cfg = load_config(a.config, a.project_root)
    check_hash(cfg['programs_config'], cfg['programs_sha256'])
    check_hash(cfg['mapping_config'], cfg['mapping_sha256'])
    depth_cfg = json.loads(cfg['depth_config'].read_text())
    check_hash(cfg['roi_config'], depth_cfg['expected_roi_sha256'])
    manifest = {'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                'mode': 'from_raw' if a.from_raw else 'from_cached_inputs', 'stages': []}
    if a.from_cached_inputs:
        root = a.cached_spatial.resolve() if a.cached_spatial else cfg['cached_spatial']
        manifest.update(status='CACHED_BASELINE_VALIDATED_RAW_STAGES_SKIPPED',
                        skipped=['public_archive_fetch', 'sparse_16_to_48_preparation',
                                 'frozen_ROI_assignment', 'ortholog_coverage', 'fixed_program_scoring'])
    else:
        if a.cached_spatial:
            raise ValueError('--cached-spatial is only valid in cached mode')
        root = cfg['baseline_cache']
        if root == cfg['depth_output'] or root in cfg['depth_output'].parents or cfg['depth_output'] in root.parents:
            raise ValueError('Raw baseline cache and new depth output must be separate directories')
        archive = checked_archive(cfg, a.download_public)
        logdir = root/'pipeline_logs'
        tasks = [
            ('spatial_prepare.py', ['--raw-tar', archive, '--out', root/'prepared48',
                '--min-counts', cfg['min_counts'], '--min-genes', cfg['min_genes']]),
            ('spatial_assign_validate_rois.py', ['--prepared48', root/'prepared48',
                '--polygons', cfg['roi_config'], '--out', root/'roi_validation']),
            ('spatial_ortholog_coverage.py', ['--programs', cfg['programs_config'],
                '--targeted-mapping', cfg['mapping_config'], '--expected-mapping-sha256', cfg['mapping_sha256'],
                '--source-metadata', cfg['mapping_sources_config'], '--prepared48', root/'prepared48', '--out', root]),
            ('spatial_programs.py', ['--prepared48', root/'prepared48',
                '--roi-metadata', root/'roi_validation/DISC_AF_NP_ROI_METADATA.tsv.gz',
                '--coverage', root/'HUMAN_MOUSE_ORTHOLOG_COVERAGE.tsv', '--programs-file', cfg['programs_config'],
                '--programs', *cfg['programs'], '--candidates', *cfg['candidates'], '--out', root]),
        ]
        for index, (script, args) in enumerate(tasks):
            manifest['stages'].append(run(script, args, logdir/f'{index:02d}_{script}.log'))
        manifest.update(status='RAW_SPATIAL_BASELINE_COMPLETED', archive_sha256=sha(archive))
    manifest.update(validate_cached(root))
    manifest['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
    out = cfg['depth_output']
    out.mkdir(parents=True, exist_ok=True)
    (out/'SPATIAL_PIPELINE_ENTRY_STATUS.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf8')
    print(json.dumps(manifest, indent=2))

if __name__ == '__main__':
    main()
