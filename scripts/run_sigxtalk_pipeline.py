#!/usr/bin/env python3
"""Own orchestration of pinned official SigXTalk; never uses legacy predictions.

All matrix/model/complete-row outputs stay under the configured private project.
The public source contains wrappers and version metadata, not upstream source.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def save(data, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')


def execute(command, log, cwd, env):
    started = now()
    with log.open('w', encoding='utf-8') as fh:
        fh.write('# New execution; started UTC ' + started + '\n')
        fh.write('# argv=' + json.dumps([str(x) for x in command]) + '\n')
        fh.flush()
        result = subprocess.run([str(x) for x in command], cwd=cwd, env=env,
                                stdout=fh, stderr=subprocess.STDOUT)
    save({'argv': [str(x) for x in command], 'exitcode': result.returncode,
          'started_utc': started, 'finished_utc': now(),
          'log_sha256': digest(log)}, log.with_suffix(log.suffix + '.execution.json'))
    if result.returncode:
        raise RuntimeError(f'New stage exited {result.returncode}; inspect {log}')


def count_csv(path):
    with path.open(newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh)
        if not {'Receptor','TF','TG','pred_label'}.issubset(reader.fieldnames or []):
            raise ValueError('Official output headers are missing')
        n, active = 0, 0
        for row in reader:
            value = float(row['pred_label'])
            if not 0 <= value <= 1:
                raise ValueError('Nonfinite/out-of-range official PA prediction')
            n += 1
            active += value > .75
    return n, active


def verify_upstream(source, lock_file):
    lock = json.loads(lock_file.read_text())
    checked = []
    for rel, expected in lock['git_blob_sha1'].items():
        raw = (source / rel).read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
        if actual != expected:
            raise ValueError(f'Pinned upstream file differs: {rel}')
        checked.append({'file': rel, 'git_blob_sha1': actual, 'sha256': hashlib.sha256(raw).hexdigest()})
    return {'official_repository': lock['repository'], 'commit': lock['commit'],
            'source_file_identity_check': 'PASS', 'files': checked}


def patch_diff(wd):
    original = wd / 'audit/PRS_calc_upstream_function.R'
    runtime = wd / 'audit/PRS_calc_runtime_resource_function.R'
    if original.exists() and runtime.exists():
        text = ''.join(difflib.unified_diff(original.read_text().splitlines(True),
                                           runtime.read_text().splitlines(True),
                                           fromfile='SigXTalkR::PRS_calc (pinned upstream)',
                                           tofile='runtime worker/RNG resource wrapper'))
        (wd / 'audit/PRS_resource_rng.patch').write_text(text, encoding='utf-8')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--mode', choices=['smoke','formal','all'], default='all')
    ap.add_argument('--stage', choices=['prepare','hgnn','prs','all'], default='all')
    args = ap.parse_args()
    config_path = args.config.resolve()
    cfg = json.loads(config_path.read_text())
    root = Path(cfg['project_root'])
    if not root.is_absolute():
        root = (config_path.parent / root).resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)
    resolve = lambda name: Path(cfg[name]) if Path(cfg[name]).is_absolute() else root / cfg[name]
    module = resolve('module_dir')
    source = resolve('official_source')
    own = Path(__file__).resolve().parent
    lock = verify_upstream(source, own / 'sigxtalk_upstream_lock.json')
    module.mkdir(parents=True, exist_ok=True)
    save(lock, module / 'PINNED_UPSTREAM_FILE_PROOF.json')
    if cfg['hgnn']['epochs'] != 100 or cfg['seed'] != 20261017:
        raise ValueError('This preregistered formal design requires 100 epochs and seed20261017')
    if (not isinstance(cfg['smoke']['epochs'], int) or
            not 1 <= cfg['smoke']['epochs'] <= 100):
        raise ValueError('Technical smoke epochs must be a predeclared integer from 1 to 100')
    smoke_name = cfg['smoke'].get('directory_name', 'technical_smoke')
    if not isinstance(smoke_name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', smoke_name) or smoke_name == 'formal_new':
        raise ValueError('smoke.directory_name must be one safe directory name distinct from formal_new')
    smoke_dir = module / smoke_name
    env = dict(os.environ)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['PYTHONHASHSEED'] = str(cfg['seed'])
    env['PYTHONPATH'] = str(source / 'pythoncodes') + os.pathsep + env.get('PYTHONPATH','')
    env['R_LIBS_USER'] = os.pathsep.join(str(root / x) if not Path(x).is_absolute() else x
                                       for x in cfg.get('r_library_paths',[]))
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
        env[key] = str(cfg.get('threads',2))
    if cfg['hgnn']['device'] == 'cpu':
        env['CUDA_VISIBLE_DEVICES'] = ''
    modes = ['smoke','formal'] if args.mode == 'all' else [args.mode]
    for mode in modes:
        wd = module / 'formal_new' if mode == 'formal' else smoke_dir
        for sub in ['audit','logs','hgnn_outputs','tables']:
            (wd / sub).mkdir(parents=True, exist_ok=True)
        stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S')
        prepare_status = wd / 'audit/R_PREPARE_STATUS.json'
        if args.stage in ['prepare','all']:
            execute([cfg['rscript'], own / '03_prepare_sigxtalk.R', config_path, mode],
                    wd / 'logs' / f'{stamp}_R_prepare.log', root, env)
        if args.stage in ['hgnn','all']:
            if mode == 'formal':
                completed = smoke_dir / 'audit/NEW_R_PYTHON_R_COMPLETION.json'
                if not completed.is_file() or not json.loads(completed.read_text()).get('complete_R_python_R'):
                    raise RuntimeError('Real-data R→HGNN→PRS smoke closure is required before formal training')
                smoke_input = json.loads((smoke_dir / 'audit/NEW_INPUT_PROVENANCE.json').read_text())
                if smoke_input['config_sha256'] != digest(config_path):
                    raise ValueError('Configuration changed after smoke; rerun bounded real-data smoke for this frozen config')
            if not prepare_status.exists() or json.loads(prepare_status.read_text())['status'] != 'DONE':
                raise RuntimeError('New raw-count R preparation did not complete')
            provenance = json.loads((wd / 'audit/NEW_INPUT_PROVENANCE.json').read_text())
            if provenance['config_sha256'] != digest(config_path):
                raise ValueError('Configuration changed after input preparation')
            runtime = wd / 'main.py'
            shutil.copy2(source / 'vignettes/main.py', runtime)
            if digest(runtime) != digest(source / 'vignettes/main.py'):
                raise ValueError('Unchanged official main.py copy identity check failed')
            # Unchanged main.py writes ./outputs relative to its own location.
            # Redirect only the filesystem directory, not the official program or math.
            outputs = wd / 'outputs'
            if not outputs.exists():
                outputs.symlink_to((wd / 'hgnn_outputs').resolve(), target_is_directory=True)
            if outputs.resolve() != (wd / 'hgnn_outputs').resolve():
                raise ValueError('Official output path must resolve to this new run hgnn_outputs')
            project = provenance['project_id']
            pred = wd / 'hgnn_outputs' / project / f'pathways_{provenance["target_type"]}.csv'
            if pred.exists():
                raise FileExistsError('A new HGNN output already exists; inspect/resume PRS rather than silently retrain or overwrite')
            epochs = cfg['hgnn']['epochs'] if mode == 'formal' else cfg['smoke']['epochs']
            log = wd / 'logs' / f'{stamp}_official_HGNN.log'
            command = [cfg['python'], '-u', runtime, '--project', project,
                       '--target_type', provenance['target_type'], '--device', cfg['hgnn']['device'],
                       '--seed', str(cfg['seed']), '--epochs', str(epochs), '--save_model', 'True']
            execute(command, log, root, env)
            text = log.read_text(errors='replace')
            if len(re.findall(r'Starting \[ Epoch \d+ of ' + str(epochs) + r' \]',text)) != epochs or 'Training finished !!' not in text:
                raise RuntimeError('Official training termination evidence is incomplete')
            n, active = count_csv(pred)
            if not n:
                raise RuntimeError('Official prediction file is empty')
            save({'status':'DONE','new_computation':True,'mode':mode,'epochs_completed':epochs,
                  'official_main_sha256':digest(runtime),'prediction_file':str(pred),
                  'prediction_sha256':digest(pred),'all_tested_rows':n,'PA_gt075_rows':active,
                  'seed':cfg['seed'],'PYTHONHASHSEED':env['PYTHONHASHSEED'],
                  'input_provenance_sha256':digest(wd / 'audit/NEW_INPUT_PROVENANCE.json'),
                  'log_sha256':digest(log),'finished_utc':now()},wd / 'audit/HGNN_STATUS.json')
        if args.stage in ['prs','all']:
            try:
                execute([cfg['rscript'], own / '05_sigxtalk_prs.R', config_path, mode],
                        wd / 'logs' / f'{stamp}_R_PRS.log', root, env)
            finally:
                patch_diff(wd)
            completed = json.loads((wd / 'audit/NEW_R_PYTHON_R_COMPLETION.json').read_text())
            if not completed['complete_R_python_R']:
                raise RuntimeError('Real new R/HGNN/PRS closure was not obtained')
            completed['formal_design'] = mode == 'formal'
            completed['run_class'] = ('FULL_NEW' if completed['fidelity_specificity_status']=='DONE' else
                                      'NEW_CORE_COMPLETE_PARTIAL_WITHOUT_FID_SPEC') if mode=='formal' else (
                                      'REAL_DATA_TECHNICAL_SMOKE_COMPLETE_NOT_FORMAL')
            completed['new_result_not_legacy_relabel'] = True
            save(completed, wd / 'NEW_OFFICIAL_FINAL_STAGE_STATUS.json')
            print(json.dumps(completed, indent=2), flush=True)


if __name__ == '__main__':
    main()
