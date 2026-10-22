#!/usr/bin/env python3
"""Optional private fitted-model hash audit; projection needs membership only."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def audit_private_models(config, project_root, require=False):
    entries=config.get('private_model_artifacts',[])
    if require and not entries:raise ValueError('No explicitly configured private fitted-model artifacts')
    records=[]
    for item in entries:
        path=Path(item['path'])
        path=path if path.is_absolute() else Path(project_root)/path
        expected=item['sha256']
        if len(expected)!=64:raise ValueError('Expected a full SHA256, not an invented or partial identity')
        actual=digest(path)
        if actual!=expected:raise ValueError('Private fitted-model artifact SHA256 mismatch: '+item['role'])
        records.append({'role':item['role'],'path':str(path.resolve()),'bytes':path.stat().st_size,
            'sha256':actual,'status':'VERIFIED_PRIVATE_INPUT_NOT_REDISTRIBUTED'})
    return {'status':'PRIVATE_ARTIFACT_HASHES_VERIFIED' if records else 'NOT_REQUESTED_MEMBERSHIP_PROJECTION_DOES_NOT_REQUIRE_FITTED_WEIGHTS',
        'records':records,'public_model_weights_required':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--project-root',required=True,type=Path)
    parser.add_argument('--require',action='store_true')
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    if args.config.suffix.lower()=='.json':config=json.loads(args.config.read_text())
    else:
        import yaml
        config=yaml.safe_load(args.config.read_text())
    proof=audit_private_models(config,args.project_root,args.require)
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof,indent=2))


if __name__=='__main__':main()
