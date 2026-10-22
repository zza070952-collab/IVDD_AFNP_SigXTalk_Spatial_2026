#!/usr/bin/env python3
"""Audit and promote completed current official SigXTalk outputs; never synthesize success.

All generated objects, barcode tables, logs and full predictions are PRIVATE runtime
artifacts. Publish this own script, not the generated project directory.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np
import pandas as pd

CANDIDATES = ['FN1','ITGB1','PLOD2','CCN2','PRSS23','THBS1','LOX','TGFBI']
ORIGIN = 'NEW_OFFICIAL_SIGXTALK_RAW_COUNTS_R_HGNN_PRS'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024**2),b''):h.update(chunk)
    return h.hexdigest()


def load(path):return json.loads(Path(path).read_text())
def read(path):return pd.read_csv(path,sep='\t',float_precision='round_trip')
def write(data,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    data.to_csv(path,sep='\t',index=False,na_rep='NA',compression='infer')
def save(data,path):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def truth(x):return str(x).lower() in ['true','1']
def require(condition,message):
    if not condition:raise ValueError(message)
def receipt(path):
    p=Path(path)
    return {'path':str(p.resolve()),'bytes':p.stat().st_size,'sha256':sha(p)}
def table_receipt(path):
    out=receipt(path);p=Path(path)
    if p.name.endswith(('.tsv','.tsv.gz','.csv','.csv.gz')):
        opener=gzip.open if p.suffix=='.gz' else open
        with opener(p,'rt',encoding='utf-8',newline='') as f:
            r=csv.reader(f,delimiter=',' if '.csv' in p.name else '\t')
            out['columns']=next(r,[]);out['data_rows']=sum(1 for _ in r)
    return out
def finite(value):return float(value) if pd.notna(value) and np.isfinite(value) else None


def verify_execution(wd,kind,expected_log_sha=None):
    found=[]
    for p in sorted((wd/'logs').glob('*.execution.json')):
        r=load(p);cmd=[str(x) for x in r.get('argv',[])]
        match=(kind=='hgnn' and any(x.endswith('main.py') for x in cmd) or
               kind=='prepare' and any(x.endswith('03_prepare_sigxtalk.R') for x in cmd) or
               kind=='prs' and any(x.endswith('05_sigxtalk_prs.R') for x in cmd))
        if not match or r.get('exitcode')!=0:continue
        log=Path(str(p)[:-len('.execution.json')])
        require(log.is_file(),'Successful execution receipt has no original log')
        require(sha(log)==r['log_sha256'],'Execution log hash mismatch')
        if expected_log_sha and r['log_sha256']!=expected_log_sha:continue
        found.append(dict(r,execution_receipt=receipt(p),log=receipt(log)))
    require(bool(found),'No current successful '+kind+' subprocess receipt')
    return found[-1]


def raw_candidate_detection(input_dir,meta,feature_map,candidates,sender,receiver):
    """Extract only declared gene counts from the actual integer coordinate matrix."""
    genes=pd.read_csv(input_dir/'genes.tsv',sep='\t',header=None,dtype=str)[0].tolist()
    target_rows={i+1:g for i,g in enumerate(genes) if g in candidates}
    group=meta.set_index('cell_id').bridge_group
    groups={sender:set(np.flatnonzero(meta.bridge_group.to_numpy()==sender)+1),
            receiver:set(np.flatnonzero(meta.bridge_group.to_numpy()==receiver)+1)}
    detected={(g,k):set() for g in candidates for k in groups};umi={(g,k):0 for g in candidates for k in groups}
    with gzip.open(input_dir/'counts.mtx.gz','rt') as f:
        header=f.readline().strip();require('coordinate' in header,'Expected coordinate raw count MatrixMarket')
        for line in f:
            if line.startswith('%'):continue
            shape=[int(x) for x in line.split()];break
        require(shape[:2]==[len(genes),len(meta)],'Raw matrix/meta dimensions disagree')
        for line in f:
            r,c,v=line.split();r=int(r)
            if r not in target_rows:continue
            c=int(c);v=float(v);g=target_rows[r]
            require(v>=0 and v.is_integer(),'Candidate count is not a nonnegative integer')
            if v>0:
                for k,ids in groups.items():
                    if c in ids:detected[g,k].add(c);umi[g,k]+=int(v)
    result={}
    for g in candidates:
        result[g]={'in_original_raw_features':g in genes}
        for k,label in [(sender,'sender'),(receiver,'receiver')]:
            n=len(groups[k]);result[g][label+'_raw_counts_detection']=len(detected[g,k])/n if g in genes and n else None
            result[g][label+'_raw_UMI_sum']=umi[g,k] if g in genes else None
    return result


def legacy_comparison(legacy_root,new_counts,new_candidates):
    rows=[];manifest=[]
    def add(metric,new,old=None,source=None,note=''):
        rows.append(dict(metric=metric,new_value=new,legacy_value=old,
            legacy_status='NOT_REQUESTED' if legacy_root is None else ('AVAILABLE' if old is not None else 'NOT_AVAILABLE'),
            legacy_source=source,legacy_method_origin='LEGACY_SIGXTALK_REANALYSIS_COMPARISON_ONLY',
            new_method_origin=ORIGIN,interpretation=note))
    old={};old_candidates={};source=None
    if legacy_root:
        lr=Path(legacy_root).resolve();require(lr.is_dir(),'Requested private legacy root does not exist')
        def find(name):
            hit=sorted(lr.rglob(name));return hit[0] if hit else None
        p=find('RTF_TARGET_ALL_TESTED.tsv.gz')
        if p:
            d=read(p);manifest.append(receipt(p));source=str(p)
            pa=next((c for c in ['pred_label','PA','original_PA'] if c in d),None)
            prs=next((c for c in ['PRS','Weight','original_PRS'] if c in d),None)
            tg=next((c for c in ['TG','Target'] if c in d),None)
            old['HGNN_tested_axes']=len(d)
            if pa:old['PA_gt075_axes']=int((pd.to_numeric(d[pa],errors='coerce')>.75).sum())
            if prs:
                z=pd.to_numeric(d[prs],errors='coerce').dropna();old['returned_positive_PRS_rows']=len(z)
                for key,func in [('PRS_min','min'),('PRS_median','median'),('PRS_max','max')]:old[key]=finite(getattr(z,func)())
            if tg:
                old['prepared_target_genes']=d[tg].nunique()
                for g in CANDIDATES:
                    q=d[d[tg]==g];old_candidates[g]={'tested':len(q),'active':int((pd.to_numeric(q[pa],errors='coerce')>.75).sum()) if pa else None,
                                                   'prs':int(pd.to_numeric(q[prs],errors='coerce').notna().sum()) if prs else None}
        p=find('SENDER_RECEIVER_INPUT.tsv')
        if p:
            d=read(p);manifest.append(receipt(p))
            nc=next((c for c in ['n_cells','Freq','cells'] if c in d),None)
            donor=next((c for c in ['donor_id','patient_id','donor'] if c in d),None)
            if nc:old['input_cells']=int(d[nc].sum())
            if donor:old['independent_patients']=d[donor].nunique()
        p=find('AF_TO_NP_CANDIDATE_LR.tsv')
        if p:
            d=read(p);manifest.append(receipt(p));old['AFNP_LR_rows']=len(d)
            if {'ligand','receptor'}.issubset(d):old['FN1_complexes']=';'.join(sorted(d.loc[d.ligand=='FN1','receptor'].dropna().unique()))
    for key,value in new_counts.items():
        add(key,value,old.get(key),source,'Differences are allowed; raw feature universe, installed versions and declared target rule can differ. No matching-to-old optimization.')
    for r in new_candidates.itertuples():
        for suffix,new in [('tested',r.n_tested_axes_as_target),('active',r.n_PA_gt075_as_target),('prs',r.n_returned_positive_PRS_as_target)]:
            add(r.gene+'_'+suffix,new,old_candidates.get(r.gene,{}).get(suffix),source,'Target testability only; an untested gene is not a negative biological effect.')
    return pd.DataFrame(rows),manifest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--legacy-root',type=Path,help='Optional PRIVATE legacy artifacts, comparison only; never copied to NEW tables')
    ap.add_argument('--environment-record',type=Path,action='append',default=[],help='Additional actual environment receipt file; repeatable')
    args=ap.parse_args();cp=args.config.resolve();cfg=load(cp)
    root=Path(cfg['project_root']);root=root if root.is_absolute() else (cp.parent/root).resolve()
    resolve=lambda k: Path(cfg[k]) if Path(cfg[k]).is_absolute() else root/cfg[k]
    module=resolve('module_dir');wd=module/'formal_new';td=wd/'tables';pd_=wd/'prs_outputs'
    status=load(wd/'NEW_OFFICIAL_FINAL_STAGE_STATUS.json');completion=load(wd/'audit/NEW_R_PYTHON_R_COMPLETION.json')
    inp=load(wd/'audit/NEW_INPUT_PROVENANCE.json');hg=load(wd/'audit/HGNN_STATUS.json')
    prep=load(wd/'audit/R_PREPARE_STATUS.json');prs_status=load(wd/'audit/R_PRS_STATUS.json')
    require(status.get('formal_design') and status.get('new_result_not_legacy_relabel') and status.get('complete_R_python_R'),'No actual current formal R/HGNN/PRS closure')
    require(completion.get('complete_R_python_R') and completion.get('no_old_pathways_or_PRS_used'),'PRS closure or new-only provenance missing')
    require(inp['config_sha256']==sha(cp),'Current config differs from formal preparation')
    require(inp['input_kind']=='real_GSE230808_raw_counts' and not inp['started_from_old_pathway_output'],'Formal inputs are not a new raw GSE230808 reconstruction')
    require(hg['status']=='DONE' and hg['epochs_completed']==100 and hg['seed']==20261017,'Formal HGNN100epochs/seed evidence missing')
    require(prep['status']=='DONE' and prs_status['status']=='DONE','R stage status not DONE')
    provfile=wd/'audit/NEW_INPUT_PROVENANCE.json'
    require(hg['input_provenance_sha256']==sha(provfile)==completion['new_input_provenance_sha256'],'Formal input provenance hash chain changed')
    pred=Path(hg['prediction_file']);pred=pred if pred.is_absolute() else root/pred
    require(sha(pred)==hg['prediction_sha256']==completion['new_HGNN_predictions_sha256'],'New prediction hash changed')
    rawfile=pd_/'PRS_RAW_OFFICIAL_NEW.tsv'
    require(sha(rawfile)==completion['new_raw_PRS_sha256'],'New raw PRS hash changed')
    executions={k:verify_execution(wd,k,hg['log_sha256'] if k=='hgnn' else None) for k in ['prepare','hgnn','prs']}
    text=Path(executions['hgnn']['log']['path']).read_text(errors='replace')
    epochs=re.findall(r'Starting \[ Epoch (\d+) of 100 \]',text)
    require(([int(x) for x in epochs] in [list(range(1,101)),list(range(100))]) and 'Training finished !!' in text,'Incomplete actual100epoch log')
    smoke_dir=module/cfg['smoke'].get('directory_name','technical_smoke')
    smoke=load(smoke_dir/'audit/NEW_R_PYTHON_R_COMPLETION.json');smoke_input=load(smoke_dir/'audit/NEW_INPUT_PROVENANCE.json')
    require(smoke.get('complete_R_python_R') and smoke_input['config_sha256']==sha(cp),'Real-data technical smoke closure is absent or a different config')
    require(smoke.get('positive_PRS_rows',0)>0 and smoke.get('no_old_pathways_or_PRS_used'),'Smoke lacks actual nonempty new PRS')
    smoke_hg=load(smoke_dir/'audit/HGNN_STATUS.json')
    require(smoke_hg['epochs_completed']==cfg['smoke']['epochs'],'Smoke epochs differ from actual completed configuration')
    smoke_exec={k:verify_execution(smoke_dir,k,smoke_hg['log_sha256'] if k=='hgnn' else None) for k in ['prepare','hgnn','prs']}
    sourceproof=load(module/'PINNED_UPSTREAM_FILE_PROOF.json')
    source=resolve('official_source')
    for r in sourceproof['files']:require(sha(source/r['file'])==r['sha256'],'Pinned official source changed')
    require(hg['official_main_sha256']==sha(source/'vignettes/main.py')==sha(wd/'main.py'),'Formal main.py is not unchanged pinned upstream')
    api=read(wd/'audit/R_API_AND_BODY_AUDIT.tsv')
    require(api.body_matches_pinned_checkout.map(truth).all(),'Installed R algorithm body differs from pinned source')
    full=read(td/'ALL_TESTED_NEW_AXES_WITH_PRS.tsv');pa=pd.read_csv(pred,float_precision='round_trip')
    keys=['Receptor','TF','TG'];require(not full.duplicated(keys).any(),'Duplicate new tested axes')
    require(len(full)==len(pa)==completion['all_tested_HGNN_rows']==hg['all_tested_rows'],'Full new tested row counts disagree')
    check=full.merge(pa[keys+['pred_label']],on=keys,validate='one_to_one',suffixes=('','_original'))
    require(np.allclose(check.pred_label.to_numpy(),check.pred_label_original.to_numpy(),rtol=0,atol=1e-14),'PA values differ beyond R text-export precision from original new predictions')
    # The canonical promoted table takes PA directly from the original Python CSV.
    canonical_pa=pa.set_index(keys).pred_label.reindex(pd.MultiIndex.from_frame(full[keys])).to_numpy()
    full['pred_label']=canonical_pa
    require(set(full.method_origin)=={ORIGIN} and not full.source_run_started_from_legacy_predictions.map(truth).any(),'Mixed historical source rows in NEW results')
    require(np.isfinite(full.pred_label).all() and full.pred_label.between(0,1).all(),'Invalid original new PA')
    raw=read(rawfile);filtered=read(pd_/'PRS_FILTERED_OFFICIAL_NEW.tsv')
    require(len(raw)==completion['positive_PRS_rows']>0 and len(filtered)==completion['filtered_PRS_rows']>0,'Actual PRS counts are empty/inconsistent')
    require((raw.Weight>0).all() and np.isfinite(raw.Weight).all(),'Returned PRS is not positive finite')
    joined=full.merge(raw.rename(columns={'SSC':'TF','Target':'TG','Weight':'PRS_original'}),on=keys,how='left',validate='one_to_one')
    require(joined.PRS.isna().equals(joined.PRS_original.isna()) and np.array_equal(joined.PRS.fillna(-1),joined.PRS_original.fillna(-1)),'Full PRS values/missingness differ from original new return')
    active=full[full.pred_label>cfg['hgnn']['activation_threshold']].copy()
    require(len(active)==completion['active_PA_rows']==hg['PA_gt075_rows'],'Active-axis counts disagree')
    meta=read(wd/'input/selected_cell_metadata.tsv');rawmeta=read(resolve('input_dir')/'meta.tsv')
    for key,p in [('raw_matrix_sha256',resolve('input_dir')/'counts.mtx.gz'),('raw_genes_sha256',resolve('input_dir')/'genes.tsv'),
                  ('cell_identity_sha256',resolve('input_dir')/'meta.tsv'),('patient_pairing_sha256',resolve('pairing_tsv'))]:
        require(inp['new_input_file_hashes'][key]==sha(p),'Actual raw input changed after preparation: '+key)
    require(not meta.cell_id.duplicated().any() and meta.cell_id.tolist()==rawmeta.cell_id.tolist(),'Selected raw metadata order/identity changed')
    require(len(meta)==inp['n_cells'] and meta.donor_id.nunique()==inp['independent_patient_count'],'Actual cell/donor counts disagree')
    barcode_audit=module/'NEW_INPUT_CELL_DONOR_BARCODE_AUDIT.tsv'
    require(barcode_audit.is_file(),'Missing actual raw barcode audit; run raw reconstruction audit first')
    ba=read(barcode_audit);require(set(ba.cell_id)==set(meta.cell_id) and len(ba)==len(meta),'Raw barcode audit does not cover exactly formal cells')
    if 'raw_barcode_found' in ba:require(ba.raw_barcode_found.map(truth).all(),'Some selected barcodes not found in raw data')
    if 'total_UMI_equal' in ba:require(ba.total_UMI_equal.map(truth).all(),'Raw/cached UMI audit mismatch')
    lr=read(td/'FOCUS_LR_COMPLEX_ALL.tsv');all_lr=read(td/'CELLCHAT_ALL_SOURCE_TARGET_LR.tsv');components=read(td/'LR_COMPLEX_COMPONENTS.tsv')
    require((lr.source==cfg['sender']).all() and (lr.target==cfg['receiver']).all(),'Focused new LR contains another sender/receiver')
    for d in [full,active,lr,all_lr,components,raw,filtered]:d['method_origin']=ORIGIN
    write(all_lr,module/'NEW_CELLCHAT_AFNP_ALL_LR.tsv.gz')
    write(lr,module/'NEW_AFNP_LR_COMPLEX.tsv');write(components,module/'NEW_RECEPTOR_COMPLEX_COMPONENTS.tsv')
    write(full,module/'NEW_HGNN_TESTED_AXES.tsv.gz');write(active,module/'NEW_HGNN_ACTIVATED_AXES.tsv.gz')
    prs_all=full.copy()
    prs_all['PRS_missing_reason']=np.where(prs_all.PRS.notna(),'returned_positive_official_PRS',np.where(prs_all.pred_label<=cfg['hgnn']['activation_threshold'],'PA_not_above_fixed_threshold_PRS_not_requested','official_PRS_not_returned_positive_importance_or_coexpression_model_constraint'))
    prs_all['NA_is_not_zero']=True
    write(prs_all,module/'NEW_PRS_ALL.tsv.gz');write(filtered,module/'NEW_PRS_FILTERED.tsv')
    write(raw,module/'NEW_PRS_RAW_RETURNED_POSITIVE.tsv')
    metrics=[]
    for label,name in [('all_new_returned_positive_PRS_rows','FID_SPE_RAW_UNIVERSE_OFFICIAL_NEW.tsv'),('new_Filter_results_relative_max_filtered_rows','FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv')]:
        p=pd_/name
        if completion['fidelity_specificity_status']=='DONE':
            require(p.is_file(),'Fid/spe marked DONE without actual table');d=read(p)
            require(len(d)>0 and np.isfinite(d[['fidelity_target_all','specificity_receptor_all']]).all().all(),'Invalid actual official fid/spe')
            d['normalization_universe']=label;d['method_origin']=ORIGIN;d['official_metric_status']='DONE';metrics.append(d)
    if metrics:write(pd.concat(metrics,ignore_index=True),module/'NEW_FIDELITY_SPECIFICITY_IF_AVAILABLE.tsv')
    else:write(pd.DataFrame([dict(official_metric_status='PRS_DONE_FID_SPEC_PARTIAL',error=completion.get('fidelity_specificity_error'),no_substitute_like_metric=True)]),module/'NEW_FIDELITY_SPECIFICITY_IF_AVAILABLE.tsv')
    paths=read(td/'GEP2_NEW_SOURCE_SUPPORTED_PATHS.tsv.gz');paths['method_origin']=ORIGIN
    write(paths,module/'NEW_GEP2_TARGET_PATHS.tsv')
    selection=read(td/'TARGET_GENE_SELECTION_ALL.tsv');gene_map=read(td/'SEURAT_FEATURE_IDENTITY_MAP.tsv')
    detections=raw_candidate_detection(resolve('input_dir'),rawmeta,gene_map,CANDIDATES,cfg['sender'],cfg['receiver'])
    candidate_rows=[]
    for gene in CANDIDATES:
        q=full[full.TG==gene];r=full[full.Receptor==gene];s=full[full.TF==gene]
        cc=components[components.receptor_component==gene];lg=lr[lr.ligand==gene]
        sl=selection[selection.gene==gene];selected=bool(len(sl) and sl.selected_target.map(truth).any())
        reason=('not_in_raw_feature_table' if not detections[gene]['in_original_raw_features'] else
                'not_selected_by_fixed_variance_rule' if not selected else 'removed_by_official_expression_or_prior_filter' if not len(q) else
                'tested_but_no_PA_gt075' if not (q.pred_label>.75).any() else 'activated_but_no_returned_positive_PRS' if not q.PRS.notna().any() else 'returned_positive_PRS')
        candidate_rows.append(dict(gene=gene,query_role='sender_ligand' if gene=='FN1' else 'receptor_complex_component' if gene=='ITGB1' else 'predeclared_target_or_readout_comparison',
            **detections[gene],selected_by_fixed_target_rule=selected,n_tested_axes_as_target=len(q),n_PA_gt075_as_target=int((q.pred_label>.75).sum()),
            n_returned_positive_PRS_as_target=int(q.PRS.notna().sum()),n_filtered_PRS_as_target=int(q.PRS_relative_max_filter_pass.map(truth).sum()),
            max_original_PA_as_target=finite(q.pred_label.max()),max_original_PRS_as_target=finite(q.PRS.max()),target_status=reason,
            n_tested_axes_as_receptor=len(r),n_tested_axes_as_database_SSC=len(s),n_AFNP_LR_as_ligand=len(lg),
            n_AFNP_LR_as_receptor_component=len(cc),original_AFNP_receptor_complexes=';'.join(sorted(cc.receptor_complex.dropna().unique())),
            n_GEP2_source_supported_paths_as_target=int((paths.TG==gene).sum()) if 'TG' in paths else 0,
            gene_not_forced_into_targets=True,missing_PRS_not_zero=True,method_origin=ORIGIN))
    candidates=pd.DataFrame(candidate_rows);write(candidates,module/'NEW_CANDIDATE_RECEPTOR_AXIS.tsv')
    counts={'input_cells':len(meta),'independent_patients':int(meta.donor_id.nunique()),'raw_features':inp['n_genes'],
        'selected_targets':inp['n_selected_targets'],'prepared_target_genes':inp['prepared_target_genes'],'HGNN_tested_axes':len(full),
        'PA_gt075_axes':len(active),'returned_positive_PRS_rows':len(raw),'filtered_PRS_rows':len(filtered),'AFNP_LR_rows':len(lr),
        'FN1_complexes':';'.join(sorted(lr.loc[lr.ligand=='FN1','receptor'].dropna().unique())),
        'PRS_min':finite(raw.Weight.min()),'PRS_median':finite(raw.Weight.median()),'PRS_max':finite(raw.Weight.max())}
    comparison,legacy_manifest=legacy_comparison(args.legacy_root,counts,candidates)
    write(comparison,module/'NEW_VS_LEGACY_AUDIT.tsv')
    direction=[]
    for source_group in sorted(meta.bridge_group.unique()):
        for target_group in sorted(meta.bridge_group.unique()):
            q=all_lr[(all_lr.source==source_group)&(all_lr.target==target_group)]
            direction.append(dict(source=source_group,target=target_group,exported_LR_rows=len(q),p005_rows=int((q.pval<.05).sum()),
                complete_zero_entries='recoverable_only_from_saved_full_CellChat_object',method_origin=ORIGIN))
    write(pd.DataFrame(direction),module/'NEW_CELLCHAT_DIRECTION_CONTEXT.tsv')
    python_probe="import json,sys,platform; from importlib import metadata; names=['SigXTalkpy','numpy','pandas','scipy','scikit-learn','torch','dhg']; out={'executable':sys.executable,'python':sys.version,'platform':platform.platform(),'packages':{}}; exec(\"for n in names:\\n try:out['packages'][n]=metadata.version(n)\\n except metadata.PackageNotFoundError:out['packages'][n]=None\"); print(json.dumps(out,indent=2))"
    actual_python=subprocess.run([cfg['python'],'-c',python_probe],capture_output=True,text=True,check=True)
    actual_python_data=json.loads(actual_python.stdout)
    save(actual_python_data,wd/'audit/PYTHON_ENVIRONMENT_ACTUAL_FINALIZE.json')
    environment_files=[wd/'audit/R_SESSION.txt',wd/'audit/PYTHON_ENVIRONMENT_ACTUAL_FINALIZE.json']+list(args.environment_record)
    for p in sorted((root/'envs').glob('*')):
        if p.is_file() and p.suffix.lower() in ['.json','.tsv','.txt','.yml','.yaml'] and p.stat().st_size<20*1024**2:environment_files.append(p)
    environment_files=list(dict.fromkeys(Path(p).resolve() for p in environment_files if Path(p).is_file()))
    artifacts=[]
    for folder in ['input','inputs','r_objects','hgnn_outputs','prs_outputs','audit','logs']:
        for p in sorted((wd/folder).rglob('*')):
            if p.is_file():artifacts.append(receipt(p))
    for p in [resolve('input_dir')/'counts.mtx.gz',resolve('input_dir')/'genes.tsv',resolve('input_dir')/'meta.tsv',resolve('pairing_tsv'),barcode_audit]:
        require(p.is_file(),'Required current input/pairing evidence file missing');artifacts.append(receipt(p))
    if cfg.get('state_thresholds_tsv') and resolve('state_thresholds_tsv').is_file():artifacts.append(receipt(resolve('state_thresholds_tsv')))
    write(pd.DataFrame(artifacts),module/'NEW_RUNTIME_ARTIFACT_HASHES.tsv')
    models=[r for r in artifacts if '/hgnn_outputs/' in r['path'].replace('\\','/') and r['path'].endswith('.pth')]
    require(bool(models) and all(r['bytes']>0 for r in models),'No actual saved new formal model weight artifact')
    patch=load(wd/'audit/PRS_RESOURCE_PATCH_AUDIT.json') if (wd/'audit/PRS_RESOURCE_PATCH_AUDIT.json').exists() else {'bounded_resource_wrapper_used':False}
    patchdiff=wd/'audit/PRS_resource_rng.patch'
    proof={'status':'OFFICIAL_NEW_END_TO_END_DONE' if completion['fidelity_specificity_status']=='DONE' else 'PRS_DONE_FID_SPEC_PARTIAL',
        'generated_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'generated_by':'finalize_new_sigxtalk.py','manually_entered_success':False,
        'source_repo':sourceproof['official_repository'],'source_commit':sourceproof['commit'],'pinned_source_file_checks':sourceproof,
        'R_api_and_body_audit':receipt(wd/'audit/R_API_AND_BODY_AUDIT.tsv'),'config':receipt(cp),'input_provenance':inp,
        'actual_input_cells':len(meta),'actual_independent_patients':int(meta.donor_id.nunique()),
        'actual_donor_group_counts':meta.groupby(['donor_id','bridge_group'],observed=True).size().rename('n_cells').reset_index().to_dict('records'),
        'barcode_identity_audit':receipt(barcode_audit),'environment':{'Rscript':cfg['rscript'],'Python':cfg['python'],
             'actual_session_and_environment_records':[receipt(p) for p in environment_files]},
        'training':{'project_id':inp['project_id'],'target_type':inp['target_type'],'seed':hg['seed'],'epochs_completed':hg['epochs_completed'],
            'started_utc':executions['hgnn']['started_utc'],'finished_utc':executions['hgnn']['finished_utc'],'exitcode':executions['hgnn']['exitcode'],
            'official_main_sha256':hg['official_main_sha256'],'pathways':receipt(pred),'tested_rows':len(full),'activated_rows':len(active)},
        'stage_executions':executions,'PRS':completion,'fidelity_specificity':{'status':completion['fidelity_specificity_status'],
            'error':completion.get('fidelity_specificity_error'),'definitions':'official matrix functions $all; separately declared raw-positive and fixed-filter universes'},
        'technical_smoke':{'completion':smoke,'input_provenance':smoke_input,'directory':str(smoke_dir),
             'completion_receipt':receipt(smoke_dir/'audit/NEW_R_PYTHON_R_COMPLETION.json'),'stage_executions':smoke_exec,
             'technical_subset_not_full_published_PBMC_reproduction':True},
        'PRS_resource_wrapper':patch,'actual_patch_diff':receipt(patchdiff) if patchdiff.exists() else None,
        'runtime_file_manifest':receipt(module/'NEW_RUNTIME_ARTIFACT_HASHES.tsv'),
        'model_artifacts':models,
        'new_outputs_derived_from_current_files_only':True,'legacy_rows_copied_into_NEW':False,
        'legacy_comparison':{'requested':args.legacy_root is not None,'files':legacy_manifest,'purpose':'private historical comparison only'},
        'original_PA_and_PRS_preserved':True,'PA_promotion_source':'original official Python pathways CSV; R cross-check tolerance1e-14 for text-export precision',
        'PRS_NA_not_zero':True,'counts':counts}
    save(proof,module/'SIGXTALK_EXECUTION_PROOF.json')
    envlines=['# Actual environment and resource patch audit','',
        'This file is generated from successful current execution receipts and their file hashes. Runtime artifacts are private.',
        '',f"Official repository: {sourceproof['official_repository']}; commit `{sourceproof['commit']}`.",
        f"R interpreter: `{cfg['rscript']}`. Python interpreter: `{cfg['python']}`.",
        'The exact actual installed R package versions are recorded in the following successful-run session file(s); no bootstrap YAML is treated as proof:','']
    envlines.extend('- '+str(p)+'; SHA256 '+sha(p) for p in environment_files)
    envlines+=['','Successful R session:','```', (wd/'audit/R_SESSION.txt').read_text(), '```','',
        'PRS resource wrapper: '+json.dumps(patch,ensure_ascii=True),
        'The worker cap, thread option and explicit task RNG seeds leave the imported response,500tree ranger regression,mtry3,permutation importance and filtering definitions unchanged. They specify the stochastic realization; historical numerical identity is not promised.',
        'The installed pinned R function bodies and unchanged official Python files were checked. SeuratObject5.4 was incompatible with internal Seurat5.2.1 slot calls; the actual successful-run versions above supersede the bootstrap recipe.',
        'The initial real PBMC2epoch attempt had zero PA>.75 paths and no PRS closure. It remains preserved. The separate same-context standard100epoch technical subset completed as recorded in SIGXTALK_EXECUTION_PROOF.json; neither attempt substitutes for the formal IVDD result.']
    (module/'OFFICIAL_ENV_LOCK_AND_PATCH_AUDIT.md').write_text('\n'.join(envlines)+'\n')
    methods=f'''# New official SigXTalk: actually executed Methods and Results

GSE230808 integer raw counts were reconstructed for {len(meta):,} frozen selected cells from {meta.donor_id.nunique()} paired patients ({inp['n_genes']:,} raw features); actual barcode, GSM, donor, state and raw-count audits are retained. A new Seurat object was log-normalized to10,000 counts per cell. Frozen NP GEP2 quartile identities were carried through at the donor level. CellChat was run on the actual groups and the complete saved object retained all directions; only {cfg['sender']}→{cfg['receiver']} was supplied to the main official Prepare_Input context. Original receptor complexes and their separately expanded components were preserved. The fixed receiver variance rule selected {inp['n_selected_targets']} targets at detection≥{cfg['target_selection']['minimum_detection']}, excluding MT-/RPL/RPS; {inp['prepared_target_genes']} remained after official expression/prior filtering. No candidate gene was inserted.

The pinned official SigXTalk implementation ({sourceproof['commit']}) generated new HGNN predictions using project `{inp['project_id']}`, seed{hg['seed']} and100 completed epochs with the official default architecture. Preparation used `{inp['actual_datatype']}`, expression threshold{cfg['prepare']['exp_threshold']}, CCC relative threshold{cfg['prepare']['CCC_threshold']} and Fisher threshold{cfg['prepare']['Fisher_threshold']}. The official PRS calculation used normalized receiver expression, coexpression cutoff{cfg['prs']['coexpression_cutoff']}, the original regression/importance definition, and a documented worker/thread/RNG resource wrapper. Original PA>{cfg['hgnn']['activation_threshold']} selected PRS inputs; Filter_results used the fixed {cfg['prs']['relative_max_threshold']} fraction-of-maximum PRS rule. Returned official fidelity/specificity matrix allocations were computed separately over all new positive PRS rows and the filtered universe when supported. These are conditional pathway allocations, not causal confidence or patient-level effects.

The new pooled model tested {len(full):,} receptor–SSC–target axes, of which {len(active):,} exceeded the fixed PA threshold. Official PRS returned {len(raw):,} positive rows and the fixed filter retained {len(filtered):,}; {int((active.PRS.isna()).sum()):,} activated rows had no returned positive PRS and remain NA. Official fidelity/specificity status was `{completion['fidelity_specificity_status']}`. There were {len(lr):,} exported AF→NP LR context rows and {len(paths):,} receptor-key-linked GEP2 context rows. The candidate table records PLOD2, CCN2, PRSS23, THBS1, LOX and TGFBI regardless of whether they were selected, activated or returned PRS. No historical pathways were used as new output.
'''
    (module/'NEW_SIGXTALK_METHODS_RESULTS.md').write_text(methods)
    (module/'NEW_SIGXTALK_LIMITATIONS.md').write_text('''# New official SigXTalk evidence boundary

This is one pooled computational model from six disease-source paired patients. Individual cells are nested observations. There is no donor-level communication disease-effect test, healthy-versus-disease inference, independent clinical validation or measured signaling transfer. Frozen NP GEP2 high/low groups are expression states, not independent patients or clinical grades.

CellChat probabilities, HGNN PA and PRS are separate algorithmic outputs. Algorithmic CellChat p values are not donor-level disease p values. Complexes are preserved; an ITGB1 component does not turn an integrin complex into a validated single receptor. Database-listed SSC intermediates are not assumed to be DNA-binding transcription factors. Joining LR context to receiver Receptor–SSC–target predictions by a receptor-component key does not demonstrate a ligand-specific four-layer causal chain.

Target-selection and expression/prior filtering restrict testability. Genes outside that universe have not been tested; PA failure and missing positive PRS are retained explicitly. Missing PRS is NA, never zero. Official fid/spe allocations depend on their declared denominator universe and do not establish causality. The fresh count-feature universe and installed versions may differ from historical processed inputs, so agreement with old row counts was not optimized. The PBMC run is a bounded technical subset and is excluded from biological results and Figure3. Wet-lab experiments remain prospective.
''')
    exported=[]
    for p in sorted(module.glob('NEW_*')):
        if p.is_file() and p.name!='NEW_DELIVERABLE_HASHES.tsv':exported.append(table_receipt(p))
    write(pd.DataFrame(exported),module/'NEW_DELIVERABLE_HASHES.tsv')
    print(json.dumps({'status':proof['status'],'proof':str(module/'SIGXTALK_EXECUTION_PROOF.json'),'counts':counts},indent=2))


if __name__=='__main__':main()
