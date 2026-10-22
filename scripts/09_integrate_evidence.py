#!/usr/bin/env python3
"""Integrate actual formal-new SigXTalk, fixed human and frozen spatial evidence.

No composite mechanism score or clinical-grade P value is created. Strict mode
requires a new formal R/HGNN/PRS completion chain. --allow-incomplete-b writes
explicitly PARTIAL drafts without reading partial/smoke/history as formal data.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd
import yaml

CANDIDATES = 'FN1 ITGB1 PLOD2 CCN2 PRSS23 THBS1 LOX TGFBI GJA1 IGFBP3'.split()
ROLES = {'FN1':'prespecified_sender_ligand_candidate',
         'ITGB1':'prespecified_integrin_receptor_complex_component',
         'PLOD2':'prespecified_primary_NP_response_readout',
         'CCN2':'program_target_or_response_readout', 'PRSS23':'program_target_or_response_readout',
         'THBS1':'comparison_ligand_or_target', 'LOX':'comparison_target_or_readout',
         'TGFBI':'comparison_target_or_readout', 'GJA1':'comparison_target_or_readout',
         'IGFBP3':'comparison_target_or_readout'}
DEFAULTS = dict(tu_candidate_table='01_human_programs/TU_CANDIDATE_GENE_BY_DONOR_IDENTITY.tsv',
    frozen_programs='01_human_programs/FROZEN_K5_PROGRAMS.tsv', spatial_input='data_cache/spatial_stage',
    human_dir='01_human_programs', sigxtalk_run='02_sigxtalk_official_new/formal_new',
    spatial_depth_dir='03_spatial_gse287499', integrated_out='04_integrated',
    manuscript_out='06_manuscript', figure_dir='05_figures', candidates=CANDIDATES,
    tu_np_identity='NP_matrix_marker_supported', tu_context_identity='immune_marker_supported')

def truth(value):
    if isinstance(value, pd.Series):
        return value.astype(str).str.lower().isin(['true','1','yes'])
    return str(value).lower() in ['true','1','yes']

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def read(path, required=True):
    path=Path(path)
    if not path.is_file():
        if required: raise FileNotFoundError(str(path))
        return pd.DataFrame()
    return pd.read_csv(path,sep='\t',float_precision='round_trip')

def jsread(path, required=False):
    path=Path(path)
    if path.is_file(): return json.loads(path.read_text(encoding='utf8'))
    if required: raise FileNotFoundError(str(path))
    return {}

def write(frame,path):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    frame.to_csv(path,sep='\t',index=False,na_rep='NA')

def text(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(value.strip()+'\n',encoding='utf8')

def fmt(value, digits=4):
    try: return f'{float(value):.{digits}g}' if np.isfinite(float(value)) else 'not estimable'
    except (ValueError,TypeError): return 'not estimable'

def joined(values):
    return ';'.join(sorted(set(str(x) for x in values if pd.notna(x) and str(x))))

def jrecords(frame,columns):
    return frame[columns].to_json(orient='records',double_precision=12)

def markdown_table(frame):
    def cell(value):
        return str(value).replace('|','/').replace('\n',' ') if pd.notna(value) else 'NA'
    rows=['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
    rows.extend('| '+' | '.join(cell(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None))
    return '\n'.join(rows)

def resolve_inputs(args):
    cfg=DEFAULTS.copy()
    if args.config:
        cfg.update(yaml.safe_load(args.config.read_text(encoding='utf8')) or {})
    root=args.project_root.resolve()
    for key in ['tu_candidate_table','spatial_input','frozen_programs','sigxtalk_config','github_proof']:
        value=getattr(args,key,None)
        if value is not None: cfg[key]=str(value)
    for key,value in list(cfg.items()):
        if key in ['tu_candidate_table','spatial_input','frozen_programs','human_dir','sigxtalk_run',
                   'spatial_depth_dir','integrated_out','manuscript_out','figure_dir','sigxtalk_config','github_proof']:
            cfg[key]=(root/value).resolve()
    return root,cfg

def formal_status(wd):
    """Completion means real formal disease counts, training and new official PRS."""
    completion=jsread(wd/'audit/NEW_R_PYTHON_R_COMPLETION.json')
    final=jsread(wd/'NEW_OFFICIAL_FINAL_STAGE_STATUS.json')
    hgnn=jsread(wd/'audit/HGNN_STATUS.json')
    provenance=jsread(wd/'audit/NEW_INPUT_PROVENANCE.json')
    reasons=[]
    if not truth(completion.get('complete_R_python_R',False)): reasons.append('new_R_python_R_closure_missing')
    if not truth(completion.get('no_old_pathways_or_PRS_used',False)): reasons.append('new_PRS_origin_unverified')
    if not truth(final.get('formal_design',False)) or not truth(final.get('new_result_not_legacy_relabel',False)):
        reasons.append('formal_new_final_status_missing')
    if hgnn.get('status')!='DONE' or not truth(hgnn.get('new_computation',False)) or hgnn.get('mode')!='formal':
        reasons.append('formal_HGNN_training_not_verified')
    if provenance.get('input_kind')!='real_GSE230808_raw_counts': reasons.append('formal_input_is_not_verified_GSE230808_counts')
    if provenance.get('started_from_old_pathway_output') is not False: reasons.append('raw_input_origin_flag_unverified')
    if completion.get('new_HGNN_predictions_sha256') != hgnn.get('prediction_sha256'):
        reasons.append('formal_prediction_audit_hashes_disagree')
    provenance_file=wd/'audit/NEW_INPUT_PROVENANCE.json'
    if provenance_file.is_file() and completion.get('new_input_provenance_sha256')!=sha(provenance_file):
        reasons.append('formal_input_provenance_hash_disagrees')
    for name in ['ALL_TESTED_NEW_AXES_WITH_PRS.tsv','FOCUS_LR_COMPLEX_ALL.tsv','LR_COMPLEX_COMPONENTS.tsv']:
        if not (wd/'tables'/name).is_file(): reasons.append('missing_'+name)
    return not reasons, reasons, completion,final,hgnn,provenance

def human_context(cfg):
    frozen=read(cfg['frozen_programs'])
    if not {'program','rank','gene'}.issubset(frozen.columns) or frozen.duplicated(['program','gene']).any():
        raise ValueError('Frozen program/rank table schema failed')
    raw=read(cfg['tu_candidate_table'])
    required={'donor_id','identity','gene','n_cells','fraction_detected',
              'mean_log1p_published_expression','feature_present'}
    if not required.issubset(raw): raise ValueError('Tu candidate table must use the disclosed donor/identity expression schema')
    if raw.duplicated(['donor_id','identity','gene']).any(): raise ValueError('Duplicate Tu donor/identity/gene')
    result=[]
    expected=set(raw.donor_id)
    stage_a=jsread(cfg['human_dir']/'STAGE_A_COMPUTATION_STATUS.json',required=True)
    if len(expected)!=int(stage_a['Tu_donors']):
        raise ValueError('Tu candidate donors do not match the actual Stage A donor count')
    for gene in cfg['candidates']:
        npdata=raw.loc[raw.gene.eq(gene)&raw.identity.eq(cfg['tu_np_identity'])].copy()
        bg=raw.loc[raw.gene.eq(gene)&raw.identity.eq(cfg['tu_context_identity']),
                   ['donor_id','mean_log1p_published_expression','fraction_detected']].rename(columns={
                       'mean_log1p_published_expression':'context_mean_log1p_published_expression',
                       'fraction_detected':'context_fraction_detected'})
        if set(npdata.donor_id)!=expected:
            raise ValueError('Tu candidate panel must include every donor for '+gene+'; unavailable features must be explicit')
        q=npdata.merge(bg,on='donor_id',how='left',validate='one_to_one')
        q['context_identity']=cfg['tu_context_identity']
        q['NP_minus_context_mean_expression']=q.mean_log1p_published_expression-q.context_mean_log1p_published_expression
        q['NP_gt_context_expression']=np.where(q.context_mean_log1p_published_expression.notna(),
                                               q.NP_minus_context_mean_expression>0,None)
        q['direction_scope']='within_donor_NP_identity_vs_context_identity;not_clinical_grade'
        q['clinical_grade_association']='UNTESTABLE_UNRESOLVED_TU_CROSSWALK'
        result.append(q)
    return frozen,raw,pd.concat(result,ignore_index=True)

def load_new_b(cfg,ready,completion):
    if not ready: return {}
    wd=cfg['sigxtalk_run'];td=wd/'tables'
    full=read(td/'ALL_TESTED_NEW_AXES_WITH_PRS.tsv')
    needed={'Receptor','TF','TG','pred_label','PRS','PRS_relative_max_filter_pass','method_origin',
            'source_project_id','source_run_started_from_legacy_predictions'}
    if not needed.issubset(full) or full.empty or full.duplicated(['Receptor','TF','TG']).any():
        raise ValueError('Formal-new complete axes schema or nonempty/identity check failed')
    if (full.method_origin.isna().any() or not full.method_origin.str.startswith('NEW_OFFICIAL_SIGXTALK').all() or
        not full.source_run_started_from_legacy_predictions.astype(str).str.lower().isin(['false','0']).all() or
        full.source_project_id.isna().any() or full.source_project_id.nunique()!=1):
        raise ValueError('Historical or unverified prediction source is forbidden in the new main analysis')
    if len(full)!=int(completion['all_tested_HGNN_rows']) or int(full.PRS.notna().sum())!=int(completion['positive_PRS_rows']):
        raise ValueError('Actual complete row tables disagree with formal completion audit')
    source_id=jsread(wd/'audit/NEW_INPUT_PROVENANCE.json',required=True).get('project_id')
    if not full.source_project_id.eq(source_id).all():
        raise ValueError('Axis source project does not match the actual formal input receipt')
    full['axis_source_row']=np.arange(len(full))+2
    full['axis_key']=full.Receptor.astype(str)+'::'+full.TF.astype(str)+'::'+full.TG.astype(str)
    lr=read(td/'FOCUS_LR_COMPLEX_ALL.tsv'); comp=read(td/'LR_COMPLEX_COMPONENTS.tsv')
    if not {'source','target','ligand','receptor','interaction_name','prob','pval'}.issubset(lr):
        raise ValueError('Real full sender/receiver LR complex schema missing')
    if set(lr.source)!={'AF_stromal'} or set(lr.target)!={'NP_GEP2_high'}:
        raise ValueError('Main LR context must be actual AF_stromal to NP_GEP2_high')
    lr['lr_source_row']=np.arange(len(lr))+2
    keys=['source','target','ligand','interaction_name']
    comp=comp.rename(columns={'receptor_complex':'receptor'})
    if not set(keys+['receptor','receptor_component']).issubset(comp):
        raise ValueError('Explicit receptor-complex component table unavailable')
    complex_keys=keys+['receptor']
    comp=comp.drop_duplicates(complex_keys+['receptor_component'])
    has_itgb1=comp.groupby(complex_keys).receptor_component.apply(lambda x:'ITGB1' in set(x)).rename('complex_contains_ITGB1').reset_index()
    lr=lr.merge(has_itgb1,on=complex_keys,how='left',validate='many_to_one')
    expanded=lr.merge(comp,on=complex_keys,how='inner',validate='one_to_many').rename(columns={'receptor_component':'Receptor'})
    cfg_b=jsread(cfg['sigxtalk_config']) if cfg.get('sigxtalk_config') else {}
    threshold=cfg_b.get('hgnn',{}).get('activation_threshold',.75)
    active=full.loc[full.pred_label>threshold]
    eligible=active.loc[active.PRS.notna()&truth(active.PRS_relative_max_filter_pass)]
    context=expanded.loc[expanded.pval<.05].merge(eligible,on='Receptor',how='inner',validate='many_to_many')
    context['join_scope']='receptor_key_context_only;not_ligand_specific_downstream_model'
    metrics=read(wd/'prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv',required=False)
    if not metrics.empty:
        metrics=metrics.rename(columns={'SSC':'TF','Target':'TG'})
        if metrics.duplicated(['Receptor','TF','TG']).any(): raise ValueError('Duplicate official fid/specificity row')
        context=context.merge(metrics,on=['Receptor','TF','TG'],how='left',validate='many_to_one',suffixes=('','_metric'))
        metric_columns=['Receptor','TF','TG']+[c for c in metrics.columns if c not in ['Receptor','TF','TG','PRS']]
        full=full.merge(metrics[metric_columns],on=['Receptor','TF','TG'],how='left',validate='one_to_one',sort=False)
    selection=read(td/'TARGET_GENE_SELECTION_ALL.tsv',required=False)
    testability=read(td/'CANDIDATE_TESTABILITY.tsv',required=False)
    prepare_lr=read(td/'PREPARE_LR_RECEIVER_FILTER_AUDIT.tsv',required=False)
    if len(prepare_lr):
        required_gate={'From','To','Weight','ligand_receiver_gene_present','receptor_receiver_gene_present'}
        if not required_gate.issubset(prepare_lr): raise ValueError('Original Prepare_Input LR gate audit schema incomplete')
        eligible=truth(prepare_lr.ligand_receiver_gene_present)&truth(prepare_lr.receptor_receiver_gene_present)
        prov=jsread(wd/'audit/NEW_INPUT_PROVENANCE.json',required=True)
        relative=prov.get('prepare_parameters',{}).get('CCC_threshold')
        if relative is None: raise ValueError('Actual official relative CCC threshold is absent from provenance')
        maximum=prepare_lr.loc[eligible,'Weight'].max()
        prepare_lr['passes_official_receiver_feature_gate']=eligible
        prepare_lr['passes_official_relative_CCC_gate']=eligible&(prepare_lr.Weight>=float(relative)*maximum)
        prepare_lr['CCC_relative_threshold_used']=relative
        prepare_lr['CCC_max_after_receiver_feature_gate']=maximum
    return dict(full=full,lr=lr,comp=comp,context=context,metrics=metrics,selection=selection,
                testability=testability,activation_threshold=threshold,config=cfg_b,
                fidelity_specificity_status=completion.get('fidelity_specificity_status','NOT_RECORDED'),
                prepare_lr=prepare_lr)

def target_evidence(gene,b,ready):
    missing=dict(n_tested_axes=np.nan,n_PA_pass_axes=np.nan,n_positive_PRS_axes=np.nan,n_filtered_PRS_axes=np.nan,
        max_original_PA=np.nan,max_original_PRS=np.nan,n_AF_LR_receptor_key_context_rows=np.nan,
        n_axes_with_official_filtered_universe_fid_spe=np.nan,
        official_fid_spe_assessability='NOT_ASSESSED_NEW_FORMAL_RUN_INCOMPLETE',
        target_testability='NOT_ASSESSED_NEW_FORMAL_RUN_INCOMPLETE',original_testability_reason='formal_new_completion_missing',
        representative_path_keys='',representative_path_details_json='[]')
    if not ready: return missing
    full=b['full'];q=full.loc[full.TG.eq(gene)]
    a=q.loc[q.pred_label>b['activation_threshold']]; p=a.loc[a.PRS.notna()]
    f=p.loc[truth(p.PRS_relative_max_filter_pass)]; ctx=b['context'].loc[b['context'].TG.eq(gene)]
    reason=''
    t=b['testability']
    if len(t) and 'gene' in t:
        own=t.loc[t.gene.eq(gene)]
        if len(own): reason=str(own.reason.iloc[0])
    selection=b['selection']; selected=None
    if len(selection) and {'gene','selected_target'}.issubset(selection):
        own=selection.loc[selection.gene.eq(gene)]
        selected=bool(truth(own.selected_target).any())
    if not len(q): status='NOT_IN_TESTED_TARGET_UNIVERSE'
    elif not len(a): status='TESTED_NO_PA_THRESHOLD_PASS'
    elif not len(p): status='PA_PASS_NO_RETURNED_OFFICIAL_PRS'
    elif not len(f): status='POSITIVE_PRS_NO_RELATIVE_FILTER_PASS'
    elif not len(ctx): status='RECEIVER_PATH_PASS_NO_AF_LR_RECEPTOR_KEY_CONTEXT'
    else: status='SOURCE_SUPPORTED_RECEPTOR_KEY_CANDIDATE_CONTEXT'
    if not len(q) and selected is False: reason=reason or 'not_selected_by_fixed_target_variance_rule'
    elif not len(q): reason=reason or 'not_tested_after_official_expression_or_prior_filter;inspect_target_and_prior_tables'
    elif not len(a): reason=reason or 'tested_but_no_original_PA_above_declared_threshold'
    elif not len(p): reason=reason or 'activated_axes_have_no_positive_official_PRS_returned;NA_not_zero'
    elif not len(f): reason=reason or 'positive_PRS_did_not_pass_official_relative_max_filter'
    else: reason=reason or 'tested;see_individual_original_axes'
    cols=['source','target','ligand','receptor','interaction_name','Receptor','TF','TG','pred_label','PRS',
          'axis_key','axis_source_row','lr_source_row','complex_contains_ITGB1']
    cols += [x for x in ['fidelity_target_all','specificity_receptor_all','normalization_universe'] if x in ctx]
    rep=ctx.sort_values(['PRS','pred_label','prob','axis_key','interaction_name'],
                        ascending=[False,False,False,True,True],kind='mergesort').head(3)
    metrics=b.get('metrics',pd.DataFrame())
    nmetrics=int(metrics.TG.eq(gene).sum()) if len(metrics) else 0
    metric_reason=('AVAILABLE_FOR_RETURNED_FILTERED_PATHS_WITH_EXPLICIT_UNIVERSE' if nmetrics else
                   'NO_RETAINED_PATH_FOR_THIS_TARGET_IN_FILTERED_OFFICIAL_METRIC_UNIVERSE' if len(metrics) else
                   'OFFICIAL_METRICS_UNAVAILABLE:'+b.get('fidelity_specificity_status','NOT_RECORDED'))
    return dict(n_tested_axes=len(q),n_PA_pass_axes=len(a),n_positive_PRS_axes=len(p),
        n_filtered_PRS_axes=len(f),max_original_PA=q.pred_label.max() if len(q) else np.nan,
        max_original_PRS=p.PRS.max() if len(p) else np.nan,n_AF_LR_receptor_key_context_rows=len(ctx),
        n_axes_with_official_filtered_universe_fid_spe=nmetrics if len(metrics) else np.nan,official_fid_spe_assessability=metric_reason,
        target_testability=status,original_testability_reason=reason,
        representative_path_keys=joined(rep.axis_key) if len(rep) else '',
        representative_path_details_json=jrecords(rep,cols) if len(rep) else '[]')

def build_matrix(cfg,frozen,tu,mouse,depth,b,ready):
    rows=[]
    for gene in cfg['candidates']:
        f=frozen.loc[frozen.program.eq('GEP2')&frozen.gene.eq(gene)]
        h=tu.loc[tu.gene.eq(gene)]; m=mouse.loc[mouse.human_gene.eq(gene)]
        np_m=m.loc[m.region.eq('NP_CORE')]; af_m=m.loc[m.region.eq('AF_CORE')]
        detected=lambda q:int((pd.to_numeric(q.n_detected_bins,errors='coerce')>0).sum())
        row=dict(gene=gene,prespecified_role=ROLES.get(gene,'candidate_comparison'),
            GEP2_frozen_member=bool(len(f)),GEP2_frozen_rank=int(f['rank'].iloc[0]) if len(f) else np.nan,
            Tu_NP_donors=len(h),Tu_NP_feature_assessable_donors=int(truth(h.feature_present).sum()),
            Tu_NP_donors_detected=int((h.fraction_detected>0).sum()),
            Tu_NP_fraction_detected_min=h.fraction_detected.min(),Tu_NP_fraction_detected_max=h.fraction_detected.max(),
            Tu_NP_gt_context_donors=int(truth(h.NP_gt_context_expression).sum()),
            Tu_NP_vs_context_assessable_donors=int(h.context_mean_log1p_published_expression.notna().sum()),
            Tu_direction_comparator=cfg['tu_context_identity'],
            Tu_per_donor_expression_json=jrecords(h,['donor_id','fraction_detected','mean_log1p_published_expression',
                'feature_present','context_fraction_detected','context_mean_log1p_published_expression','NP_gt_context_expression']),
            clinical_grade_test='UNTESTABLE_NO_GRADE_P_VALUE',
            mouse_NP_discs_assessable=int(truth(np_m.native_feature_assessable).sum()),
            mouse_AF_discs_assessable=int(truth(af_m.native_feature_assessable).sum()),
            mouse_NP_discs_detected=detected(np_m),mouse_AF_discs_detected=detected(af_m),
            mouse_NP_detected_bins=int(np_m.n_detected_bins.sum()),mouse_AF_detected_bins=int(af_m.n_detected_bins.sum()),
            mouse_NP_total_UMI=np_m.total_raw_UMI.sum(min_count=1),mouse_AF_total_UMI=af_m.total_raw_UMI.sum(min_count=1),
            mouse_per_disc_region_json=jrecords(m,['capture_id','disc_id','region','n_bins','fraction_detected_bins',
                'n_detected_bins','total_raw_UMI','native_feature_assessable']),
            GEP2_depth_outcome=depth['outcome'],GEP2_depth_original_eligible_mean=depth['mean_paired_delta_original_same_eligible_bins'],
            GEP2_depth100_mean=depth['mean_paired_delta_depth100'],candidate_expression_depth_corrected=False,
            genotype_inference='NOT_PERFORMED_CAPTURE_GENOTYPE_CONFOUNDED',
            biological_unit='human_donor;mouse_disc_nested_in_capture;cells_bins_not_biological_n')
        row.update(target_evidence(gene,b,ready))
        row.update(AF_ligand_LR_rows=np.nan,AF_ligand_LR_p005_rows=np.nan,
                   receptor_complex_component_LR_rows=np.nan,FN1_ITGB1_complexes='',
                   role_specific_joined_context_rows=np.nan,role_specific_AF_ligands='',
                   role_specific_LR_receptor_complexes='',role_specific_receiver_Receptor_keys='',
                   role_specific_database_SSC_symbols='',role_specific_receiver_targets='',
                   role_specific_path_keys='',official_receiver_feature_CCC_gate_pairs=np.nan,
                   official_LR_gate_passing_receptor_keys='',
                   official_receiver_feature_CCC_gate_reason='NOT_ASSESSED_NEW_FORMAL_RUN_INCOMPLETE')
        negatives=['clinical_grade_crosswalk_unresolved','one_spatial_capture_per_genotype',
                   'GEP2_spatial_depth_sensitive' if depth.get('observed_mean_difference_attenuated') else 'single_limited_depth_check']
        if ready:
            ligand=b['lr'].loc[b['lr'].ligand.eq(gene)]
            component=b['comp'].loc[b['comp'].receptor_component.eq(gene)]
            row.update(AF_ligand_LR_rows=len(ligand),AF_ligand_LR_p005_rows=int((ligand.pval<.05).sum()),
                       receptor_complex_component_LR_rows=len(component))
            row['FN1_ITGB1_complexes']=joined(b['lr'].loc[b['lr'].ligand.eq('FN1')&truth(b['lr'].complex_contains_ITGB1),'receptor'])
            role_context=(b['context'].loc[b['context'].ligand.eq('FN1')] if gene=='FN1' else
                          b['context'].loc[truth(b['context'].complex_contains_ITGB1)] if gene=='ITGB1' else
                          b['context'].loc[b['context'].TG.eq(gene)])
            row.update(role_specific_joined_context_rows=len(role_context),
                       role_specific_AF_ligands=joined(role_context.ligand),
                       role_specific_LR_receptor_complexes=joined(role_context.receptor),
                       role_specific_receiver_Receptor_keys=joined(role_context.Receptor),
                       role_specific_database_SSC_symbols=joined(role_context.TF),
                       role_specific_receiver_targets=joined(role_context.TG),
                       role_specific_path_keys=joined(role_context.axis_key))
            gate=b['prepare_lr']
            if len(gate):
                own_gate=gate.loc[gate.From.eq(gene)] if gene!='ITGB1' else gate.loc[gate.To.eq(gene)]
                passed=own_gate.loc[truth(own_gate.passes_official_relative_CCC_gate)]
                row['official_receiver_feature_CCC_gate_pairs']=len(passed)
                row['official_LR_gate_passing_receptor_keys']=joined(passed.To)
                row['official_receiver_feature_CCC_gate_reason']=('PASSES_LR_RECEIVER_FEATURE_AND_RELATIVE_CCC_GATES;not_full_pathway_causality' if len(passed) else
                    'NOT_RETURNED_IN_SIGNIFICANT_LR_HELPER' if own_gate.empty else
                    'REMOVED_BY_OFFICIAL_RECEIVER_FEATURE_OR_RELATIVE_CCC_GATE')
            else: row['official_receiver_feature_CCC_gate_reason']='ORIGINAL_PREPARE_LR_GATE_AUDIT_UNAVAILABLE'
        if gene=='FN1': negatives.append('FN1_detected_in_NP_not_AF_exclusive')
        if gene=='ITGB1': negatives.append('integrin_complex_component_not_single_receptor_or_complete_functional_complex')
        if gene=='PRSS23': negatives.append('sparse_mouse_core_detection_does_not_establish_absence_in_human')
        if gene=='CCN2' and detected(af_m)>detected(np_m): negatives.append('mouse_AF_weighted_detection')
        if not ready: level='PENDING_NEW_FORMAL_SIGXTALK;EXPRESSION_CONTEXT_ONLY'
        elif row['role_specific_joined_context_rows']:
            level='PREDICTED_AF_RECEPTOR_KEY_CONTEXT_WITH_RECEIVER_PATH'
        elif gene=='FN1' and row['AF_ligand_LR_p005_rows']:
            level='AF_LIGAND_LR_PREDICTION_WITHOUT_TARGET_PATH_CLAIM'
        elif gene=='ITGB1' and row['receptor_complex_component_LR_rows']:
            level='RECEPTOR_COMPLEX_COMPONENT_CONTEXT'
        else: level='EXPRESSION_CONTEXT_WITHOUT_JOINED_NEW_FORMAL_PATH'
        row['evidence_level']=level
        row['key_counterevidence']=';'.join(negatives)
        row['experimental_feasibility']='future_factorial_FN1_exposure_by_ITGB1_knockdown;requires_reagent_and_viability_checks' if gene in ['FN1','ITGB1'] else 'future_response_readout;no_experimental_result_generated'
        rows.append(row)
    return pd.DataFrame(rows)

def verified_release(path):
    proof=jsread(path) if path else {}
    good=(proof.get('verified_public') is True and
          re.fullmatch(r'https://github\.com/[^/\s]+/[^/\s]+',str(proof.get('repo_url',''))) and
          re.fullmatch(r'[0-9a-fA-F]{40}',str(proof.get('head_sha',''))) and proof.get('verified_utc'))
    return bool(good),proof

def manuscript(cfg,matrix,ready,reasons,completion,final,hgnn,provenance,b,depth):
    out=cfg['manuscript_out'];integ=cfg['integrated_out'];a=cfg['human_dir'];figs=cfg['figure_dir']
    ast=jsread(a/'STAGE_A_COMPUTATION_STATUS.json',required=True)
    brief=read(a/'ALL_PROGRAMS_BRIEF.tsv')
    support=read(a/'GEP2_DONOR_SUPPORT.tsv')
    tu_support=support.loc[support.cohort.astype(str).str.startswith('Tu')]
    gep2=brief.loc[brief.program.eq('GEP2')].iloc[0]
    intro=("We evaluated a frozen human nucleus pulposus (NP) transcriptional program using independent donor-level "
           "expression, newly inferred annulus fibrosus (AF)-to-NP communication context, and fixed adult mouse "
           "anatomical expression cores. These evidence layers address different questions and were not combined "
           "into a numerical mechanism score. No new wet-laboratory result was generated.")
    human_methods=("The five original 25-gene cNMF programs were retained without refitting or replacement of their "
      "gene names. Human data were obtained from GSE230808 and GSE165722. The independent Tu matrix contains "
      "published normalized expression, which was retained as processed expression rather than relabeled as raw "
      "UMI counts. Sparse parsing, QC, and independent marker-support labels were applied using the executed "
      "human workflow. QC required at least 200 detected features, positive expression sum and mitochondrial "
      "fraction no greater than 0.20. Authoritative HGNC identities were accepted only when a unique released "
      "feature matched; absent or ambiguous genes remained unavailable. Official decoupler 2.2.0 AUCell used "
      "n_up=1500 and tmin=5, with implementation-specific feature ranking, unavailable-target pruning and effective "
      "denominators retained in the audit. Discovery UCell scores used the recorded original scoring workflow. "
      "AUCell and UCell numerical scales were not compared across cohorts. Donor medians and within-donor identity "
      "comparisons were descriptive; individual cells were not independent patients. The unresolved GEO I–IV versus "
      "paper Pfirrmann II–V crosswalk precluded clinical-grade testing.")
    human_methods+=(f" The executed summary recorded discovery score origin {ast.get('Swahn_score_origin','not recorded')} "
                    f"and Tu score method {ast.get('Tu_score_method','not recorded')}. Cached score reuse is distinguished "
                    "from new cNMF training; neither frozen projection nor table regeneration refits cNMF.")
    b_methods="New formal SigXTalk Methods are pending completion of the verified formal R/HGNN/PRS chain; technical smoke and historical predictions are excluded from these main results."
    b_results="The new formal SigXTalk analysis was not complete at integration time. No new formal pathway count, positive candidate decision, or Figure3 result is claimed. Pending checks: "+'; '.join(reasons)+'.'
    if ready:
        prep=provenance.get('prepare_parameters',{})
        b_methods=(f"The newly executed official SigXTalk workflow started from real GSE230808 raw counts, with "
          f"{provenance.get('n_cells','not recorded')} selected cells and {provenance.get('n_receiver','not recorded')} "
          f"NP_GEP2_high receivers across {provenance.get('independent_patient_count','not recorded')} donor identifiers. "
          "CellChat inference retained full LR results and explicitly restricted the main sender/receiver context to "
          "AF_stromal→NP_GEP2_high. Receptor complexes were preserved and component membership was recorded separately. "
          f"Prepare_Input used datatype {provenance.get('actual_datatype','not recorded')} and recorded parameters "
          f"{json.dumps(prep,sort_keys=True)}. Target selection followed the recorded rule: {provenance.get('target_rule','not recorded')}. "
          f"Gene filtering left {provenance.get('RTF_rows_after_gene_filter','not recorded')} RTF and "
          f"{provenance.get('TFT_rows_after_gene_filter','not recorded')} TFT prior rows, with "
          f"{provenance.get('prepared_target_genes','not recorded')} prepared target genes. The unchanged official Python "
          f"HGNN entry point completed {hgnn.get('epochs_completed','not recorded')} epochs, seed {hgnn.get('seed','not recorded')}. "
          f"Original activation probability (PA) greater than {b['activation_threshold']} selected input axes for official PRS. "
          f"The official PRS expression/coexpression cutoff was {completion.get('coexpression_cutoff','not recorded')}, "
          f"and Filter_results used the recorded relative-maximum threshold {completion.get('optional_relative_max_filter_threshold','not recorded')}. "
          "Missing/unreturned PRS remained NA rather than zero. Official fidelity and specificity, where available, "
          "used explicitly recorded raw or filtered pathway normalization universes. LR probability, PA, PRS and "
          "fidelity/specificity were kept as separate quantities and were not interpreted as patient-level disease "
          "P values. Receiver Receptor–SSC–target predictions were linked to source LR records by receptor component "
          "keys only. Database-listed SSC intermediates were not automatically treated as experimentally verified "
          "transcription factors, and a receptor-key join does not assign ligand-specific downstream causality.")
        session=cfg['sigxtalk_run']/'audit/R_SESSION.txt'
        if session.is_file():
            s=session.read_text(encoding='utf8',errors='replace')
            version_matches=re.findall(r'\b(SeuratObject|Seurat|CellChat|SigXTalkR|Matrix|ranger|tidymodels|tidyfst|NMF)_([0-9][A-Za-z0-9.\-]*)',s)
            versions=dict(version_matches)
            r_version=next((line for line in s.splitlines() if line.startswith('R version ')),'R version recorded in session receipt')
            b_methods+=' Runtime session: '+r_version+'; '+', '.join(k+' '+v for k,v in sorted(versions.items()))+'.'
        patch=jsread(cfg['sigxtalk_run']/'audit/PRS_RESOURCE_PATCH_AUDIT.json')
        if patch:
            b_methods+=(f" The declared bounded-resource PRS wrapper used {patch.get('workers','not recorded')} workers "
                        f"and explicit seed {patch.get('seed','not recorded')}; the source-attributed patch changed worker/thread "
                        "limits and random-seed specification while retaining the recorded official forest response, "
                        "tree count, mtry, importance definition and collection threshold. Numerical identity to historical "
                        "stochastic PRS values is not assumed.")
        b_results=(f"The new formal R→HGNN→PRS chain generated {len(b['full']):,} tested receiver axes, "
          f"{completion['active_PA_rows']:,} axes above the declared PA threshold, {completion['positive_PRS_rows']:,} "
          f"positive returned official PRS rows, and {completion['filtered_PRS_rows']:,} relative-filtered rows. "
          f"The original focused AF-to-NP LR table contained {len(b['lr']):,} records, of which "
          f"{int((b['lr'].pval<.05).sum()):,} met the recorded CellChat p<0.05 display rule. The run class was "
          f"{final.get('run_class')}; official fidelity/specificity status was {completion.get('fidelity_specificity_status')}. "
          "All tested and negative/NA rows were retained. No historical pathways or PRS values were substituted.")
    spatial_methods=("GSE287499 supplied two six-month-old mouse captures. Original 16-µm counts were aggregated "
      "sparsely in non-overlapping global 3×3 grid groups to 48 µm, requiring all nine positions, all nine count "
      "barcodes and all nine eligible in-tissue flags. Original cKO positions outside the cropped H&E image were "
      "flagged and excluded from aggregation without moving coordinates. QC required at least 20 UMIs and 20 "
      "detected symbols. Eight provisional disc structures had eight NP and sixteen bilateral AF core polygons "
      "frozen before program or candidate expression review; no CEP region was assigned. Capture-specific 48-µm "
      "half-diagonal margins and artifact exclusions were used. Independent Krt19/T/Krt8/Acan and Chad/Clec3a "
      "panels assessed directional structural-marker compatibility after freezing, without expert pathological "
      "validation. Strict one-to-one human/mouse mapping used a frozen source-attributed MGI snapshot and targeted "
      "Ensembl stable IDs or exact MGI symbols. Program scores were means of log(1+CP10k) over the same fixed "
      "available orthologs, including zeros. The predeclared depth check first quantified score–UMI associations, "
      "then sampled exactly 100 molecules without replacement across all genes in eligible bins with ≥100 original "
      "UMIs, one fixed draw per bin, seed 20261022. At least five retained bins per compartment and three paired "
      "discs per capture were required. No threshold search or additional seeds were used. Original count, image, "
      "ROI, coordinate and score checksums/reproductions were audited. Discs remained nested in capture; one "
      "capture per genotype prevented genotype inference, and bins were not independent animals.")
    nnp=int(matrix.mouse_NP_discs_assessable.max());naf=int(matrix.mouse_AF_discs_assessable.max())
    spatial_results=(f"The original anatomical NP-favored GEP2 pattern was sensitive to sequencing depth. "
      f"{depth['adequate_paired_discs']} disc pairs met the fixed coverage safeguard and "
      f"{depth['NP_greater_AF_after_depth100']} retained NP>AF after sampling. Among these same adequately covered "
      f"discs, the mean NP-minus-AF contrast on original eligible bins was "
      f"{fmt(depth['mean_paired_delta_original_same_eligible_bins'],5)} and became "
      f"{fmt(depth['mean_paired_delta_depth100'],5)} at 100 UMIs. The declared outcome was {depth['outcome']}. "
      "This attenuation weakens the original NP-localization interpretation; it does not establish absence of the "
      "program or independent disease replication. Original candidate measurements were not depth-corrected.")
    original_programs=read(cfg['spatial_input']/'NP_PROGRAM_BY_REGION.tsv',required=False)
    if len(original_programs) and {'program','disc_id','region','mean_program_score'}.issubset(original_programs):
        g5=original_programs.loc[original_programs.program.eq('GEP5')].pivot(
            index=['capture_id','disc_id'],columns='region',values='mean_program_score')
        if {'NP_CORE','AF_CORE'}.issubset(g5):
            g5=g5.dropna(subset=['NP_CORE','AF_CORE'])
            spatial_results+=(f" Original GEP5 scores showed AF>NP in {int((g5.AF_CORE>g5.NP_CORE).sum())}/{len(g5)} "
                              "disc pairs. This opposite compartment context is retained and is not relabeled as NP replication.")
    association=read(cfg['spatial_depth_dir']/'GEP2_SCORE_UMI_ASSOCIATION.tsv',required=False)
    if len(association):
        allqc=association.loc[association.scope.eq('ALL_ORIGINAL_QC_BINS')]
        spatial_results+=' Descriptive score–UMI Spearman rho across all original QC bins was '+', '.join(
            str(r.capture_id)+' '+fmt(r.spearman_rho_score_vs_UMI) for r in allqc.itertuples())+'. No bin-level P value was calculated.'
    retention=read(cfg['spatial_depth_dir']/'DEPTH_RETENTION_BY_CAPTURE_REGION.tsv',required=False)
    if len(retention):
        spatial_results+=' Retained bins were '+', '.join(
            f'{r.capture_id} {r.region}: {int(r.retained_n_bins)}/{int(r.original_n_bins)} ({100*r.retention_fraction:.1f}%)'
            for r in retention.itertuples())+'. Eligibility therefore retained a smaller subset of low-depth AF than NP.'
    gene_results=[]
    for gene in ['FN1','ITGB1','PLOD2','CCN2','PRSS23','THBS1','LOX','TGFBI']:
        row=matrix.set_index('gene').loc[gene]
        gene_results.append(f"{gene}: Tu NP detection in {int(row.Tu_NP_donors_detected)}/{int(row.Tu_NP_donors)} donors "
          f"(fraction {fmt(row.Tu_NP_fraction_detected_min)}–{fmt(row.Tu_NP_fraction_detected_max)}), mouse NP "
          f"{int(row.mouse_NP_discs_detected)}/{int(row.mouse_NP_discs_assessable)} and AF "
          f"{int(row.mouse_AF_discs_detected)}/{int(row.mouse_AF_discs_assessable)} detected cores; new target status "
          f"{row.target_testability}.")
        if gene=='ITGB1':
            gene_results[-1]+=' This status refers to its target role; ITGB1 was modeled as a receptor-complex component.'
    human_results=(f"The regenerated human tables contained {ast.get('Tu_donors','not recorded')} Tu donors and "
      f"{ast.get('Swahn_donors','not recorded')} discovery donors, with {ast.get('Tu_cells_qc_pass','not recorded')} "
      f"Tu QC cells. GEP2 had {int(gep2.mapped_genes_min)}–{int(gep2.mapped_genes_max)} mapped genes per Tu sample "
      f"and NP donor-median scores from {fmt(gep2.NP_donor_median_min)} to {fmt(gep2.NP_donor_median_max)}. "
      f"NP donor medians exceeded the declared within-donor identity context in {int(truth(tu_support.NP_gt_context).sum())}/"
      f"{int(tu_support.context_donor_median.notna().sum())} assessable Tu donors (of {len(tu_support)} total). Expression review status was {gep2.expression_interpretation}; "
      "clinical-grade association remained untestable. These within-cohort observations do not compare absolute "
      "AUCell versus UCell scales or infer a cross-series disease effect.")
    refs=("Method sources: [SigXTalk, Nature Communications 2025](https://doi.org/10.1038/s41467-025-61149-7); "
      "[official decoupler AUCell](https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.aucell.html); "
      "[UCell](https://doi.org/10.1016/j.csbj.2021.06.043). Dataset accessions and exact source/commit/version "
      "receipts accompany the computational audit. No method is presented as newly invented here.")
    text(out/'METHODS_COMPUTATION_ACTUALLY_EXECUTED.md', '# Computational Methods: actually executed scope\n\n'+intro+'\n\n'+human_methods+'\n\n'+b_methods+'\n\n'+spatial_methods+'\n\n'+refs)
    text(out/'RESULTS_COMPUTATION_ACTUALLY_EXECUTED.md','# Computational Results\n\n'+human_results+'\n\n'+b_results+'\n\n'+spatial_results+'\n\n'+'\n\n'.join(gene_results)+'\n\nNo FN1→ITGB1→PLOD2 causal chain, new wet-laboratory validation or Figure5 experimental result was demonstrated.')
    limitations=("Clinical-grade mapping in Tu remains unresolved, and only processed normalized expression was "
      "available there. Cohorts differ in platform, tissue sampling, age and clinical composition; a pooled "
      "cross-series disease effect is not estimable from this design. Some annotation markers overlap GEP5, "
      "reducing independent specificity assessment for that program. Scores from AUCell and UCell are method-specific. "
      "Communication inference pools the disclosed disease-internal sender/receiver cells and does not provide "
      "donor-specific intervention effects. FN1 is expressed in NP as well as AF; ITGB1 is an integrin-complex "
      "component. PA, PRS and official fidelity/specificity are model quantities, not causal confidence or clinical "
      "P values. Receptor-key joins do not establish which ligand induced a receiver pathway. Mouse core "
      "measurements show sparse Prss23 and AF-weighted Ccn2, and the GEP2 contrast weakens after fixed-depth "
      "sampling. Low AF retention and one fixed random draw limit the depth check. Geometry is provisional and "
      "marker-compatible without expert pathology validation; there are no CEP ROIs or animal-to-disc labels. "
      "Genotype is completely confounded with capture. Neither spatial adjacency, directional molecular transfer "
      "nor functional FN1–integrin signaling was measured. Proposed perturbations remain future work.")
    if not ready: limitations+=' New formal SigXTalk remains pending; the integrated candidate decision cannot be finalized.'
    text(out/'EVIDENCE_LIMITATIONS.md','# Evidence limitations\n\n'+limitations)
    text(integ/'NEGATIVE_AND_DISCORDANT_EVIDENCE.md','# Negative and discordant evidence\n\n'+limitations+'\n\n'+markdown_table(matrix[['gene','target_testability','original_testability_reason','key_counterevidence']]))
    release_ok,proof=verified_release(cfg.get('github_proof'))
    if release_ok:
        code=(f"The computational scripts are publicly accessible at [{proof['repo_url']}]({proof['repo_url']}) "
              f"(commit `{proof['head_sha']}`, publicly verified {proof['verified_utc']}). The repository documents "
              "the executed cached-input workflow and the available raw-data entry points; exact reproduction "
              "limits and third-party dependencies are stated separately.")
    else:
        code=("CODE_RELEASE_STATUS=BLOCKED_UNVERIFIED. A public repository URL and commit have not been independently "
              "verified in the supplied release receipt. The prepared code package is not evidence of public release. "
              "This is a draft availability statement and must be replaced only after verified public access.")
    text(out/'CODE_AVAILABILITY_FOR_MANUSCRIPT.md','# Code Availability\n\n'+code)
    data=("The reused human single-cell datasets are available from GEO under "
      "[GSE230808](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230808) and "
      "[GSE165722](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165722). The adult mouse spatial dataset, "
      "including counts, positions and H&E images, is available under "
      "[GSE287499](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287499). GSE229711 is retained only where "
      "explicitly identified in frozen discovery provenance and does not constitute an additional disease-effect "
      "comparison in these analyses. The computational review package contains derived donor/disc summaries, "
      +("complete new official prediction tables, " if ready else "explicitly pending new-formal status rather than substituted predictions, ")+
      "figure source data and provenance manifests. Large source matrices "
      "and model intermediates are not duplicated in the public code repository; public source links and preparation "
      "instructions are provided. No new sequencing or wet-laboratory dataset was generated in this computation.\n\n"
      "Repository/citation actions: the review package is a local audit artifact, not yet a permanently deposited "
      "data accession. A stable public location for any newly deposited derived dataset must be supplied before "
      "submission; no new accession or data DOI is claimed here.\n\n中文核对：GEO 原始来源为公开重用；本地审核包不等于正式数据存储库。新衍生数据的长期存储位置仍需按真实发布状态填写；不虚构 DOI。")
    text(out/'DATA_AVAILABILITY_FOR_MANUSCRIPT.md','# Data Availability\n\n'+data)
    legends=[]
    for title,patterns in [('Figures1–2',['**/*LEGEND*.md','**/*legend*.md']),('Figure3',['**/*Figure3*LEGEND*.md','**/*FIGURE3*LEGEND*.md']),('Figure4',['Figure4_LEGEND.md'])]:
        paths=[]
        if title=='Figure3' and not ready:
            legends.append('## Figure3\n\nPending verified new formal output. No new formal communication figure is claimed.');continue
        for pattern in patterns:
            paths.extend(figs.glob(pattern))
        paths=sorted(set(p for p in paths if p.is_file() and
            ((title=='Figure4' and 'Figure4' in p.name) or (title=='Figure3' and ('Figure3' in p.name or 'FIGURE3' in p.name)) or
             (title=='Figures1–2' and 'Figure3' not in p.name and 'FIGURE3' not in p.name and 'Figure4' not in p.name))))
        if paths:
            legends.append('## '+title+'\n\n'+'\n\n'.join(p.read_text(encoding='utf8') for p in paths))
        else:
            legends.append('## '+title+'\n\nActual panel legend file not located in the supplied figure directory; copy the verified panel legend before submission. No panel details are invented here.')
    legends.append('## Figure5\n\nReserved for future perturbation experiments; no experimental figure was generated.')
    text(out/'FIGURE_LEGENDS_COMPUTATION.md','# Computational figure legends\n\n'+'\n\n'.join(legends))
    text(out/'SUPPLEMENTARY_TABLES_AND_DATA_STRUCTURE.md',
      '# Supplementary tables and source data\n\n'
      '01_human_programs provides frozen program provenance, gene mappings, donor and identity summaries, effective AUCell coverage and cohort/grade ambiguity. '
      '02_sigxtalk_official_new/formal_new retains raw-count input audits, all focused LR complexes, explicit complex components, every tested HGNN axis, activated axes, returned and filtered official PRS, original fidelity/specificity and completion/version/hash receipts. '
      '03_spatial_gse287499 provides fixed ROI provenance, per-disc compartment/depth summaries, all original and sampled values, retention, score–UMI association and unadjusted candidate detectability. '
      '04_integrated links candidate evidence with exact receptor/SSC/target row keys, official unmodified values and counterevidence. '
      '05_figures contains figure-specific TSV sources and matching PNG/PDF/SVG. Cells and bins are technical observations; human donors and discs nested in captures define the disclosed descriptive units. '
      'Private absolute-path runtime logs and large count matrices are not public-code payloads. No missing PRS is converted to zero and no model or synthetic technical fixture is reported as a scientific observation.')
    write(pd.DataFrame([
      {'new_figure':'Figure1','purpose':'source/cohort/QC/program provenance','old_material_policy':'frozen computational provenance only; old experiments excluded'},
      {'new_figure':'Figure2','purpose':'fixed-program independent donor expression context','old_material_policy':'new regenerated summaries; no clinical grade replication'},
      {'new_figure':'Figure3','purpose':'new formal official LR and receiver pathway predictions' if ready else 'PENDING_NEW_FORMAL_RESULT','old_material_policy':'historical pathway/PRS rows forbidden in main results'},
      {'new_figure':'Figure4','purpose':'unchanged geometry/original context plus single depth sensitivity','old_material_policy':'frozen original counts/ROIs retained; attenuation now explicit'},
      {'new_figure':'Figure5','purpose':'future experiment only; not generated','old_material_policy':'old THBS1 experiments are not new validation'}]),out/'OLD_TO_NEW_FIGURE_CROSSWALK.tsv')
    return release_ok

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,required=True)
    ap.add_argument('--config',type=Path)
    for flag in ['tu-candidate-table','spatial-input','frozen-programs','sigxtalk-config','github-proof']:
        ap.add_argument('--'+flag,type=Path)
    ap.add_argument('--allow-incomplete-b',action='store_true',help='Write explicitly partial drafts; never consume partial B as formal results')
    args=ap.parse_args();root,cfg=resolve_inputs(args)
    if not set(CANDIDATES).issubset(set(cfg['candidates'])):
        raise ValueError('All ten prespecified comparison candidates must remain in the evidence matrix')
    ready,reasons,completion,final,hgnn,provenance=formal_status(cfg['sigxtalk_run'])
    previous=jsread(cfg['integrated_out']/'INTEGRATION_STATUS.json')
    if not ready and previous.get('status')=='DONE_NEW_FORMAL_INTEGRATION':
        raise RuntimeError('Refuse to overwrite completed integration with a pending draft; inspect changed formal inputs')
    if not ready and not args.allow_incomplete_b:
        raise RuntimeError('New formal SigXTalk completion is not verified: '+'; '.join(reasons)+
                           '. Use --allow-incomplete-b only for explicitly pending drafts.')
    frozen,tu_all,tu=human_context(cfg)
    depth=jsread(cfg['spatial_depth_dir']/'DEPTH_ANALYSIS_STATUS.json',required=True)
    if depth.get('status')!='DONE_SINGLE_PREDECLARED_DEPTH_CHECK' or not depth.get('protected_old_inputs_unchanged'):
        raise ValueError('Real frozen spatial depth analysis is incomplete or old inputs changed')
    mouse=read(cfg['spatial_depth_dir']/'CANDIDATE_GENES_PER_DISC_REGION.tsv')
    expected_mouse_keys=mouse[['capture_id','disc_id','region']].drop_duplicates()
    if mouse.duplicated(['capture_id','disc_id','region','human_gene']).any():
        raise ValueError('Repeated mouse candidate disc-region identity')
    for gene in cfg['candidates']:
        if len(mouse.loc[mouse.human_gene.eq(gene)])!=len(expected_mouse_keys):
            raise ValueError('Mouse candidate coverage rows incomplete: '+gene)
    b=load_new_b(cfg,ready,completion)
    matrix=build_matrix(cfg,frozen,tu,mouse,depth,b,ready)
    out=cfg['integrated_out'];out.mkdir(parents=True,exist_ok=True)
    write(matrix,out/'MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv')
    write(tu,out/'TU_CANDIDATE_DONOR_CONTEXT.tsv')
    write(mouse.loc[mouse.human_gene.isin(cfg['candidates'])],out/'SPATIAL_CANDIDATE_DISC_REGION_SOURCE.tsv')
    if ready:
        write(b['full'].loc[b['full'].TG.isin(cfg['candidates'])],out/'CANDIDATE_ALL_TESTED_NEW_AXES_WITH_PRS.tsv.gz')
        write(b['context'].loc[b['context'].TG.isin(cfg['candidates'])],out/'CANDIDATE_SOURCE_SUPPORTED_RECEPTOR_KEY_PATHS.tsv.gz')
        if len(b['prepare_lr']): write(b['prepare_lr'],out/'OFFICIAL_PREPARE_LR_GATE_SOURCE.tsv')
    by=matrix.set_index('gene')
    context=b.get('context',pd.DataFrame())
    fnctx=context.loc[context.ligand.eq('FN1')&truth(context.complex_contains_ITGB1)] if ready and len(context) else pd.DataFrame()
    passing_fn1_keys=set(by.loc['FN1','official_LR_gate_passing_receptor_keys'].split(';'))-{''}
    plod_supported=bool(len(fnctx.loc[fnctx.TG.eq('PLOD2')&fnctx.Receptor.isin(passing_fn1_keys)])) if len(fnctx) else False
    if not ready:
        decision='PENDING_NEW_FORMAL_RUN; do not finalize the FN1 exposure × ITGB1 knockdown nomination from historical or smoke results.'
    elif plod_supported:
        decision=('RETAIN_AS_TESTABLE_HYPOTHESIS. New source-supported FN1/integrin-complex receptor-key context includes PLOD2. '
          'A factorial FN1 exposure × ITGB1 knockdown experiment with PLOD2 as a prespecified response readout remains reasonable; '
          'complex co-members, FN1 non-exclusivity and depth-sensitive spatial support must remain explicit. This is not a predicted experimental outcome.')
    else:
        alternatives=matrix.loc[matrix.gene.isin(['CCN2','PRSS23','THBS1','LOX','TGFBI'])&
                                 matrix.n_AF_LR_receptor_key_context_rows.gt(0),'gene'].tolist()
        decision=('DO_NOT_PRIORITIZE_AS_AN_ESTABLISHED_FN1_ITGB1_PLOD2_AXIS. The new formal source-supported FN1/ITGB1-containing '
          'complex context does not supply a passing PLOD2 receiver path. PLOD2 can remain an exploratory biological readout, '
          'but should not be presented as computationally validated downstream of FN1. '
          +('Alternative passing receiver target contexts among the prespecified comparators are '+', '.join(alternatives)+
            '; inspect their exact ligand/complex/receptor/SSC/target rows and molecular feasibility before choosing one. ' if alternatives else
            'None of the prespecified alternative targets supplies a passing source-supported receiver context; no replacement mechanism is invented. ')
          +'Expression and compartment context alone do not rescue a missing formal path.')
    text(out/'FN1_INTEGRIN_CONTEXT_REAL_EVIDENCE.md','# FN1–integrin context: real evidence\n\n'+
      f"New formal completion verified: {ready}. FN1 ligand LR rows: {by.loc['FN1','AF_ligand_LR_rows']}; "
      f"CellChat p<0.05 rows: {by.loc['FN1','AF_ligand_LR_p005_rows']}. ITGB1-containing FN1 complexes: "
      f"{by.loc['FN1','FN1_ITGB1_complexes'] or 'not established by completed inputs'}. "
      f"Passing FN1/ITGB1-containing receptor-key context rows across targets: {len(fnctx)}.\n\n"+
      f"FN1 LR receiver-feature/relative-CCC gate pairs: {by.loc['FN1','official_receiver_feature_CCC_gate_pairs']}; "
      f"status {by.loc['FN1','official_receiver_feature_CCC_gate_reason']}. These are the original official LR-stage gates, not proof of a ligand-specific receiver path.\n\n"+
      'FN1 is a candidate sender ligand; ITGB1 is a receptor-complex component, not a single-molecule FN1 receptor claim. '
      'CellChat retains the original complex, while HGNN/PRS describe receiver Receptor–SSC–target predictions. '
      'Their join is by receptor key and does not show ligand-specific induction. Original row keys and official values '
      'are preserved in source tables. FN1 is detected in NP as well as AF; no AF-exclusive source or directional molecular transfer is established.\n\n'+decision)
    text(out/'FINAL_CANDIDATE_AND_WETLAB_DECISION.md','# Candidate and future-experiment decision\n\n'+decision+
      '\n\nNo joint LR×PA×PRS score was constructed. Clinical-grade evidence remains unavailable. Mouse spatial support is weakened by depth sensitivity. '
      'All prespecified targets, including unavailable/untested/PA-failing/PRS-missing targets, remain in the evidence matrix. '
      'No qPCR, western blot, ELISA or perturbation outcome was simulated or reused as a new result.')
    text(out/'FUTURE_EXPERIMENT_PLAN.md','# Future experiment plan: not executed\n\n'+decision+
      '\n\nIf pursued after reviewing the complete formal evidence, use a prespecified factorial design crossing FN1 exposure with '
      'ITGB1 knockdown and appropriate vehicle/non-targeting controls. Confirm knockdown, viability, receptor-complex context and '
      'exposure conditions before interpreting an NP response. Biological replicates, reagent feasibility, timing, doses and power '
      'require independent experimental specification. PLOD2 is an a priori response candidate; retain alternative readouts and '
      'null/opposite responses. No expected qPCR/WB/ELISA values are supplied.')
    claims=[
      dict(figure='Figure1',claim='frozen programs and disclosed cohort/identity coverage',source='01_human_programs/COHORT_AND_DONOR_MANIFEST.tsv;frozen program audit',status='EXPRESSION_PROVENANCE',boundary='no cross-series disease effect'),
      dict(figure='Figure2',claim='independent within-donor fixed-program expression context',source='01_human_programs/PROGRAM_SCORE_DONOR_SUMMARY.tsv;ALL_PROGRAMS_BRIEF.tsv',status='DONOR_DESCRIPTIVE_EXPRESSION',boundary='no Tu clinical-grade test; AUCell/UCell scales differ'),
      dict(figure='Figure3',claim='new AF LR and receiver Receptor–SSC–target predictions',source='02_sigxtalk_official_new/formal_new/tables/ALL_TESTED_NEW_AXES_WITH_PRS.tsv;FOCUS_LR_COMPLEX_ALL.tsv',status='NEW_FORMAL_PREDICTION' if ready else 'PENDING_NO_NEW_FORMAL_CLAIM',boundary='complex-component receptor-key join;not ligand-specific causality'),
      dict(figure='Figure4',claim='original adult mouse context is depth-sensitive',source='03_spatial_gse287499/DEPTH_DISC_PAIR_COMPARISON.tsv;CANDIDATE_GENES_PER_DISC_REGION.tsv',status=depth['outcome'],boundary='disc nested in capture;genotype confounded;provisional anatomy;no human spatial validation'),
      dict(figure='Figure5',claim='future perturbation only',source='04_integrated/FUTURE_EXPERIMENT_PLAN.md',status='NOT_EXECUTED_NO_FIGURE',boundary='old THBS1 experiments and synthetic fixtures are not new validation')]
    write(pd.DataFrame(claims),out/'FIGURE_CLAIM_EVIDENCE_MAP.tsv')
    release_ok=manuscript(cfg,matrix,ready,reasons,completion,final,hgnn,provenance,b,depth)
    inputs=[cfg['tu_candidate_table'],cfg['frozen_programs'],cfg['human_dir']/'PROGRAM_SCORE_DONOR_SUMMARY.tsv',
            cfg['human_dir']/'GEP2_DONOR_SUPPORT.tsv',cfg['human_dir']/'STAGE_A_COMPUTATION_STATUS.json',
            cfg['human_dir']/'ALL_PROGRAMS_BRIEF.tsv',cfg['spatial_depth_dir']/'DEPTH_ANALYSIS_STATUS.json',
            cfg['spatial_depth_dir']/'CANDIDATE_GENES_PER_DISC_REGION.tsv']
    if ready: inputs.extend([cfg['sigxtalk_run']/'audit/NEW_R_PYTHON_R_COMPLETION.json',
        cfg['sigxtalk_run']/'NEW_OFFICIAL_FINAL_STAGE_STATUS.json',cfg['sigxtalk_run']/'tables/ALL_TESTED_NEW_AXES_WITH_PRS.tsv',
        cfg['sigxtalk_run']/'audit/NEW_INPUT_PROVENANCE.json',cfg['sigxtalk_run']/'audit/HGNN_STATUS.json',
        cfg['sigxtalk_run']/'tables/FOCUS_LR_COMPLEX_ALL.tsv',cfg['sigxtalk_run']/'tables/LR_COMPLEX_COMPONENTS.tsv'])
    optional=[cfg['spatial_input']/'NP_PROGRAM_BY_REGION.tsv',
              cfg['spatial_depth_dir']/'GEP2_SCORE_UMI_ASSOCIATION.tsv',
              cfg['spatial_depth_dir']/'DEPTH_RETENTION_BY_CAPTURE_REGION.tsv']
    if ready:
        optional.extend([cfg['sigxtalk_run']/name for name in [
            'tables/TARGET_GENE_SELECTION_ALL.tsv','tables/CANDIDATE_TESTABILITY.tsv','tables/PREPARE_LR_RECEIVER_FILTER_AUDIT.tsv',
            'prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv','audit/R_SESSION.txt','audit/PRS_RESOURCE_PATCH_AUDIT.json']])
    optional.extend(cfg['figure_dir'].rglob('*LEGEND*.md'))
    for key in ['sigxtalk_config','github_proof']:
        if cfg.get(key): optional.append(cfg[key])
    if args.config: optional.append(args.config.resolve())
    inputs=list(dict.fromkeys(inputs+[p for p in optional if p.is_file()]))
    write(pd.DataFrame([dict(input_name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in inputs]),out/'INTEGRATION_INPUT_MANIFEST.tsv')
    status=dict(status='DONE_NEW_FORMAL_INTEGRATION' if ready else 'PARTIAL_NEW_B_PENDING',
        finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(),candidate_rows=len(matrix),formal_new_B_verified=ready,
        B_pending_reasons=reasons,B_run_class=final.get('run_class','NOT_COMPLETE'),
        B_epochs_completed=hgnn.get('epochs_completed'),old_SigXTalk_rows_used=False,
        depth_outcome=depth['outcome'],depth_attenuation=depth.get('observed_mean_difference_attenuated'),
        FN1_ITGB1_complex_context_with_PLOD2=plod_supported,no_composite_score=True,
        no_genotype_or_clinical_grade_P_values=True,wetlab_executed=False,
        github_code_release='VERIFIED_PUBLIC' if release_ok else 'BLOCKED_UNVERIFIED',
        source_script_sha256=sha(Path(__file__)),
        configuration_sha256=sha(args.config) if args.config else None,
        input_manifest_sha256=sha(out/'INTEGRATION_INPUT_MANIFEST.tsv'))
    text(out/'INTEGRATION_STATUS.json',json.dumps(status,indent=2));print(json.dumps(status,indent=2))

if __name__=='__main__': main()
