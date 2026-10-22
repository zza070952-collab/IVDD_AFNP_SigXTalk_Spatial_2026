#!/usr/bin/env python3
"""Summarize only the newly computed official outputs, with explicit layer joins."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def read(p):
    return pd.read_csv(p,sep='\t',float_precision='round_trip')


def write(d,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(p,sep='\t',index=False,na_rep='NA')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--mode',choices=['smoke','formal'],default='formal')
    ap.add_argument('--plot-official-metrics',action='store_true')
    a=ap.parse_args();cp=a.config.resolve();cfg=json.loads(cp.read_text())
    root=Path(cfg['project_root']);root=root if root.is_absolute() else (cp.parent/root).resolve()
    resolve=lambda key: Path(cfg[key]) if Path(cfg[key]).is_absolute() else root/cfg[key]
    wd=resolve('module_dir')/('formal_new' if a.mode=='formal' else cfg['smoke'].get('directory_name','technical_smoke'))
    status=json.loads((wd/'audit/NEW_R_PYTHON_R_COMPLETION.json').read_text())
    if not status['complete_R_python_R']:raise ValueError('New official R/HGNN/PRS did not complete')
    tables=wd/'tables';full=read(tables/'ALL_TESTED_NEW_AXES_WITH_PRS.tsv')
    lr=read(tables/'FOCUS_LR_COMPLEX_ALL.tsv');comp=read(tables/'LR_COMPLEX_COMPONENTS.tsv')
    inp=read(tables/'DONOR_GROUP_INPUT_COUNTS.tsv')
    active=full[full.pred_label>cfg['hgnn']['activation_threshold']].copy()
    write(inp,tables/'SENDER_RECEIVER_INPUT.tsv')
    lr['method_origin']='NEW_OFFICIAL_SIGXTALK';lr['not_donor_level_disease_pvalue']=True
    write(lr,tables/'AF_TO_NP_CANDIDATE_LR.tsv')
    write(active,tables/'RTF_TARGET_RESULTS_IF_AVAILABLE.tsv')
    signatures=read(resolve('frozen_programs_tsv'))
    if 'program' not in signatures and 'signature' in signatures:
        signatures['program']=signatures.signature.str.extract(r'^(GEP[1-5])')[0]
        signatures=signatures[signatures.signature.str.endswith('_top25')]
    assert {'gene','program'}.issubset(signatures)
    signatures=signatures[['program','gene']].drop_duplicates()
    rows=[]
    for z in signatures.itertuples(index=False):
        q=full[full.TG==z.gene];qa=q[q.pred_label>.75]
        rows.append({'program':z.program,'gene':z.gene,'n_tested_axes_as_target':len(q),
                     'n_activated_axes_as_target':len(qa),'n_returned_positive_PRS_axes_as_target':int(qa.PRS.notna().sum()),
                     'max_original_PA_as_target':q.pred_label.max() if len(q) else np.nan,
                     'max_original_PRS_as_target':q.PRS.max() if len(q) else np.nan,
                     'target_status':'not_tested' if not len(q) else 'tested',
                     'method_origin':'NEW_OFFICIAL_SIGXTALK',
                     'independent_GEP2_validation':False})
    bridge=pd.DataFrame(rows);write(bridge,tables/'PROGRAM_COMMUNICATION_BRIDGE.tsv')
    passing=active.PRS_relative_max_filter_pass.astype(str).str.lower().isin(['true','1'])
    eligible=active[passing & active.PRS.notna() & active.TG.isin(signatures[signatures.program=='GEP2'].gene)]
    comp=comp.rename(columns={'receptor_component':'Receptor','receptor_complex':'receptor'})
    keys=['source','target','ligand','receptor','interaction_name']
    expanded=lr.merge(comp,on=keys,how='inner',validate='one_to_many')
    joined=expanded[expanded.pval<.05].merge(eligible,on='Receptor',how='inner',validate='many_to_many')
    joined['receptor_key_join_not_ligand_specific_model']=True
    joined['SSC_not_verified_transcription_factor']=True
    joined['no_composite_score_created']=True
    joined=joined.sort_values(['PRS','pred_label','prob','TG','TF','ligand','interaction_name'],
                              ascending=[False,False,False,True,True,True,True],kind='mergesort')
    write(joined,tables/'GEP2_NEW_SOURCE_SUPPORTED_PATHS.tsv.gz')
    pool=joined.drop_duplicates('TG');reps=pool.head(3).copy()
    write(pool,tables/'NEW_REPRESENTATIVE_TARGET_POOL.tsv')
    write(reps,tables/'NEW_REPRESENTATIVE_PATHS.tsv')
    metrics_path=wd/'prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv'
    metric_rows=pd.DataFrame()
    if len(reps) and metrics_path.exists():
        metrics=read(metrics_path).rename(columns={'SSC':'TF','Target':'TG'})
        metric_rows=reps.merge(metrics,on=['Receptor','TF','TG'],validate='one_to_one',suffixes=('','_official_metric'))
        write(metric_rows,tables/'NEW_REPRESENTATIVE_OFFICIAL_FID_SPE.tsv')
    if a.plot_official_metrics and len(metric_rows):
        # Own visualization of returned official numbers; no R plotting or recreated metric.
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
                                     'font.size':8,'pdf.fonttype':42,'svg.fonttype':'none',
                                     'axes.spines.top':False,'axes.spines.right':False})
        fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),sharey=True)
        y=np.arange(len(metric_rows));labels=[f'{r.Receptor} / {r.TF} / {r.TG}' for r in metric_rows.itertuples()]
        for ax,col,title,color in zip(axs,['fidelity_target_all','specificity_receptor_all'],
                                      ['Official fidelity: target universe','Official specificity: receptor universe'],
                                      ['#678eaa','#9781ad']):
            ax.barh(y,metric_rows[col],color=color,height=.58)
            ax.set_xlabel(title);ax.set_yticks(y,labels);ax.invert_yaxis()
            ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2))
        fig.suptitle('New official pathway allocations within the declared filtered universe',fontsize=10)
        fig.text(.01,.01,'Database-listed intermediates (SSC); allocations are conditional predictions, not causal confidence.',fontsize=7)
        fig.subplots_adjust(left=.29,right=.98,wspace=.36,bottom=.28,top=.80)
        for ext in ['png','pdf','svg']:
            fig.savefig(wd/'tables'/('SUPP_NEW_OFFICIAL_FID_SPE.'+ext),dpi=400)
        plt.close(fig)
    result={'method_origin':'NEW_OFFICIAL_SIGXTALK','mode':a.mode,
            'not_a_legacy_reanalysis':True,'all_tested_rows':len(full),'PA_gt075_rows':len(active),
            'positive_PRS_rows':int(full.PRS.notna().sum()),'LR_focus_rows':len(lr),
            'LR_p005_rows':int((lr.pval<.05).sum()),'GEP2_tested_target_genes':int((bridge[bridge.program=='GEP2'].n_tested_axes_as_target>0).sum()),
            'representative_paths':len(reps),'representative_targets':reps.TG.tolist(),
            'official_fid_spe_display_rows':len(metric_rows),
            'molecular_chain_experimentally_verified':False,
            'old_rows_used':False}
    (tables/'NEW_REPORT_STATUS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
