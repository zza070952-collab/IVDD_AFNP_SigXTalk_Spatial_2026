#!/usr/bin/env python3
"""Figure3 exclusively from completed new official IVDD outputs, never legacy/demo."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches


def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')


def write(d,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(p,sep='\t',index=False,na_rep='NA')


def title(ax,letter,text,y=1.08):
    ax.text(-.02,y,letter,transform=ax.transAxes,fontsize=11,fontweight='bold',ha='right')
    ax.text(.025,y,text,transform=ax.transAxes,fontsize=9.3,fontweight='bold',ha='left')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--out',type=Path)
    a=ap.parse_args();cp=a.config.resolve();cfg=json.loads(cp.read_text())
    root=Path(cfg['project_root']);root=root if root.is_absolute() else (cp.parent/root).resolve()
    module=Path(cfg['module_dir']);module=module if module.is_absolute() else root/module
    wd=module/'formal_new';tables=wd/'tables'
    status=json.loads((wd/'NEW_OFFICIAL_FINAL_STAGE_STATUS.json').read_text())
    provenance=json.loads((wd/'audit/NEW_INPUT_PROVENANCE.json').read_text())
    report=json.loads((tables/'NEW_REPORT_STATUS.json').read_text())
    if not (status['formal_design'] and status['new_result_not_legacy_relabel'] and
            status['complete_R_python_R'] and report['mode']=='formal' and
            provenance['input_kind']=='real_GSE230808_raw_counts' and
            not provenance['started_from_old_pathway_output']):
        raise ValueError('Figure3 requires new real-IVDD formal outputs, not PBMC/import/legacy results')
    out=(a.out or root/'05_figures').resolve();out.mkdir(parents=True,exist_ok=True);sd=out/'source_data/Figure3_new'
    inp=read(tables/'SENDER_RECEIVER_INPUT.tsv');lr=read(tables/'AF_TO_NP_CANDIDATE_LR.tsv')
    bridge=read(tables/'PROGRAM_COMMUNICATION_BRIDGE.tsv')
    reps=read(tables/'NEW_REPRESENTATIVE_PATHS.tsv')
    metfile=tables/'NEW_REPRESENTATIVE_OFFICIAL_FID_SPE.tsv'
    if metfile.is_file():reps=read(metfile)
    counts=inp.pivot(index='donor_id',columns='group',values='Freq').fillna(0).sort_index()
    gp2=bridge[bridge.program=='GEP2'].copy();assert len(gp2)==25
    top=lr[lr.pval<.05].sort_values(['prob','interaction_name'],ascending=[False,True]).head(8)
    write(counts.reset_index(),sd/'Figure3a_new_donor_input.tsv')
    write(top,sd/'Figure3b_new_LR_context.tsv')
    write(reps,sd/'Figure3c_new_path_metrics.tsv')
    write(gp2,sd/'Figure3d_new_GEP2_coverage.tsv')
    colors={'AF':'#528e83','high':'#8d78aa','low':'#d0d5d8','ink':'#25323c','muted':'#59636b','line':'#a6afb5'}
    matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
      'font.size':8,'axes.labelsize':8,'xtick.labelsize':7.4,'ytick.labelsize':7.4,
      'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
      'text.color':colors['ink'],'axes.labelcolor':colors['ink'],'pdf.fonttype':42,'svg.fonttype':'none',
      'figure.facecolor':'white','savefig.facecolor':'white','legend.frameon':False})
    fig=plt.figure(figsize=(8,7.5))
    fig.text(.075,.965,'New official SigXTalk: AF-source context and NP program hypotheses',fontsize=11.7,fontweight='bold')
    fig.text(.075,.935,f'Raw-count R → HGNN → PRS | one fixed pooled model | {len(counts)} paired GSE230808 patients',fontsize=8,color=colors['muted'])
    ax=fig.add_axes([.105,.645,.285,.225]);title(ax,'a','Actual new model input')
    y=np.arange(len(counts));left=np.zeros(len(counts))
    for group,label,key in [('AF_stromal','AF sender','AF'),('NP_GEP2_high','NP-high','high'),('NP_GEP2_low','NP-low','low')]:
        v=counts[group].to_numpy();ax.barh(y,v,left=left,color=colors[key],height=.67,edgecolor='white',linewidth=.5,label=label)
        for yy,x,l in zip(y,v,left):
            if x>0:ax.text(l+x/2,yy,str(int(x)),ha='center',va='center',fontsize=7,color='white' if key!='low' else colors['ink'])
        left+=v
    ax.set_yticks(y,counts.index);ax.invert_yaxis();ax.set_xlim(0,max(left)*1.05)
    ax.set_xlabel('Cells nested within patients');ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
    ax.legend(loc='lower left',bbox_to_anchor=(-.02,1.005),ncol=3,fontsize=7,handlelength=1.1,handletextpad=.3,columnspacing=.7,borderaxespad=0)
    ax=fig.add_axes([.665,.645,.285,.225]);title(ax,'b','New AF→NP LR context')
    if len(top):
        ax.barh(np.arange(len(top)),top.prob,color=colors['AF'],height=.68)
        ax.set_yticks(np.arange(len(top)),[f'{r.ligand} – {str(r.receptor).replace("_","+")}' for r in top.itertuples()])
        ax.invert_yaxis();ax.set_xlim(0,max(top.prob)*1.1)
    else:ax.text(.5,.5,'No original LR p<0.05 rows',ha='center',transform=ax.transAxes)
    ax.set_xlabel('Original CellChat probability');ax.tick_params(axis='y',length=0,labelsize=7)
    ax.spines['left'].set_visible(False)
    ax.text(.99,-.29,f'Top {len(top)} of {int((lr.pval<.05).sum())} algorithmic p<0.05 rows',ha='right',transform=ax.transAxes,fontsize=7.2,color=colors['muted'])
    ax=fig.add_axes([.08,.275,.87,.255]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    title(ax,'c','New source-linked receiver hypotheses')
    for x,label in [(.065,'AF ligand'),(.285,'Complex / HGNN node'),(.505,'NP SSC'),(.675,'GEP2 target'),(.89,'Original pathway values')]:
        ax.text(x,.96,label,ha='center',va='top',fontsize=7.5,color=colors['muted'])
    if not len(reps):
        ax.text(.5,.48,'No GEP2 target retained a joined LR and positive filtered PRS under the fixed rules.',ha='center',va='center',fontsize=9)
    for yy,row in zip([.73,.46,.19],reps.head(3).itertuples()):
        entries=[(.065,.105,row.ligand,'#e6f0ed'),(.285,.215,str(row.receptor).replace('_','+')+f'\nHGNN: {row.Receptor}','#eee9f4'),
                 (.505,.11,row.TF,'#e7eaf0'),(.675,.13,row.TG,'#eee9f4')]
        for x,w,text,color in entries:
            ax.add_patch(patches.FancyBboxPatch((x-w/2,yy-.085),w,.17,boxstyle='round,pad=0.009,rounding_size=.014',facecolor=color,edgecolor=colors['line'],linewidth=.7))
            ax.text(x,yy,text,ha='center',va='center',fontsize=7.8,linespacing=1.3)
        for x1,x2 in [(.13,.166),(.403,.44),(.571,.6)]:
            ax.annotate('',xy=(x2,yy),xytext=(x1,yy),arrowprops={'arrowstyle':'->','linestyle':(0,(3,2)),'color':colors['line'],'lw':.9,'mutation_scale':9})
        text=f'PA {row.pred_label:.3f} | PRS {row.PRS:.3g}'
        if hasattr(row,'fidelity_target_all'):
            text+=f'\nFid {row.fidelity_target_all:.2g} | Spe {row.specificity_receptor_all:.2g}'
        else:text+='\nFid / Spe: not available'
        ax.text(.89,yy,text,ha='center',va='center',fontsize=7.5,linespacing=1.55)
    ax.text(0,-.035,'SSC denotes a database-listed intermediate, not a confirmed transcription factor.',fontsize=7.3,color=colors['muted'])
    ax.text(0,-.115,'PA/PRS and official Fid/Spe are distinct quantities; joins do not identify ligand-specific causal effects.',fontsize=7.1,color=colors['muted'])
    ax=fig.add_axes([.08,.101,.87,.056]);tested=int(gp2.n_tested_axes_as_target.gt(0).sum())
    title(ax,'d',f'Frozen GEP2 target coverage: {tested} of 25 tested',y=1.95)
    ax.set_xlim(-.5,24.5);ax.set_ylim(0,1);ax.set_yticks([]);ax.set_xticks(range(25),gp2.gene,rotation=60,ha='right',rotation_mode='anchor')
    ax.tick_params(axis='x',length=0,pad=3,labelsize=7.4)
    for spine in ax.spines.values():spine.set_visible(False)
    for i,row in enumerate(gp2.itertuples()):
        key='high' if row.n_tested_axes_as_target>0 else 'low'
        ax.add_patch(patches.Rectangle((i-.43,.14),.86,.61,facecolor=colors[key],edgecolor='white',linewidth=.4))
    ax.legend(handles=[patches.Patch(color=colors['high'],label='Tested target (not necessarily activated)'),
                       patches.Patch(color=colors['low'],label='Untested target')],loc='lower left',bbox_to_anchor=(0,1.05),ncol=2,fontsize=7.2,borderaxespad=0)
    fig.text(.08,.012,'New official predictions; no donor-level communication test, clinical-grade replication or measured AF→NP signaling.',fontsize=7,color=colors['muted'])
    base=out/'Figure3_New_Official_SigXTalk'
    for ext in ['png','pdf','svg']:fig.savefig(str(base)+'.'+ext,dpi=400)
    plt.close(fig)
    legend=f'''Fig. 3 | Newly computed official SigXTalk predictions from disease-source AF and NP raw counts. a, Actual donor/group contributions to one pooled GSE230808 model; cells are nested observations, not independent patients. b, Up to eight highest original-probability AF→NP LR complex rows among those with algorithmic CellChat p<0.05. Complex gene names are preserved, with '+' replacing '_' for display only. c, Up to three distinct frozen GEP2 targets linked through original receptor components to LR context, selected from actual new PA>0.75 and returned positive PRS rows passing the fixed official relative-max filter. The original PA and PRS are shown separately. Fid and Spe are values returned through the official fidelity and specificity matrix functions on the new filtered PRS universe: fidelity normalizes by all retained receptor–SSC paths ending at that target; specificity normalizes by all retained SSC–target paths beginning at that receptor. Their denominators are conditional on the declared universe and do not represent causal confidence. Database-listed SSCs are not assumed to be DNA-binding transcription factors. Dashed receptor-key-linked hypotheses do not establish a ligand-specific downstream model or measured signaling transfer. d, All25 fixed GEP2 genes; {tested} were actually tested targets. Tested coverage does not mean every gene passed PA/PRS criteria, and untested genes are not negative biological effects. The new formal run used seed20261017 and100epochs; only the separate technical real-data smoke used shortened training. Neither historical October17 pathways nor PBMC smoke results were substituted into this figure. No patient-level disease-effect test or new experiment is shown. Source Data: Figure3a_new_donor_input.tsv, Figure3b_new_LR_context.tsv, Figure3c_new_path_metrics.tsv and Figure3d_new_GEP2_coverage.tsv. Complete outputs and new-run provenance are retained separately.
'''
    smoke_kind='PBMC' if cfg['smoke']['input_kind']=='pbmc_rds' else 'IVDD'
    legend=legend.replace('only the separate technical real-data smoke used shortened training.',
        f"the separate real-data {smoke_kind} technical subset used {cfg['smoke']['epochs']} epochs and was not a full published PBMC reproduction.")
    (out/'Figure3_New_Official_LEGEND.md').write_text(legend,encoding='utf-8')
    (out/'Figure3_New_Official_QA.json').write_text(json.dumps({'rendered':True,'width_inches':8,'height_inches':7.5,'png_dpi':400,
      'new_formal_counts_only':True,'patients':len(counts),'cells':int(counts.sum().sum()),'display_paths':min(3,len(reps)),
      'GEP2_tested_targets':tested,'visual_review':'PENDING_ACTUAL_IMAGE_INSPECTION'},indent=2)+'\n')
    print(base)


if __name__=='__main__':main()
