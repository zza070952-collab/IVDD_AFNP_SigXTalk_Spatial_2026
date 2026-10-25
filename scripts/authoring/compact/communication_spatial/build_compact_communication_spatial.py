#!/usr/bin/env python3
"""Build compact supplement artifacts from private frozen caches.

No training, resampling, ROI changes, experimental source data, or manuscript
authoring. All numeric study results are read or derived at execution.
"""
import json,hashlib,argparse,os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,ListedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.ticker import FuncFormatter

def env_path(name):
    value = os.environ.get(name)
    return Path(value).expanduser() if value else None

def separate_output(output, *inputs):
    output = output.expanduser().resolve()
    for source in inputs:
        source = source.expanduser().resolve()
        if output == source or output in source.parents or source in output.parents:
            raise ValueError("Output must be separate from every private input cache: " + str(source))
    return output

A=argparse.ArgumentParser(description=__doc__)
A.add_argument('--v4', type=Path, default=env_path('IVDD_V4_ROOT'))
A.add_argument('--v3-review', type=Path, default=env_path('IVDD_V3_ROOT'))
A.add_argument('--roi-json', type=Path, default=env_path('IVDD_ROI_JSON'))
A.add_argument('--output-root', type=Path, default=env_path('IVDD_OUTPUT_ROOT'))
A.add_argument('--main-tables-root', type=Path, default=env_path('IVDD_MAIN_TABLES_ROOT'), help='Main communication source-table directory; defaults to the private V4 cache')
a=A.parse_args()
if not a.v4 or not a.v3_review or not a.output_root:
    A.error('Set --v4/IVDD_V4_ROOT, --v3-review/IVDD_V3_ROOT and --output-root/IVDD_OUTPUT_ROOT')
P=a.v4.expanduser().resolve(); V=a.v3_review.expanduser().resolve()
a.roi_json=(a.roi_json or V/'07_public_repo/configs/frozen_spatial_rois.json').expanduser().resolve()
T=(a.main_tables_root or P/'03_tables/communication_spatial').expanduser().resolve()
O=separate_output(a.output_root, P, V, a.roi_json, T); O.mkdir(parents=True,exist_ok=True)
Q=O/'source_tables';Q.mkdir(exist_ok=True);F=O/'figures';F.mkdir(exist_ok=True)
ST=V/'02_sigxtalk_official_new/formal_new/tables';SP=V/'03_spatial_gse287499'
SOURCES={};EXPORTS=[]
def read(path,id,portable):
    SOURCES[id]={'relative_path':portable,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}
    return pd.read_csv(path,sep='\t')
def write(d,name):d.to_csv(Q/name,sep='\t',index=False,na_rep='NA')
def k(row):return str(row['Receptor'])+'::'+str(row['TF'])+'::'+str(row['TG'])
NP='#268C92';AF='#C97B3A';GREY='#A8ADB3';C={'CTR':'#426D9A','cKO':'#A45A88'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False})
def label(ax,ch,title):
    ax.set_title(title,loc='left',pad=9);ax.text(-.15,1.07,ch,transform=ax.transAxes,fontsize=13,fontweight='bold',va='bottom')
def clean(ax):
    for s in ax.spines.values():s.set_visible(False)
    ax.tick_params(length=0)
def export(fig,name):
    fig.canvas.draw();ren=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height();skip=set();bad=[]
    for ax in fig.axes:
        for axis,lim in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            lo,hi=sorted(lim)
            for tick in list(axis.get_major_ticks())+list(axis.get_minor_ticks()):
                if not ax.axison or not lo-1e-10<=tick.get_loc()<=hi+1e-10:skip.update([tick.label1,tick.label2])
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text() or t in skip:continue
        b=t.get_window_extent(ren)
        if b.x0<0 or b.y0<0 or b.x1>w+1 or b.y1>h+1:bad.append(t.get_text())
    for ext in ['pdf','svg','png']:fig.savefig(F/(name+'.'+ext),dpi=300)
    EXPORTS.append({'figure':name,'width_mm':fig.get_figwidth()*25.4,'height_mm':fig.get_figheight()*25.4,'minimum_font_pt':min(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()),'text_outside_canvas':bad})
    plt.close(fig)

axes=read(ST/'ALL_TESTED_NEW_AXES_WITH_PRS.tsv','AXES','V3/02_sigxtalk_official_new/formal_new/tables/ALL_TESTED_NEW_AXES_WITH_PRS.tsv')
axes['axis_key']=axes.apply(k,axis=1)
assert axes.axis_key.is_unique and not axes.source_run_started_from_legacy_predictions.any()
lr=read(ST/'AF_TO_NP_CANDIDATE_LR.tsv','LR','V3/02_sigxtalk_official_new/formal_new/tables/AF_TO_NP_CANDIDATE_LR.tsv')
alllr=read(ST/'CELLCHAT_ALL_SOURCE_TARGET_LR.tsv','ALL_LR','V3/02_sigxtalk_official_new/formal_new/tables/CELLCHAT_ALL_SOURCE_TARGET_LR.tsv')
comp=read(ST/'LR_COMPLEX_COMPONENTS.tsv','LR_COMPONENTS','V3/02_sigxtalk_official_new/formal_new/tables/LR_COMPLEX_COMPONENTS.tsv')
raw=read(V/'02_sigxtalk_official_new/formal_new/prs_outputs/FID_SPE_RAW_UNIVERSE_OFFICIAL_NEW.tsv','FID_RAW','V3/02_sigxtalk_official_new/formal_new/prs_outputs/FID_SPE_RAW_UNIVERSE_OFFICIAL_NEW.tsv')
fm=read(V/'02_sigxtalk_official_new/formal_new/prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv','FID_FILTER','V3/02_sigxtalk_official_new/formal_new/prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv')
don=read(ST/'DONOR_GROUP_INPUT_COUNTS.tsv','DONORS','V3/02_sigxtalk_official_new/formal_new/tables/DONOR_GROUP_INPUT_COUNTS.tsv')
frozen=read(V/'01_human_programs/FROZEN_K5_PROGRAMS.tsv','FROZEN','V3/01_human_programs/FROZEN_K5_PROGRAMS.tsv')
targets=read(ST/'TARGET_GENE_SELECTION_ALL.tsv','TARGETS','V3/02_sigxtalk_official_new/formal_new/tables/TARGET_GENE_SELECTION_ALL.tsv')
proofpath=V/'02_sigxtalk_official_new/SIGXTALK_EXECUTION_PROOF.json'
SOURCES['PROOF']={'relative_path':'V3/02_sigxtalk_official_new/SIGXTALK_EXECUTION_PROOF.json','sha256':hashlib.sha256(proofpath.read_bytes()).hexdigest(),'bytes':proofpath.stat().st_size}
proof=json.loads(proofpath.read_text())
assert int(targets.selected_target.sum())==proof['counts']['selected_targets']
assert axes.TG.nunique()==proof['counts']['prepared_target_genes']
candidate=read(V/'04_integrated/MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv','CROSS','V3/04_integrated/MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv')
mainpaths=read(T/'Figure3_D_representative_paths.tsv','MAIN3D','V4/03_tables/communication_spatial/Figure3_D_representative_paths.tsv')
mainmetrics=read(T/'Figure3_E_official_fidelity_specificity.tsv','MAIN3E','V4/03_tables/communication_spatial/Figure3_E_official_fidelity_specificity.tsv')
matrix=read(T/'Figure3_C_receptor_target_max_PRS.tsv','MAIN3C','V4/03_tables/communication_spatial/Figure3_C_receptor_target_max_PRS.tsv')
genes7=['FN1','ITGB1','PLOD2','CCN2','PRSS23','LOX','THBS1']
active=axes[axes.PA_threshold_pass];pos=axes[axes.positive_official_PRS_returned];fixed=axes[axes.PRS_relative_max_filter_pass]
assert not axes.empty and set(fixed.index).issubset(pos.index) and set(pos.index).issubset(active.index)
for d in [raw,fm]:d['axis_key']=d.Receptor+'::'+d.SSC+'::'+d.Target
rawix=raw.set_index('axis_key');fixix=fm.set_index('axis_key');axix=axes.set_index('axis_key',drop=False)
selected={}
def add(q,reason):
    for row in q.to_dict('records'):
        key=k(row);selected.setdefault(key,[])
        if reason not in selected[key]:selected[key].append(reason)
add(mainpaths,'main_Figure3D_all12_fixed_keys');add(mainmetrics,'main_Figure3E_all3_metric_keys')
for z in matrix.itertuples():
    q=pos[(pos.Receptor==z.Receptor)&(pos.TG==z.Target)].sort_values(['PRS','Receptor','TF','TG'],ascending=[False,True,True,True])
    add(q.head(1),'main_Figure3C_max_positive_PRS_per_displayed_pair')
for gene in genes7:
    add(fixed[fixed.TG==gene].sort_values(['PRS','Receptor','TF','TG'],ascending=[False,True,True,True]).head(3),'top3_fixed_filter_PRS_per_discussed_target')
    add(active[(active.TG==gene)&active.PRS.isna()].sort_values(['pred_label','Receptor','TF','TG'],ascending=[False,True,True,True]).head(1),'highest_PA_active_NA_per_discussed_target')
add(fixed[fixed.Receptor.eq('ITGB1')].sort_values(['PRS','Receptor','TF','TG'],ascending=[False,True,True,True]).head(3),'top3_fixed_filter_PRS_ITGB1_receptor_role')
add(axes[(axes.Receptor=='ITGB1')&(axes.TF=='HIF1A')&(axes.TG=='PLOD2')],'explicit_main_text_predicted_axis_not_causal')
mandatory_count=len(selected)
for row in fixed.sort_values(['PRS','Receptor','TF','TG'],ascending=[False,True,True,True]).to_dict('records'):
    if len(selected)>=80:break
    if k(row) in selected:continue
    add(pd.DataFrame([row]),'fill_to80_by_fixed_filter_PRS_desc_then_receptor_SSC_target')
assert len(selected)==80 and mandatory_count<=80
rep=axix.loc[list(selected)].copy();rep['selection_rule']=[';'.join(selected[x]) for x in rep.index]
for metric,ix in [('raw',rawix),('filtered',fixix)]:
    rep['fidelity_'+metric]=rep.axis_key.map(ix.fidelity_target_all);rep['specificity_'+metric]=rep.axis_key.map(ix.specificity_receptor_all)
rep['NA_reason']=[('PRS_not_returned_for_PA_active_axis; metrics_not_defined; ' if pd.isna(r.PRS) else ('fixed_filter_metrics_NA_outside_filtered_universe; ' if not r.PRS_relative_max_filter_pass else ''))+'ligand_not_fitted_for_receiver_axis' for r in rep.itertuples()]
write(rep.reset_index(drop=True),'S3_80_REPRESENTATIVE_PATHS_NOT_ALL.tsv')
LR=lr.sort_values(['prob','interaction_name'],ascending=[False,True]).head(18)
main18keys=set(LR.interaction_name)
LR=pd.concat([LR,lr[lr.ligand.eq('FN1')]]).drop_duplicates('interaction_name').copy()
LR['selection_rule']=[';'.join((['main_Figure3B_top18_probability'] if r.interaction_name in main18keys else [])+(['all_observed_FN1_complexes'] if r.ligand=='FN1' else [])) for r in LR.itertuples()]
assert LR.selection_rule.str.contains('main_Figure3B_top18_probability').sum()==18
write(LR,'S3_MAIN18_PLUS_FN1_LR_COMPLEXES.tsv')
role=[]
for gene in genes7:
    q=axes[axes.TG==gene];r=axes[axes.Receptor==gene];pp=q[q.PA_threshold_pass]
    role.append({'gene':gene,'GEP2_member':bool(((frozen.program=='GEP2')&(frozen.gene==gene)).any()),'LR_ligand_records':int((lr.ligand==gene).sum()),'LR_component_records':int((comp.receptor_component==gene).sum()),'receptor_tested_axes':len(r),'target_tested_axes':len(q),'active_target_axes':len(pp),'target_positive_PRS':int(q.positive_official_PRS_returned.sum()),'target_active_PRS_NA':int((pp.PRS.isna()).sum()),'target_fixed_paths':int(q.PRS_relative_max_filter_pass.sum()),'target_status':'TESTED' if len(q) else 'NOT_SELECTED_TARGET; not a biological negative'})
role=pd.DataFrame(role);write(role,'S3_SEVEN_CANDIDATE_ROLES_AND_NA.tsv')
counts=pd.DataFrame([('HGNN tested',len(axes)),('PA > 0.75',len(active)),('Positive PRS',len(pos)),('Fixed filter',len(fixed)),('Active, PRS NA',int(active.PRS.isna().sum()))],columns=['state','n'])
write(counts,'S2A_OFFICIAL_COUNTS.tsv')
# Four actually computed metric examples, including the explicitly discussed path.
mkeys=[k(x) for x in mainmetrics.to_dict('records')]+['ITGB1::HIF1A::PLOD2'];m=rep.loc[mkeys].copy();write(m.reset_index(drop=True),'S2D_OFFICIAL_METRICS_RAW_AND_FILTERED.tsv')
fig=plt.figure(figsize=(180/25.4,180/25.4));aa=fig.add_axes([.19,.66,.29,.25]);bb=fig.add_axes([.66,.66,.28,.25]);cc=fig.add_axes([.15,.20,.30,.28]);dd=fig.add_axes([.72,.20,.22,.28])
aa.barh(range(5),counts.n,color=['#657B87','#426D9A',NP,AF,GREY],height=.6);aa.set_yticks(range(5));aa.set_yticklabels(counts.state);aa.invert_yaxis();aa.set_xscale('log');aa.set_xlim(3000,700000);aa.set_xticks([10000,100000]);aa.set_xticklabels(['10,000','100,000']);aa.set_xlabel('Axes (log scale)')
for i,n in enumerate(counts.n):aa.text(n*1.08,i,f'{n:,}',va='center',fontsize=8)
label(aa,'A','Distinct official output states')
for i,row in enumerate(role.itertuples()):
    if row.active_target_axes:
        frac=row.target_positive_PRS/row.active_target_axes;bb.barh(i,frac*100,color=NP,height=.6);bb.barh(i,(1-frac)*100,left=frac*100,color=GREY,height=.6);bb.text(102,i,f'{row.target_positive_PRS}/{row.target_active_PRS_NA}',va='center',fontsize=8)
    else:bb.text(50,i,'Target not tested',ha='center',va='center',fontsize=8)
bb.set(yticks=range(7),yticklabels=genes7,xlim=(0,155),xticks=[0,50,100],xlabel='Among PA-active target axes (%)');bb.invert_yaxis();label(bb,'B','Returned PRS and missing PRS')
bb.legend(handles=[Line2D([0],[0],color=NP,lw=5,label='Positive'),Line2D([0],[0],color=GREY,lw=5,label='NA')],loc='upper center',bbox_to_anchor=(.38,-.28),ncol=2,handlelength=1.2,columnspacing=1)
rolecols=['GEP2_member','LR_ligand_records','LR_component_records','target_tested_axes'];arr=np.column_stack([role[z].astype(bool).astype(int) for z in rolecols]);cc.imshow(arr,cmap=ListedColormap(['#ECEEF0','#5B96A0']),vmin=0,vmax=1,aspect='auto')
for j in range(7):
    for i in range(4):cc.text(i,j,'yes' if arr[j,i] else 'no',ha='center',va='center',fontsize=8,color='white' if arr[j,i] else '#555555')
cc.set(yticks=range(7),yticklabels=genes7,xticks=range(4));cc.set_xticklabels(['GEP2','LR ligand','LR component','Tested target'],rotation=40,ha='right');clean(cc);label(cc,'C','Separate candidate roles')
metriclabels=[]
for i,(key,row) in enumerate(m.iterrows()):
    metriclabels.append(row.Receptor+' / '+row.TF+'\n'+row.TG)
    for off,name,col,mark in [(-.11,'fidelity',NP,'o'),(.11,'specificity',AF,'^')]:
        xx=[row[name+'_raw'],row[name+'_filtered']];dd.plot(xx,[i+off]*2,color=col,lw=.7);dd.scatter(xx[0],i+off,facecolors='white',edgecolors=col,marker=mark,s=32,linewidths=.9);dd.scatter(xx[1],i+off,c=col,marker=mark,s=32,linewidths=.5)
dd.set(yticks=range(4),yticklabels=metriclabels,xlim=(-.004,.185),xlabel='Official fraction');dd.invert_yaxis();dd.set_xticks([0,.05,.10,.15]);label(dd,'D','Actual official metrics')
fig.legend(handles=[Line2D([0],[0],marker='o',color=NP,ls='',label='Fidelity'),Line2D([0],[0],marker='^',color=AF,ls='',label='Specificity'),Line2D([0],[0],marker='o',markerfacecolor='white',markeredgecolor='#555555',ls='',label='Raw-positive'),Line2D([0],[0],marker='o',color='#555555',ls='',label='Fixed-filter')],loc='lower center',bbox_to_anchor=(.74,.050),ncol=2,handletextpad=.4,columnspacing=.8)
ng=frozen.loc[frozen.program.eq('GEP2'),'gene']; nt=int(ng.isin(axes.TG).sum())
fig.text(.15,.030,f'{nt}/{len(ng)} GEP2 targets tested;\n{len(ng)-nt} not tested.\nITGB1:{int(axes.Receptor.eq("ITGB1").sum()):,} receptor axes;\nreceptor/target roles separate.',fontsize=8,va='bottom')
export(fig,'Supplementary_Figure_S2_Official_Prediction_Boundaries')

# Geometry-only compact S3: no H&E pixels or third-party figure raster.
roi=json.loads(a.roi_json.read_text());SOURCES['ROI']={'relative_path':'V3/07_public_repo/configs/frozen_spatial_rois.json','sha256':hashlib.sha256(a.roi_json.read_bytes()).hexdigest(),'bytes':a.roi_json.stat().st_size}
caps=sorted(roi['captures'],key=lambda x:x['genotype']!='CTR');discs=[(cap,ds) for cap in caps for ds in cap['discs']]
sp=read(P/'01_sources_readonly/spatial_reference/SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz','SPATIAL_BINS','V4/01_sources_readonly/spatial_reference/SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz')
dep=read(SP/'DEPTH_DISC_PAIR_COMPARISON.tsv','DEPTH','V3/03_spatial_gse287499/DEPTH_DISC_PAIR_COMPARISON.tsv')
ret=read(SP/'DISC_LEVEL_PROGRAM_AND_DEPTH.tsv','RETENTION','V3/03_spatial_gse287499/DISC_LEVEL_PROGRAM_AND_DEPTH.tsv')
gn=read(SP/'CANDIDATE_GENES_PER_DISC_REGION.tsv','SPATIAL_GENES','V3/03_spatial_gse287499/CANDIDATE_GENES_PER_DISC_REGION.tsv')
mk=read(P/'01_sources_readonly/spatial_reference/roi_validation/INDEPENDENT_MARKERS_BY_DISC_REGION.tsv','MARKERS','V4/01_sources_readonly/spatial_reference/roi_validation/INDEPENDENT_MARKERS_BY_DISC_REGION.tsv')
pg=read(P/'02_analysis/spatial_allprograms/NP_PROGRAM_BY_REGION.tsv','SPATIAL_PROGRAMS','V4/02_analysis/spatial_allprograms/NP_PROGRAM_BY_REGION.tsv')
elig=dep[dep.adequate_fixed_depth_pair].copy();before=elig.original_eligible_NP_minus_AF.mean();after=elig.depth100_NP_minus_AF.mean();atten=100*(1-after/before)
assert len(dep)==len(discs) and len(elig)>0 and np.isfinite(atten)
npositive=int((elig.depth100_NP_minus_AF>0).sum())
fig=plt.figure(figsize=(180/25.4,237/25.4));geometry=[]
for i,(cap,ds) in enumerate(discs):
    ax=fig.add_axes([.04+(i%4)*.24,.81-(i//4)*.16,.22,.12]);vv=np.vstack([r['vertices'] for r in ds['regions']]);x0,x1,y0,y1=vv[:,0].min()-45,vv[:,0].max()+45,vv[:,1].min()-45,vv[:,1].max()+45
    for reg in ds['regions']:
        ax.add_patch(Polygon(reg['vertices'],closed=True,fill=False,edgecolor=NP if reg['region']=='NP_CORE' else AF,lw=1))
        for j,(x,y) in enumerate(reg['vertices']):geometry.append({'capture':cap['capture_id'],'disc':ds['disc_id'],'region':reg['region'],'flank':reg['side'],'vertex':j,'x_image_pixel':x,'y_image_pixel':y,'microns_per_pixel':cap['microns_per_pixel']})
    d=sp[(sp.capture_id==cap['capture_id'])&(sp.disc_id==ds['disc_id'])&sp.anatomical_region.isin(['NP_CORE','AF_CORE'])]
    ax.scatter(d.x_fullres,d.y_fullres,c=[NP if r=='NP_CORE' else AF for r in d.anatomical_region],s=4,alpha=.6,linewidths=0)
    ax.set(xlim=(x0,x1),ylim=(y1,y0));ax.set_aspect('equal');ax.axis('off');ax.set_title(ds['disc_id'],fontsize=9,pad=3)
    yy=y1-.06*(y1-y0);xx=x0+.08*(x1-x0);length=200/cap['microns_per_pixel'];ax.plot([xx,xx+length],[yy,yy],c='#333333',lw=1.5);ax.text(xx+length/2,yy-.03*(y1-y0),'200 µm',ha='center',va='bottom',fontsize=8)
fig.text(.035,.965,'A',fontsize=13,fontweight='bold');fig.text(.09,.965,f'All {len(discs)} frozen core ROIs: original coordinate geometry',fontsize=9)
write(pd.DataFrame(geometry),'S3A_ALL8_ROI_COORDINATE_VERTICES.tsv')
b=fig.add_axes([.12,.405,.36,.16]);c=fig.add_axes([.70,.405,.25,.16]);markers=['Krt19','T','Krt8','Acan','Chad','Clec3a'];discorder=[ds['disc_id'] for cp,ds in discs];dr=[(di,r) for di in discorder for r in ['NP_CORE','AF_CORE']]
mx=mk.copy();mx['disc_region']=mx.disc_id+' '+mx.histology_core_region.str.replace('_CORE','',regex=False);mm=mx.pivot(index='marker',columns='disc_region',values='mean_log1p_per_10000').reindex(index=markers,columns=[di+' '+r.replace('_CORE','') for di,r in dr]);im=b.imshow(mm,aspect='auto',cmap='magma',vmin=0,vmax=mm.max().max());b.set(yticks=range(6),yticklabels=markers,xticks=range(16));b.set_xticklabels([x.replace('_',' ') for x in mm.columns],rotation=90);label(b,'B','Independent structural markers')
bc=fig.add_axes([.51,.405,.012,.16]);cb=fig.colorbar(im,cax=bc);cb.set_label('Mean log(1 + CP10k)')
gg=gn[gn.human_gene.isin(genes7)].groupby(['genotype','region','human_gene'],sort=False).agg(mean_expression=('mean_log1p_CP10k','mean'),mean_detection=('fraction_detected_bins','mean'),raw_UMI=('total_raw_UMI','sum')).reset_index();cols=[('CTR','NP_CORE'),('CTR','AF_CORE'),('cKO','NP_CORE'),('cKO','AF_CORE')]
for j,gene in enumerate(genes7):
    for i,(gt,r) in enumerate(cols):
        z=gg[(gg.genotype==gt)&(gg.region==r)&(gg.human_gene==gene)].iloc[0]
        if z.mean_detection==0:c.scatter(i,j,s=12,facecolors='white',edgecolors='#888888',lw=.7)
        else:c.scatter(i,j,s=10+140*z.mean_detection,c=[z.mean_expression],vmin=0,vmax=mm.max().max(),cmap='magma',edgecolors='white',lw=.3)
c.set(yticks=range(7),yticklabels=genes7,xticks=range(4),ylim=(6.6,-.6));c.set_xticklabels(['CTR\nNP','CTR\nAF','cKO\nNP','cKO\nAF']);clean(c);label(c,'C','Candidate detection coverage')
fig.text(.68,.34,'Area:detected fraction; hollow:0.\nCoverage from original cached counts.',fontsize=8)
d=fig.add_axes([.12,.105,.39,.18]);e=fig.add_axes([.70,.105,.25,.18])
for i,r in enumerate(elig.itertuples()):
    off=(i-3)*.016;d.plot([off,1+off],[r.original_eligible_NP_minus_AF,r.depth100_NP_minus_AF],c=C[r.genotype],marker='o' if r.genotype=='CTR' else '^',lw=.85,ms=4,alpha=.8)
d.plot([0,1],[before,after],'D-',c='#20252A',lw=2,ms=5);d.axhline(0,c='#888888',lw=.6,ls='--');d.set(xticks=[0,1],xticklabels=['Same eligible bins\noriginal counts','One fixed draw\n100 UMI/bin'],ylabel='GEP2 NP − AF contrast',xlim=(-.15,1.15),ylim=(-.06,.35));label(d,'D',f'{len(elig)} eligible structures: depth counterevidence')
d.text(.5,.97,f'{atten:.1f}% attenuation',ha='center',transform=d.transAxes,fontsize=10,fontweight='bold');d.text(.5,-.37,f'{before:.6f} → {after:.6f}; {npositive}/{len(elig)} remain NP > AF',ha='center',transform=d.transAxes,fontsize=8)
tab=[]
for di in discorder:
    r=dep[dep.disc_id==di].iloc[0];tab.append([di,str(r.NP_retained_bins),str(r.AF_retained_bins),'yes' if r.adequate_fixed_depth_pair else 'NO'])
e.axis('off');table=e.table(cellText=tab,colLabels=['Structure','NP n','AF n','Pass'],cellLoc='center',bbox=[-.03,0,1.08,1],colWidths=[.42,.22,.22,.22]);table.auto_set_font_size(False);table.set_fontsize(8)
for (i,j),cell in table.get_celld().items():cell.set_edgecolor('white');cell.set_facecolor('#E9EFF2' if i==0 else ('#F5F6F7' if i%2 else 'white'))
label(e,'E','Eligibility for every structure')
excluded=dep.loc[~dep.adequate_fixed_depth_pair,'disc_id'].astype(str).tolist()
fig.text(.70,.030,'≥5 bins/core required.\nExcluded: '+', '.join(excluded)+'\nNo animal/genotype P value.',fontsize=8,va='bottom')
export(fig,'Supplementary_Figure_S3_ROIs_Coverage_Depth')
write(gg,'S3C_SEVEN_GENE_DISPLAY_VALUES.tsv');write(dep,'S3DE_ALL8_DEPTH_ELIGIBILITY.tsv');write(mk,'S3B_ALL_MARKER_VALUES.tsv')

# Compact static workbook payloads; source records are portable and hash-specific.
H3=['record_type','record_key','ligand','receptor_complex','receptor_node','SSC','target','PA','PRS','fidelity_raw','specificity_raw','fidelity_filtered','specificity_filtered','LR_probability','LR_P','status','NA_reason','selection_rule','source_id','value_or_n']
R3=[]
def row3(**kw):R3.append({z:kw.get(z,None) for z in H3})
ng=frozen.loc[frozen.program.eq('GEP2'),'gene']; tested_gep=int(ng.isin(axes.TG).sum())
group_counts=don.groupby('group').Freq.sum()
summary=[('paired_donors',don.donor_id.nunique()),('input_cells',int(don.Freq.sum())),('AF_cells',int(group_counts.get('AF_stromal',0))),('NP_high_cells',int(group_counts.get('NP_GEP2_high',0))),('NP_background_cells',int(group_counts.get('NP_GEP2_low',0))),('AF_NP_LR',len(lr)),('LR_algorithm_P_lt005',int((lr.pval<.05).sum())),('tested_axes',len(axes)),('PA_gt075',len(active)),('positive_PRS_returned',len(pos)),('active_PRS_NA',int(active.PRS.isna().sum())),('fixed_filter_paths',len(fixed)),('target_selected',int(targets.selected_target.sum())),('target_prepared',axes.TG.nunique()),('GEP2_tested_genes',tested_gep),('GEP2_not_tested_genes',len(ng)-tested_gep),('epochs',proof['training']['epochs_completed']),('seed',proof['training']['seed']),('PA_threshold',.75),('relative_filter_not_P',.01)]
for key,n in summary:
    src='PROOF'
    if key in ['paired_donors','input_cells','AF_cells','NP_high_cells','NP_background_cells']:src='DONORS'
    elif key in ['AF_NP_LR','LR_algorithm_P_lt005']:src='LR'
    elif key in ['tested_axes','PA_gt075','positive_PRS_returned','active_PRS_NA','fixed_filter_paths','target_prepared']:src='AXES'
    elif key=='target_selected':src='TARGETS'
    elif key.startswith('GEP2_'):src='AXES;FROZEN'
    status=f'NEW_OFFICIAL_RUN; pooled {don.donor_id.nunique()} donors'
    if key=='GEP2_not_tested_genes':status+='; untested='+','.join(g for g in ng if g not in set(axes.TG))
    row3(record_type='RUN_SUMMARY',record_key=key,value_or_n=n,status=status,source_id=src)
row3(record_type='RUN_METADATA',record_key='upstream_source_commit',status=proof['source_commit'],NA_reason='not an effect estimate; official resource/RNG wrapper provenance read from private execution proof',source_id='PROOF')
for donor,g in don.groupby('donor_id'):
    row3(record_type='DONOR_INPUT',record_key=donor,status=';'.join(f'{r.group}={r.Freq}' for r in g.itertuples()),value_or_n=int(g.Freq.sum()),source_id='DONORS')
for target in ['NP_GEP2_high','NP_GEP2_low']:
    q=alllr[(alllr.source=='AF_stromal')&(alllr.target==target)];row3(record_type='MAIN_FIG3A_EDGE',record_key='AF_stromal→'+target,value_or_n=int((q.pval<.05).sum()),status='edge width=count algorithmP<.05; total_LR='+str(len(q)),source_id='ALL_LR')
for r in LR.itertuples():row3(record_type='LR_REPRESENTATIVE',record_key=r.interaction_name,ligand=r.ligand,receptor_complex=r.receptor,LR_probability=r.prob,LR_P=r.pval,status='AF_stromal→NP_GEP2_high; receptor complex retained',NA_reason='PA/PRS/official path metrics not LR quantities',selection_rule=r.selection_rule,source_id='LR')
for r in role.itertuples():row3(record_type='CANDIDATE_ROLE',record_key=r.gene,ligand=r.gene if r.LR_ligand_records else None,receptor_node=r.gene if r.receptor_tested_axes else None,target=r.gene if r.target_tested_axes else None,status=f'GEP2={r.GEP2_member}; ligand_LR={r.LR_ligand_records}; component_LR={r.LR_component_records}; receptor_axes={r.receptor_tested_axes}; target_axes={r.target_tested_axes}; active_targets={r.active_target_axes}; positivePRS_targets={r.target_positive_PRS}; NA_targets={r.target_active_PRS_NA}; filtered_targets={r.target_fixed_paths}',NA_reason=r.target_status+'; role summary is not a single path',value_or_n=r.target_active_PRS_NA,source_id='AXES;LR;LR_COMPONENTS;FROZEN')
for r in rep.itertuples():row3(record_type='REPRESENTATIVE_PATH_NOT_ALL',record_key=r.axis_key,receptor_node=r.Receptor,SSC=r.TF,target=r.TG,PA=r.pred_label,PRS=r.PRS,fidelity_raw=r.fidelity_raw,specificity_raw=r.specificity_raw,fidelity_filtered=r.fidelity_filtered,specificity_filtered=r.specificity_filtered,status=f'PApass={r.PA_threshold_pass}; positivePRS={r.positive_official_PRS_returned}; fixedfilter={r.PRS_relative_max_filter_pass}',NA_reason=r.NA_reason,selection_rule=r.selection_rule,source_id='AXES;FID_RAW;FID_FILTER')
for r in matrix.itertuples():row3(record_type='MAIN_FIG3C_MATRIX',record_key=r.Receptor+'::'+r.Target,receptor_node=r.Receptor,target=r.Target,PRS=r.max_returned_positive_PRS,status=r.display_state,NA_reason='maximum positive PRS over SSC; untested is not zero; not a single ligand fit',value_or_n=r.tested_axes,source_id='MAIN3C')

H4=['record_type','capture','disc','entity','NP_mean','AF_mean','NP_detected_fraction','AF_detected_fraction','NP_n','AF_n','NP_raw_UMI','AF_raw_UMI','eligible','paired_contrast','NA_or_boundary','source_id'];R4=[]
def row4(**kw):R4.append({z:kw.get(z,None) for z in H4})
for r in dep.itertuples():
    re=ret[ret.disc_id.eq(r.disc_id)].set_index('region');row4(record_type='ROI_QC',capture=r.capture_id,disc=r.disc_id,entity='NP core and bilateral AF cores',NP_mean=re.loc['NP_CORE','original_UMI_median'],AF_mean=re.loc['AF_CORE','original_UMI_median'],NP_detected_fraction=re.loc['NP_CORE','retention_fraction'],AF_detected_fraction=re.loc['AF_CORE','retention_fraction'],NP_n=int(re.loc['NP_CORE','original_n_bins']),AF_n=int(re.loc['AF_CORE','original_n_bins']),eligible='YES' if r.adequate_fixed_depth_pair else 'NO',NA_or_boundary='means=medianUMI/bin; fractions=retention>=100UMI; provisionalROI; animalID_unverified;1capture/genotype',source_id='RETENTION;ROI')
    for state,npv,afv,nnp,naf in [('GEP2_ORIGINAL',r.original_all_NP_CORE,r.original_all_AF_CORE,int(re.loc['NP_CORE','original_n_bins']),int(re.loc['AF_CORE','original_n_bins'])),('GEP2_SAME_ELIGIBLE_ORIGINAL',r.original_eligible_NP_CORE,r.original_eligible_AF_CORE,r.NP_retained_bins,r.AF_retained_bins),('GEP2_SINGLE100UMI',r.depth100_NP_CORE,r.depth100_AF_CORE,r.NP_retained_bins,r.AF_retained_bins)]:
        row4(record_type=state,capture=r.capture_id,disc=r.disc_id,entity='GEP2',NP_mean=npv,AF_mean=afv,NP_n=nnp,AF_n=naf,eligible='YES' if r.adequate_fixed_depth_pair else 'NO',paired_contrast=npv-afv,NA_or_boundary=f'meanlog1pCP10k; {r.disc_id} excluded from eligible paired mean' if not r.adequate_fixed_depth_pair else 'meanlog1pCP10k; sameeligiblepair; descriptive_notanimal',source_id='DEPTH')
for (capture,disc,gene),q in gn[gn.human_gene.isin(genes7)].groupby(['capture_id','disc_id','human_gene']):
    q=q.set_index('region');n=q.loc['NP_CORE'];f=q.loc['AF_CORE'];row4(record_type='CANDIDATE_ORTHOLOG',capture=capture,disc=disc,entity=gene,NP_mean=n.mean_log1p_CP10k,AF_mean=f.mean_log1p_CP10k,NP_detected_fraction=n.fraction_detected_bins,AF_detected_fraction=f.fraction_detected_bins,NP_n=int(n.n_detected_bins),AF_n=int(f.n_detected_bins),NP_raw_UMI=n.total_raw_UMI,AF_raw_UMI=f.total_raw_UMI,NA_or_boundary='mouse='+str(n.mouse_symbol)+'; feature_present='+str(n.native_feature_assessable)+'; zeros measured; meanlog1pCP10k; no genotype inference',source_id='SPATIAL_GENES')
for (capture,disc,marker),q in mk.groupby(['capture_id','disc_id','marker']):
    q=q.set_index('histology_core_region');n=q.loc['NP_CORE'];f=q.loc['AF_CORE'];row4(record_type='INDEPENDENT_MARKER',capture=capture,disc=disc,entity=marker,NP_mean=n.mean_log1p_per_10000,AF_mean=f.mean_log1p_per_10000,NP_detected_fraction=n.detected_bin_fraction,AF_detected_fraction=f.detected_bin_fraction,NP_n=int(n.n_detected_bins),AF_n=int(f.n_detected_bins),NP_raw_UMI=n.total_raw_UMI,AF_raw_UMI=f.total_raw_UMI,NA_or_boundary='mouse marker; meanlog1pCP10k; annotation compatibility, not pathology certification',source_id='MARKERS')
for (capture,disc),q in pg[pg.program.eq('GEP5')].groupby(['capture_id','disc_id']):
    q=q.set_index('region');n=q.loc['NP_CORE'];f=q.loc['AF_CORE'];row4(record_type='GEP5_DISCORDANT_PROGRAM',capture=capture,disc=disc,entity='GEP5',NP_mean=n.mean_program_score,AF_mean=f.mean_program_score,NP_n=int(n.n_bins),AF_n=int(f.n_bins),paired_contrast=n.mean_program_score-f.mean_program_score,NA_or_boundary='AF-favored context retained; fixedprogram; meanlog1pCP10k; no refit',source_id='SPATIAL_PROGRAMS')
for label0,val,unit in [('same_eligible_original_mean',before,'NP−AF meanlog1pCP10k'),('single100UMI_mean',after,'NP−AF meanlog1pCP10k'),('attenuation_percent',atten,'percent;100*(1-after/before)')]:row4(record_type='DEPTH_SUMMARY',entity=label0,NP_n=len(elig),AF_n=npositive,paired_contrast=val,eligible=f'{len(elig)}eligible;{npositive}positive',NA_or_boundary=unit+';same eligible structures; one cached fixed draw; no animal P',source_id='DEPTH')

H5=['record_type','gene','human_program_role','Tu_NP_detected','Tu_NP_gt_immune','LR_role','receiver_role','fixed_target_paths','mouse_NP_UMI','mouse_AF_UMI','mouse_NP_detected_discs','mouse_AF_detected_discs','experimental_role_and_release','negative_or_missing','source_id'];R5=[]
experimental={'FN1':'Plasma-FN coating exposure; not AF secretome; source unresolved','ITGB1':'siRNA perturbation and mRNA readout; protein not validated; source unresolved','PLOD2':'Primary qPCR endpoint; supplied interactions not established; source unresolved','CCN2':'Secondary qPCR transcript; source unresolved','PRSS23':'Secondary qPCR transcript; source unresolved','LOX':'NOT_ASSAYED in current cell experiment','THBS1':'NOT_ASSAYED in current cell experiment'}
# Counterevidence is read from the private integrated evidence table.
for gene in genes7:
    z=candidate.set_index('gene').loc[gene];r=role.set_index('gene').loc[gene];pr=frozen[frozen.gene.eq(gene)];R5.append(dict(zip(H5,['CANDIDATE',gene,';'.join(pr.program+' rank'+pr['rank'].astype(str)) if len(pr) else 'Not a frozen top25member',f'{z.Tu_NP_donors_detected}/{z.Tu_NP_donors}',f'{z.Tu_NP_gt_context_donors}/{z.Tu_NP_vs_context_assessable_donors}; comparatorNA={int(z.Tu_NP_donors-z.Tu_NP_vs_context_assessable_donors)}',f'ligandLR={r.LR_ligand_records}; receptorComponentLR={r.LR_component_records}',f'targetAxes={r.target_tested_axes}; receptorAxes={r.receptor_tested_axes}',int(z.n_filtered_PRS_axes) if r.target_tested_axes else None,z.mouse_NP_total_UMI,z.mouse_AF_total_UMI,f'{z.mouse_NP_discs_detected}/{z.mouse_NP_discs_assessable}',f'{z.mouse_AF_discs_detected}/{z.mouse_AF_discs_assessable}',experimental[gene],str(z.key_counterevidence),'CROSS;AXES;LR_COMPONENTS;SPATIAL_GENES;FROZEN'])))

def source_rows(headers,rows,kind,ids):
    for id in ids:
        ss=SOURCES[id]
        if kind==3:row3(record_type='SOURCE_REFERENCE',record_key=id,status=ss['relative_path'],NA_reason='sha256='+ss['sha256'],source_id=id)
        elif kind==4:row4(record_type='SOURCE_REFERENCE',entity=id,NA_or_boundary=ss['relative_path']+';sha256='+ss['sha256'],source_id=id)
        else:rows.append({z:('SOURCE_REFERENCE' if z=='record_type' else id if z=='gene' or z=='source_id' else ss['relative_path']+';sha256='+ss['sha256'] if z=='negative_or_missing' else None) for z in headers})
source_rows(H3,R3,3,['AXES','LR','ALL_LR','LR_COMPONENTS','FID_RAW','FID_FILTER','DONORS','MAIN3C','PROOF','TARGETS','FROZEN'])
source_rows(H4,R4,4,['ROI','DEPTH','RETENTION','SPATIAL_GENES','MARKERS','SPATIAL_PROGRAMS'])
source_rows(H5,R5,5,['CROSS','AXES','LR_COMPONENTS','SPATIAL_GENES','FROZEN'])
def primitive(v):
    if v is None:return None
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,np.floating):return float(v)
    if isinstance(v,np.bool_):return bool(v)
    return v
def sheet(name,headers,rows,numeric,source,description):
    rr=[[primitive(q.get(z)) for z in headers] for q in rows];data=pd.DataFrame(rr,columns=headers);write(data,name+'.tsv')
    return {'name':name,'headers':headers,'rows':rr,'types':['number' if z in numeric else 'text' for z in headers],'description':description,'source':SOURCES[source]['relative_path'],'source_sha256':SOURCES[source]['sha256']}
sheets=[sheet('S3_Official_SigXTalk',H3,R3,['PA','PRS','fidelity_raw','specificity_raw','fidelity_filtered','specificity_filtered','LR_probability','LR_P','value_or_n'],'AXES','New official run only.80representative paths, not all paths. MainFigure3LR18/path12/metric3/maxPRS matrix retained. Metrics raw/filter universes distinct; missingPRS not zero; gene role rows distinguish receptor and target. SOURCE_REFERENCE rows give portable paths/hashes. Counts and PA/PRS are not donor-level inference.'),sheet('S4_Spatial_Regions',H4,R4,['NP_mean','AF_mean','NP_detected_fraction','AF_detected_fraction','NP_n','AF_n','NP_raw_UMI','AF_raw_UMI','paired_contrast'],'SPATIAL_GENES','All cached structure pairs and their eligibility/positive states. Record type defines metric: ROI_QC means=medianUMI and fractions=retention; expression/program/marker means=log1pCP10k; n=detected bins for genes/markers and available bins for scores/QC. GEP5discordance retained. Noanimal/genotypeP. SOURCE_REFERENCE rows portable.'),sheet('S5_Cross_Evidence',H5,R5,['fixed_target_paths','mouse_NP_UMI','mouse_AF_UMI'],'CROSS','Seven discussed genes; separate roles, original indicators and negative/not-tested states, no composite score. Human detection and comparator denominators are read from the source table. Experiments remain source-unresolved/internal or NOT_ASSAYED. SOURCE_REFERENCE rows portable.')]
payload={'sheets':sheets,'meta':{'representative_path_count':80,'mandatory_path_keys_before_filling':mandatory_count,'selection_rule':'Union allmain3D12/main3E3; main3CwinningmaxPRSaxes; top3fixedPRSperdiscussedtarget; highestPAactiveNAperdiscussedtarget; top3ITGB1receptorpaths; explicitITGB1/HIF1A/PLOD2. Fillto80withfixedPRSdescending,tiesreceptor/SSC/target. Deduplicateexact3geneaxiskeys.','representative_paths_are_not_all_paths':True,'sources':SOURCES,'no_models_refit':True,'no_new_depth_draw':True,'old_V4_and_main_figures_modified':False,'new_S3_reuses_HE_raster_pixels':False}}
(O/'COMPACT_WORKBOOK_PAYLOAD.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
(O/'SOURCE_HASHES.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
(O/'FIGURE_EXPORT_CHECKS.json').write_text(json.dumps(EXPORTS,indent=2)+'\n')
(O/'CRITICAL_RESULTS.json').write_text(json.dumps({'tested':len(axes),'PAactive':len(active),'positivePRS':len(pos),'missingPRS':int(active.PRS.isna().sum()),'fixedfilter':len(fixed),'paired_donors':don.donor_id.nunique(),'cells':int(don.Freq.sum()),'LR':len(lr),'representative_paths':len(rep),'mainFig3D_all12_keys_in_subset':set(mainpaths.apply(k,axis=1)).issubset(rep.index),'mainFig3E_all3_keys_in_subset':set(mainmetrics.apply(k,axis=1)).issubset(rep.index),'all8spatial_structures':len(dep),'eligible_structures':len(elig),'positive_structures':int((elig.depth100_NP_minus_AF>0).sum()),'same_set_before':before,'same_set_after':after,'attenuation_percent':atten,'no_HE_pixels':True,'sheets':[{x['name']:len(x['rows'])} for x in sheets]},indent=2)+'\n')
print(json.dumps({'sheets':[(x['name'],len(x['rows']),len(x['headers'])) for x in sheets],'paths':len(rep),'mandatory':mandatory_count,'figures':EXPORTS,'attenuation':atten},indent=2))
