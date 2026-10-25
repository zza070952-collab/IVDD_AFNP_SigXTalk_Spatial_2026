#!/usr/bin/env python3
"""Build manuscript communication/spatial panels from frozen aggregate outputs.

No model training, ROI fitting, additional depth draws, or inferential bin tests.
Use --review-root and --project-root for private, read-only cache trees.
All generated artifacts go to a separate --output-root. Private H&E pixels and
ROI coordinates are loaded only at execution and are never bundled here.
"""
import argparse, json, hashlib, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from PIL import Image

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

AP=argparse.ArgumentParser(description=__doc__)
AP.add_argument('--review-root', type=Path, default=env_path('IVDD_V3_ROOT'))
AP.add_argument('--project-root', type=Path, default=env_path('IVDD_V4_ROOT'))
AP.add_argument('--roi-json', type=Path, default=env_path('IVDD_ROI_JSON'))
AP.add_argument('--output-root', type=Path, default=env_path('IVDD_OUTPUT_ROOT'))
args=AP.parse_args()
if not args.review_root or not args.project_root or not args.output_root:
    AP.error('Set --review-root/IVDD_V3_ROOT, --project-root/IVDD_V4_ROOT and --output-root/IVDD_OUTPUT_ROOT')
V=args.review_root.expanduser().resolve(); P=args.project_root.expanduser().resolve()
args.roi_json=(args.roi_json or V/'07_public_repo/configs/frozen_spatial_rois.json').expanduser().resolve()
W=separate_output(args.output_root, V, P, args.roi_json)
T=W/'03_tables/communication_spatial'; R=W/'08_review/communication_spatial'
T.mkdir(parents=True,exist_ok=True);R.mkdir(parents=True,exist_ok=True)
S=P/'01_sources_readonly/spatial_reference'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
 'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
 'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white'})
NP='#268C92'; AF='#C97B3A'; CTR='#426D9A'; CKO='#A45A88'; GREY='#77818B'
capcol={'CTR':CTR,'cKO':CKO}; regcol={'NP_CORE':NP,'AF_CORE':AF}
panels=[]; manifest=[]; export_checks=[]

def read(path):
    d=pd.read_csv(path,sep='\t');manifest.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size,'role':'read_only_frozen_source'})
    return d
def save(d,name):
    path=T/name;d.to_csv(path,sep='\t',index=False,na_rep='NA');return path
def panel(figid,pid,name,data,transform,unit):
    path=save(data,figid+'_'+pid+'_'+name+'.tsv')
    panels.append({'figure':figid,'panel_id':pid,'source_data_file':str(path.relative_to(W)),'n_rows':len(data),'transform':transform,'unit':unit,'independent_n':f'{donors.donor_id.nunique()} pooled human donors' if figid in ['Figure3','S3'] else f'structures nested in {len(caps)} captures; animal IDs unavailable','statistical_test':'none; descriptive'})
    return path
def lab(ax,letter,title):
    ax.set_title(title,loc='left',pad=9,fontweight='normal')
    ax.text(-.13,1.06,letter,transform=ax.transAxes,fontsize=13,fontweight='bold',va='bottom')
def clean(ax):
    for sp in ax.spines.values():sp.set_visible(False)
    ax.tick_params(length=0)
def export(fig,name,folder):
    folder.mkdir(parents=True,exist_ok=True)
    fig.canvas.draw()
    ren=fig.canvas.get_renderer();outside=[]
    fw,fh=fig.canvas.get_width_height()
    nonrendered_ticks=set()
    for ax in fig.axes:
        for axis,limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            lo,hi=sorted(limits)
            for tick in list(axis.get_major_ticks())+list(axis.get_minor_ticks()):
                if not ax.axison or not lo-1e-10 <= tick.get_loc() <= hi+1e-10:
                    nonrendered_ticks.update([tick.label1,tick.label2])
    for tx in fig.findobj(matplotlib.text.Text):
        if not tx.get_visible() or not tx.get_text():continue
        if tx in nonrendered_ticks:continue
        if tx.axes is not None and not tx.axes.axison and tx in tx.axes.get_xticklabels()+tx.axes.get_yticklabels():continue
        box=tx.get_window_extent(ren)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fw+1 or box.y1 > fh+1:
            outside.append(tx.get_text())
    for ext in ['pdf','svg','png','tiff']:
        kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
        fig.savefig(folder/(name+'.'+ext),dpi=300,**kw)
    export_checks.append({'figure':name,'width_mm':round(fig.get_figwidth()*25.4,3),'height_mm':round(fig.get_figheight()*25.4,3),'dpi':300,'minimum_font_pt':min(x.get_fontsize() for x in fig.findobj(matplotlib.text.Text) if x.get_visible() and x.get_text()),'text_outside_canvas':outside})
    plt.close(fig)

st=V/'02_sigxtalk_official_new/formal_new/tables'
axes=read(st/'ALL_TESTED_NEW_AXES_WITH_PRS.tsv')
lr=read(st/'AF_TO_NP_CANDIDATE_LR.tsv')
alllr=read(st/'CELLCHAT_ALL_SOURCE_TARGET_LR.tsv')
components=read(st/'LR_COMPLEX_COMPONENTS.tsv')
targets=read(st/'TARGET_GENE_SELECTION_ALL.tsv')
donors=read(st/'DONOR_GROUP_INPUT_COUNTS.tsv')
frozen=read(V/'01_human_programs/FROZEN_K5_PROGRAMS.tsv')
testability=read(st/'CANDIDATE_TESTABILITY.tsv')
metrics=read(st/'NEW_REPRESENTATIVE_OFFICIAL_FID_SPE.tsv')
rawmetrics=read(V/'02_sigxtalk_official_new/formal_new/prs_outputs/FID_SPE_RAW_UNIVERSE_OFFICIAL_NEW.tsv')
filteredmetrics=read(V/'02_sigxtalk_official_new/formal_new/prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv')
candidate=read(V/'04_integrated/MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv')
active=axes[axes.PA_threshold_pass].copy();positive=axes[axes.positive_official_PRS_returned].copy();filtered=axes[axes.PRS_relative_max_filter_pass].copy()
assert not axes.empty and not axes.duplicated(['Receptor','TF','TG']).any()
assert set(filtered.index).issubset(positive.index)
assert set(positive.index).issubset(active.index)
assert donors.Freq.ge(0).all()
run_counts=pd.DataFrame([('tested_axes',len(axes)),('PA_gt_0.75',len(active)),('positive_PRS_returned',len(positive)),('active_PRS_NA',active.PRS.isna().sum()),('fixed_relative_filter',len(filtered)),('AF_NP_LR_records',len(lr)),('AF_NP_LR_algorithm_P_lt_0.05',(lr.pval<.05).sum())],columns=['state','n_records'])
save(run_counts,'S5_S6_COMMUNICATION_COUNTS.tsv')
recsummary=axes.groupby('Receptor').agg(tested_axes=('TG','size'),PA_active_axes=('PA_threshold_pass','sum'),positive_PRS_axes=('positive_official_PRS_returned','sum'),filtered_axes=('PRS_relative_max_filter_pass','sum')).reset_index()
recsummary=recsummary.merge(components.groupby('receptor_component').size().rename('LR_component_records'),left_on='Receptor',right_index=True,how='left');recsummary['LR_component_records']=recsummary.LR_component_records.fillna(0).astype(int)
save(recsummary,'S6_RECEPTOR_NODE_AXES_AND_LR_COMPONENTS.tsv')
for source,name in [(lr,'S5_ALL_AF_TO_NP_LR.tsv'),(alllr,'S5_ALL_SOURCE_TARGET_LR.tsv'),(targets,'S5_FULL_TARGET_SELECTION_UNIVERSE.tsv'),(components,'S5_LR_COMPLEX_COMPONENTS.tsv'),(donors,'S5_DONOR_GROUP_COUNTS.tsv'),(testability,'S6_CANDIDATE_TESTABILITY.tsv'),(metrics,'S6_REPRESENTATIVE_OFFICIAL_METRICS.tsv')]:save(source,name)

geps=[]
for gene in frozen.loc[frozen.program.eq('GEP2')].sort_values('rank').gene:
    a=axes[axes.TG.eq(gene)]
    tr=targets[targets.gene.eq(gene)]
    geps.append({'gene':gene,'selected_target':bool(tr.selected_target.any()),'tested_axes':len(a),'active_axes':int(a.PA_threshold_pass.sum()),'positive_PRS_axes':int(a.positive_official_PRS_returned.sum()),'active_PRS_NA_axes':int((a.PA_threshold_pass & a.PRS.isna()).sum()),'filtered_axes':int(a.PRS_relative_max_filter_pass.sum()),'status':'not_tested' if len(a)==0 else 'tested'})
gep=pd.DataFrame(geps);save(gep,'S6_ALL_25_GEP2_TESTABILITY_AND_PRS_STATES.tsv')

# Figure 3: measured inputs, transparent metric definitions and finite networks.
fig=plt.figure(figsize=(180/25.4,237/25.4))
a=fig.add_axes([.07,.738,.22,.205]);b=fig.add_axes([.56,.735,.33,.21])
c=fig.add_axes([.10,.405,.35,.215]);d=fig.add_axes([.59,.405,.35,.215])
e=fig.add_axes([.12,.095,.33,.20]);f=fig.add_axes([.64,.095,.30,.20])
edge=alllr[(alllr.source=='AF_stromal') & alllr.target.isin(['NP_GEP2_high','NP_GEP2_low'])].groupby('target').agg(n_LR=('prob','size'),n_algorithm_P005=('pval',lambda x:int((x<.05).sum())),sum_CellChat_probability=('prob','sum')).reset_index()
panel('Figure3','A','sender_receiver_edges',edge,'edge width = count of CellChat algorithm P<0.05 records; pooled input','pooled LR records')
pos={'AF_stromal':(.17,.5),'NP_GEP2_high':(.83,.77),'NP_GEP2_low':(.83,.23)}
for row in edge.itertuples():
    ar=FancyArrowPatch(pos['AF_stromal'],pos[row.target],arrowstyle='-|>',mutation_scale=14,lw=.6+row.n_algorithm_P005/20,color=GREY,shrinkA=20,shrinkB=20)
    a.add_patch(ar);a.text(.5,.68 if 'high' in row.target else .30,str(row.n_algorithm_P005),ha='center',fontsize=9)
for node,(xx,yy) in pos.items():
    a.scatter(xx,yy,s=650,c=AF if node=='AF_stromal' else NP,edgecolor='white',zorder=3)
    a.text(xx,yy+(.20 if node!='NP_GEP2_low' else -.20),{'AF_stromal':'AF','NP_GEP2_high':'NP GEP2-high','NP_GEP2_low':'NP background'}[node].replace('NP ','NP\n'),ha='center',va='center')
a.set(xlim=(0,1),ylim=(-.05,1.1));a.axis('off');lab(a,'A','Pooled LR context')
a.text(.5,-.08,'Width: LR records with P<0.05',ha='center',transform=a.transAxes)
sel=lr.sort_values(['prob','interaction_name'],ascending=[False,True]).head(18).copy()
sel['display_label']=sel.ligand+' → '+sel.receptor
panel('Figure3','B','top18_LR',sel,'top18 descending actual CellChat probability; bubble area = probability; color = untransformed algorithm P','LR record')
b.scatter(sel.prob,np.arange(len(sel)),s=sel.prob/sel.prob.max()*200,c=sel.pval,cmap='viridis_r',vmin=0,vmax=1,edgecolors='#444444',linewidths=.3)
b.set_yticks(np.arange(len(sel)));b.set_yticklabels(sel.display_label,fontsize=8);b.invert_yaxis();b.set_xlabel('CellChat probability');b.set_xlim(-.008,.215);lab(b,'B','Actual LR complexes')
b.grid(axis='x',lw=.35,color='#dddddd');cb=fig.colorbar(ScalarMappable(norm=colors.Normalize(0,1),cmap='viridis_r'),ax=b,fraction=.035,pad=.03);cb.set_label('Algorithm P',fontsize=8)
rec=['ITGB1','ITGA3','ITGA5','ITGAV','CD44','SDC4'];tg=['FN1','PLOD2','CCN2','PRSS23','COL1A1']
hm=[]
for rr in rec:
    for gg in tg:
        q=axes[(axes.Receptor==rr)&(axes.TG==gg)];pp=q[q.positive_official_PRS_returned]
        hm.append({'Receptor':rr,'Target':gg,'tested_axes':len(q),'positive_PRS_axes':len(pp),'active_PRS_NA_axes':int((q.PA_threshold_pass&q.PRS.isna()).sum()),'max_returned_positive_PRS':pp.PRS.max() if len(pp) else np.nan,'display_state':'not_tested' if not len(q) else 'tested_no_positive' if not len(pp) else 'positive_PRS'})
hm=pd.DataFrame(hm);panel('Figure3','C','receptor_target_max_PRS',hm,'maximum returned positive PRS over unique receptor/SSC/target keys; tested-no-positive shown0-colored square, untestedX','receptor × target; SSC axes pooled')
mat=hm.pivot(index='Receptor',columns='Target',values='max_returned_positive_PRS').reindex(index=rec,columns=tg)
cmap=plt.get_cmap('magma').copy();cmap.set_bad('#ebedef');im=c.imshow(mat,cmap=cmap,aspect='auto',vmin=0,vmax=hm.max_returned_positive_PRS.max())
for j,rr in enumerate(rec):
    for i,gg in enumerate(tg):
        z=hm[(hm.Receptor==rr)&(hm.Target==gg)].iloc[0]
        if z.display_state=='not_tested':c.text(i,j,'×',ha='center',va='center',color='#555555')
        elif z.display_state=='tested_no_positive':c.text(i,j,'0+',ha='center',va='center',color='#555555')
c.set_xticks(range(len(tg)));c.set_xticklabels(tg,rotation=45,ha='right');c.set_yticks(range(len(rec)));c.set_yticklabels(rec);lab(c,'C','Receiver target context')
cb=fig.colorbar(im,ax=c,fraction=.04,pad=.04);cb.set_label('Max positive PRS',fontsize=8)
# Twelve paths selected transparently by PRS; one per SSC/target, receptor keys retained.
paths=filtered.sort_values(['PRS','Receptor','TF','TG'],ascending=[False,True,True,True]).drop_duplicates(['TF','TG']).head(12).copy()
panel('Figure3','D','representative_paths',paths,'top12 fixed-filter PRS paths after one per SSC/target; linewidth proportional to PRS; no ligand layer','unique receptor–SSC–target key')
cols=[list(dict.fromkeys(paths.Receptor)),list(dict.fromkeys(paths.TF)),list(dict.fromkeys(paths.TG))]
ymaps=[dict(zip(names,np.linspace(.96,.04,len(names)))) for names in cols]
for row in paths.itertuples():
    ys=[ymaps[0][row.Receptor],ymaps[1][row.TF],ymaps[2][row.TG]]
    d.plot([0,1,2],ys,lw=.7+1.4*row.PRS/paths.PRS.max(),color='#639AA5',alpha=.5,zorder=1)
for xx,names in enumerate(cols):
    for name in names:
        yy=ymaps[xx][name];d.scatter(xx,yy,s=16,c=NP if xx==0 else '#657B87',zorder=3)
        d.text(xx+(-.06 if xx==0 else .06),yy,name,ha='right' if xx==0 else 'left',va='center',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','pad':.15})
d.set(xlim=(-.63,2.63),ylim=(-.03,1.05));d.axis('off');lab(d,'D','Representative predicted paths')
for xx,title in enumerate(['Receptor','SSC','Target']):d.text(xx,1.06,title,ha='center',fontsize=8)
metrics=metrics.copy();metrics['path_label']=metrics.Receptor+' / '+metrics.TF+' / '+metrics.TG
panel('Figure3','E','official_fidelity_specificity',metrics,'official Calculate_Fidelity_Matrix$all and Calculate_Specificity_Matrix$all; fixed-filter universe only','predicted path; no patient uncertainty')
y=np.arange(len(metrics));e.scatter(metrics.fidelity_target_all,y-.08,c=NP,s=35,label='Fidelity');e.scatter(metrics.specificity_receptor_all,y+.08,c=AF,s=35,label='Specificity')
e.set_yticks(y);e.set_yticklabels(metrics.path_label.str.replace(' / ','\n',regex=False),fontsize=8);e.invert_yaxis();e.set_xlim(0,.14);e.set_xlabel('Official normalized fraction');e.legend(loc='lower right',bbox_to_anchor=(1.07,-.39),ncol=2,handletextpad=.3,columnspacing=.7);lab(e,'E','Official metrics')
statecount=pd.DataFrame({'state':['HGNN tested','Not tested','Any positive PRS','Any active PRS NA','Fixed-filter path'],'n_genes':[int(gep.tested_axes.astype(bool).sum()),int((gep.tested_axes==0).sum()),int(gep.positive_PRS_axes.astype(bool).sum()),int(gep.active_PRS_NA_axes.astype(bool).sum()),int(gep.filtered_axes.astype(bool).sum())]})
panel('Figure3','F','GEP2_target_states',statecount,'nonexclusive counts among all25 fixed GEP2 members; individual-gene states inS3',f'gene; total {len(gep)} frozen genes')
f.barh(range(len(statecount)),statecount.n_genes,color=[CTR,'#D7D9DC',NP,'#B6B8BD',AF],height=.6)
f.set_yticks(range(len(statecount)));f.set_yticklabels(statecount.state);f.invert_yaxis();f.set_xlim(0,len(gep)+1);f.set_xticks(np.arange(0,len(gep)+1,5));f.set_xlabel('GEP2 members (nonexclusive states)')
for i,n in enumerate(statecount.n_genes):f.text(n+.4,i,str(n),va='center')
lab(f,'F',f'All {len(gep)} frozen members')
export(fig,'Figure3_Predicted_Communication',W/'04_figures/Figure3')

# Spatial original measurements and the already completed single fixed draw.
roi=json.loads(args.roi_json.read_text());caps=sorted(roi['captures'],key=lambda x:x['genotype']!='CTR')
sp=read(S/'SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz')
programs=read(S/'NP_PROGRAM_BY_REGION.tsv');allprograms=read(P/'02_analysis/spatial_allprograms/NP_PROGRAM_BY_REGION.tsv')
chk=programs.merge(allprograms,on=['capture_id','disc_id','region','program'],suffixes=('_old','_new'))
assert len(chk)==len(programs) and np.allclose(chk.mean_program_score_old,chk.mean_program_score_new,atol=1e-12,rtol=0)
depth=read(V/'03_spatial_gse287499/DEPTH_DISC_PAIR_COMPARISON.tsv')
depthvals=read(V/'03_spatial_gse287499/DEPTH_SENSITIVITY_VALUES.tsv')
retention=read(V/'03_spatial_gse287499/DISC_LEVEL_PROGRAM_AND_DEPTH.tsv')
genes=read(V/'03_spatial_gse287499/CANDIDATE_GENES_PER_DISC_REGION.tsv')
markers=read(S/'roi_validation/INDEPENDENT_MARKERS_BY_DISC_REGION.tsv')
marker_validation=read(S/'roi_validation/INDEPENDENT_MARKER_VALIDATION_BY_DISC.tsv')
coverage=read(S/'PROGRAM_GENE_SPATIAL_COVERAGE.tsv')
eligible=depth[depth.adequate_fixed_depth_pair].copy()
before=eligible.original_eligible_NP_minus_AF.mean();after=eligible.depth100_NP_minus_AF.mean();att=100*(1-after/before)
assert len(depth)==sum(len(c['discs']) for c in caps) and len(eligible)>0
npositive=int((eligible.depth100_NP_minus_AF>0).sum())
assert np.isfinite(att)
save(allprograms,'S8_ALL_FIVE_PROGRAMS_DISC_REGION.tsv')
for d0,nm in [(markers,'S7_ALL_INDEPENDENT_MARKERS.tsv'),(marker_validation,'S7_MARKER_COMPATIBILITY.tsv'),(retention,'S8_DEPTH_RETENTION_EACH_DISC_REGION.tsv'),(depth,'S8_DEPTH_DISC_PAIRS.tsv'),(depthvals,'S8_SINGLE_FIXED_DRAW_BIN_VALUES.tsv'),(genes,'S8_CANDIDATE_DISC_REGION_VALUES.tsv'),(coverage,'S7_ALL_PROGRAM_ORTHOLOG_COVERAGE.tsv')]:save(d0,nm)
images={}
image_rows=[]
for cap in caps:
    ip=S/'prepared48'/cap['capture_id']/'source_images'/cap['source_image_filename'];h=hashlib.sha256(ip.read_bytes()).hexdigest()
    assert h==cap['source_image_sha256']
    images[cap['capture_id']]=np.asarray(Image.open(ip))
    image_rows.append({'capture_id':cap['capture_id'],'source_GSM':cap['source_GSM'],'original_image_sha256':h,'width_px':cap['source_image_pixels'][0],'height_px':cap['source_image_pixels'][1],'microns_per_pixel':cap['microns_per_pixel'],'image_adjustment':'none','crop_main':'full original image','source':'GSE287499 supplementary original H&E; author ROI overlay','reuse':'Figure4A and SupplementaryS4A; same originals; no local contrast edits'})
save(pd.DataFrame(image_rows),'S7_IMAGE_INTEGRITY_PROVENANCE.tsv')
def outlines(ax,cap,which=None,numbers=True):
    for i,disc in enumerate(cap['discs']):
        if which and disc['disc_id']!=which:continue
        for reg in disc['regions']:ax.add_patch(Polygon(reg['vertices'],closed=True,fill=False,edgecolor=regcol[reg['region']],linewidth=.9))
        if numbers:
            vv=np.array(disc['regions'][0]['vertices']);ax.text(vv[:,0].mean(),vv[:,1].min()-22,disc['disc_id'].split('_')[-1],fontsize=8,ha='center',va='bottom',bbox={'facecolor':'white','edgecolor':'none','pad':.3})
def mapstyle(ax,cap,window=None,bar=500):
    w,h=cap['source_image_pixels'];x0,x1,y0,y1=window or (0,w,0,h)
    ax.set(xlim=(x0,x1),ylim=(y1,y0));ax.set_aspect('equal');ax.axis('off')
    start=x0+.07*(x1-x0);yy=y1-.07*(y1-y0);length=bar/cap['microns_per_pixel']
    ax.plot([start,start+length],[yy,yy],color='#222222',lw=1.5)
    ax.text(start+length/2,yy-.02*(y1-y0),str(bar)+' µm',ha='center',va='bottom',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','pad':.3})
fig=plt.figure(figsize=(180/25.4,237/25.4))
maps=[fig.add_axes([xx,.705,.205,.24]) for xx in [.045,.285,.525,.765]]
vmax=sp.GEP2_mean_log1p_CP10k.max()
for i,cap in enumerate(caps):
    aa=maps[i];aa.imshow(images[cap['capture_id']],origin='upper');outlines(aa,cap);mapstyle(aa,cap);aa.set_title(cap['genotype'],fontsize=9,pad=5)
    bb=maps[i+2];ss=sp[sp.capture_id.eq(cap['capture_id'])]
    bb.scatter(ss.x_fullres,ss.y_fullres,c=ss.GEP2_mean_log1p_CP10k,cmap='viridis',vmin=0,vmax=vmax,s=.5,marker='s',linewidths=0,rasterized=True);outlines(bb,cap);mapstyle(bb,cap);bb.set_title(cap['genotype'],fontsize=9,pad=5)
fig.text(.035,.97,'A',fontweight='bold',fontsize=13);fig.text(.10,.97,'Frozen H&E anatomy',fontsize=9)
fig.text(.515,.97,'B',fontweight='bold',fontsize=13);fig.text(.58,.97,'Original GEP2 score',fontsize=9)
fig.legend(handles=[Line2D([0],[0],color=NP,lw=1.5,label='NP core'),Line2D([0],[0],color=AF,lw=1.5,label='AF core')],loc='upper left',bbox_to_anchor=(.07,.697),ncol=2,handlelength=1.2,columnspacing=1.2)
cax=fig.add_axes([.58,.68,.32,.012]);cb=fig.colorbar(ScalarMappable(norm=colors.Normalize(0,vmax),cmap='viridis'),cax=cax,orientation='horizontal');cb.set_label('Mean log(1 + CP10k)',labelpad=2)
panel('Figure4','A','H_and_E_provenance',pd.DataFrame(image_rows),'full original images; frozen polygons; calibrated500µm scale bars',f'{len(caps)} captures; {len(depth)} structures; frozen core polygons')
vertices=read(V/'05_figures/source_data/Figure4ab_frozen_HE_vertices.tsv');panel('Figure4','A2','ROI_vertices',vertices,'unchanged frozen full-resolution coordinates','polygon vertex')
panel('Figure4','B','original_spatial_coordinates',sp[['spot_id','capture_id','genotype','disc_id','anatomical_region','x_fullres','y_fullres','n_counts','n_genes','GEP2_mean_log1p_CP10k']],f'all {len(sp)} actual QC-bin centers; common score color limits; no smoothing','48µm bin nested in capture')
c=fig.add_axes([.13,.425,.31,.17]);d=fig.add_axes([.61,.425,.32,.17])
mg=markers.groupby(['genotype','histology_core_region','marker'],sort=False).agg(mean_expression=('mean_log1p_per_10000','mean'),mean_detected_fraction=('detected_bin_fraction','mean'),n_disc_structures=('disc_id','nunique')).reset_index()
panel('Figure4','C','independent_markers',mg,'equal mean over4structures/capture; bubble area = detectionfraction; color=meanlog1pCP10k','marker × capture × core')
ms=['Krt19','T','Krt8','Acan','Chad','Clec3a'];cols=[('CTR','NP_CORE'),('CTR','AF_CORE'),('cKO','NP_CORE'),('cKO','AF_CORE')]
for j,g in enumerate(ms):
    for i,(gt,reg) in enumerate(cols):
        q=mg[(mg.genotype==gt)&(mg.histology_core_region==reg)&(mg.marker==g)].iloc[0]
        c.scatter(i,j,s=10+85*q.mean_detected_fraction,c=[q.mean_expression],vmin=0,vmax=mg.mean_expression.max(),cmap='magma',edgecolors='white',linewidths=.3)
c.set(xticks=range(4),yticks=range(6),ylim=(5.6,-.6));c.set_xticklabels(['CTR\nNP','CTR\nAF','cKO\nNP','cKO\nAF']);c.set_yticklabels(ms,fontstyle='italic');clean(c);lab(c,'C','Independent structural markers')
cb=fig.colorbar(ScalarMappable(norm=colors.Normalize(0,mg.mean_expression.max()),cmap='magma'),ax=c,fraction=.04,pad=.04);cb.set_label('Mean log(1 + CP10k)')
panel('Figure4','D','eight_disc_pairs',depth,'originalallcoremeans; each lineone structure; no P value',f'{len(depth)} structures nested in {len(caps)} captures')
for i,row in enumerate(depth.itertuples()):
    xx=np.array([0,1])+(i%4-1.5)*.055+(0 if row.genotype=='CTR' else .12)
    d.plot(xx,[row.original_all_AF_CORE,row.original_all_NP_CORE],color=capcol[row.genotype],marker='o' if row.genotype=='CTR' else '^',ms=4,lw=.8,alpha=.8)
d.set(xticks=[0,1],xticklabels=['AF core','NP core'],ylabel='GEP2 mean log(1 + CP10k)',xlim=(-.25,1.35));lab(d,'D',f'All {len(depth)} structures')
d.legend(handles=[Line2D([0],[0],color=CTR,marker='o',label=f'CTR ({int((depth.genotype=="CTR").sum())})'),Line2D([0],[0],color=CKO,marker='^',label=f'cKO ({int((depth.genotype=="cKO").sum())})')],loc='upper left',bbox_to_anchor=(-.05,-.20),ncol=2,handlelength=1.0,columnspacing=.6,handletextpad=.3)
e=fig.add_axes([.12,.105,.38,.215]);f=fig.add_axes([.68,.105,.21,.215])
plotrows=[]
for i,row in enumerate(eligible.itertuples()):
    off=(i-3)*.018;e.plot([0+off,1+off],[row.original_eligible_NP_minus_AF,row.depth100_NP_minus_AF],c=capcol[row.genotype],lw=.85,marker='o' if row.genotype=='CTR' else '^',ms=4,alpha=.8)
    plotrows.extend([{'capture_id':row.capture_id,'disc_id':row.disc_id,'stage':'original_same_eligible_bins','contrast':row.original_eligible_NP_minus_AF},{'capture_id':row.capture_id,'disc_id':row.disc_id,'stage':'single100UMIdraw','contrast':row.depth100_NP_minus_AF}])
e.plot([0,1],[before,after],'D-',c='#20252A',lw=2,ms=5,zorder=5);e.axhline(0,c='#888888',lw=.6,ls='--');e.set(xticks=[0,1],xticklabels=['Original eligible\nbins (≥100 UMI)','Single fixed draw\n100 UMI/bin'],ylabel='GEP2 NP − AF contrast',xlim=(-.17,1.17),ylim=(-.065,.39));lab(e,'E',f'Depth sensitivity: same {len(eligible)} structures')
e.text(.5,.93,f'{att:.1f}% attenuation',transform=e.transAxes,ha='center',fontweight='bold',fontsize=10)
e.text(.5,-.31,f'Mean: {before:.6f} → {after:.6f}; {npositive}/{len(eligible)} remain NP > AF',transform=e.transAxes,ha='center',fontsize=8)
dp=pd.DataFrame(plotrows);dp['descriptive_group_mean']=dp.stage.map({'original_same_eligible_bins':before,'single100UMIdraw':after});dp['attenuation_percent_same_eligible_set']=att
panel('Figure4','E','same_eligible_depth_contrasts',dp,'exact same eligible bins; one cached fixed-depth draw; black diamonds mean structure contrasts',f'{len(eligible)} eligible structures, {npositive} positive; no animal inference')
gn=['FN1','ITGB1','PLOD2','CCN2','PRSS23'];gg=genes[genes.human_gene.isin(gn)].groupby(['genotype','region','human_gene'],sort=False).agg(mean_expression=('mean_log1p_CP10k','mean'),mean_detection=('fraction_detected_bins','mean'),total_raw_UMI=('total_raw_UMI','sum'),n_disc_structures=('disc_id','nunique')).reset_index()
panel('Figure4','F','candidate_orthologs',gg,'disc-equal meanexpressionanddetection; originalunadjustedcounts; exactmousemappinginS7','gene×capture×core;4structures each')
for j,g in enumerate(gn):
    for i,(gt,reg) in enumerate(cols):
        q=gg[(gg.genotype==gt)&(gg.region==reg)&(gg.human_gene==g)].iloc[0]
        if q.mean_detection==0:f.scatter(i,j,s=15,facecolors='white',edgecolors='#888888',lw=.6)
        else:f.scatter(i,j,s=12+110*q.mean_detection,c=[q.mean_expression],vmin=0,vmax=gg.mean_expression.max(),cmap='magma',edgecolors='white',lw=.3)
f.set(xticks=range(4),yticks=range(5),ylim=(4.7,-.6));f.set_xticklabels(['CTR\nNP','CTR\nAF','cKO\nNP','cKO\nAF']);f.set_yticklabels(gn,fontstyle='italic');clean(f);lab(f,'F','Candidate orthologs')
cb=fig.colorbar(ScalarMappable(norm=colors.Normalize(0,gg.mean_expression.max()),cmap='magma'),ax=f,fraction=.055,pad=.05);cb.set_label('Mean log(1 + CP10k)')
f.text(.5,-.31,f'Hollow:0 detected; PRSS23:{int(genes.loc[genes.human_gene.eq("PRSS23"),"total_raw_UMI"].sum())} UMI total',ha='center',transform=f.transAxes,fontsize=8)
export(fig,'Figure4_Anatomy_and_Depth',W/'04_figures/Figure4')

# S3: complete LR matrix, 25-member status, full normalization universes.
fig=plt.figure(figsize=(180/25.4,215/25.4))
a=fig.add_axes([.24,.20,.71,.73]);recs=sorted(lr.receptor.unique());ligs=sorted(lr.ligand.unique())
for q in lr.itertuples():a.scatter(recs.index(q.receptor),ligs.index(q.ligand),s=8+100*q.prob/lr.prob.max(),c=[q.pval],cmap='viridis_r',vmin=0,vmax=1,edgecolors='#666666',lw=.2)
a.set(xticks=range(len(recs)),yticks=range(len(ligs)),xlim=(-.6,len(recs)-.4),ylim=(len(ligs)-.4,-.6));a.set_xticklabels(recs,rotation=65,ha='right');a.set_yticklabels(ligs,fontstyle='italic');clean(a);lab(a,'A',f'All {len(lr)} AF→NP GEP2-high LR records')
panel('S3','A','all_LR_matrix',lr,'allLRrecords; areaCellChatprobability; coloralgorithmP; no truncation','LRrecord')
fig.text(.24,.025,'Bubble area:CellChat probability; color:algorithm P (0 yellow → 1 purple).',fontsize=8)
export(fig,f'S3_Complete_{len(lr)}_LR_Complexes',W/'06_supplementary/S3')
fig=plt.figure(figsize=(180/25.4,190/25.4))
b=fig.add_axes([.17,.20,.37,.68]);cm=gep.set_index('gene')[['tested_axes','active_axes','positive_PRS_axes','active_PRS_NA_axes','filtered_axes']]
im=b.imshow(np.log10(1+cm),aspect='auto',cmap='Blues',vmin=0,vmax=np.log10(1+cm.to_numpy().max()))
for j,g in enumerate(gep.gene):
    if not gep.iloc[j].tested_axes:b.text(2,j,'NOT TESTED',ha='center',va='center',fontsize=8,color='#555555')
b.set(xticks=range(5),yticks=range(len(gep)));b.set_xticklabels(['Tested','PA>.75','PRS+','PRS NA*','Filtered'],rotation=45,ha='right');b.set_yticklabels(gep.gene,fontsize=8);lab(b,'B','All GEP2 member states')
gax=fig.add_axes([.17,.055,.37,.013]);gcb=fig.colorbar(im,cax=gax,orientation='horizontal');gcb.set_label('log10(1 + number of axes)',labelpad=1)
panel('S3','B','all25_target_status',gep,'heatmaplog10(1+axiscount); zerosareactualcounts; targetnot-tested explicit; NA*=activePRSNA',f'{len(gep)} frozen GEP2 genes')
c=fig.add_axes([.74,.58,.22,.31]);c.barh(range(len(run_counts)),run_counts.n_records,color=[GREY,CTR,NP,'#B7B8BD',AF,GREY,GREY]);c.set_yticks(range(len(run_counts)));c.set_yticklabels(['Tested','PA>.75','PRS+','Active PRS NA','Fixed filter','AF→NP LR','LR P<.05'],fontsize=8);c.invert_yaxis();c.set_xscale('log');c.set_xlim(50,400000);c.set_xticks([100,10000]);c.set_xticklabels(['100','10,000']);c.set_xlabel('Records (log scale)');lab(c,'C','Different denominators')
panel('S3','C','run_counts',run_counts,'actualdistinctdenominators,logscaleaxis','record')
d=fig.add_axes([.74,.20,.22,.22]);m0=rawmetrics[['Receptor','SSC','Target','fidelity_target_all','specificity_receptor_all']].merge(filteredmetrics[['Receptor','SSC','Target','fidelity_target_all','specificity_receptor_all']],on=['Receptor','SSC','Target'],suffixes=('_raw','_filtered'),validate='one_to_one')
for i,(title,col) in enumerate([('Fidelity','fidelity_target_all'),('Specificity','specificity_receptor_all')]):
    for j,(lab0,cc) in enumerate([('Raw-positive','_raw'),('Fixed-filter','_filtered')]):
        yy=m0[col+cc].median();d.scatter(i+(j-.5)*.15,yy,c=NP if j==0 else AF,s=35,label=lab0 if i==0 else None)
d.set_xticks([0,1]);d.set_xticklabels(['Fidelity','Specificity'],rotation=25,ha='right');d.set_ylabel('Median fraction');d.legend(loc='lower center',bbox_to_anchor=(.4,-.75),ncol=1);lab(d,'D','Separate universes')
panel('S3','D','official_universe_comparison',m0,f'{len(m0)} same keys; official normalization separately raw-positive and filtered; dots medians','path; no biological error bar')
export(fig,'S3_Communication_Completeness',W/'06_supplementary/S3')

# Separate S3 continued: finite paths, source-receptor mapping and evidence values.
fig=plt.figure(figsize=(180/25.4,175/25.4))
a=fig.add_axes([.12,.60,.80,.32]);cc=components.groupby(['ligand','receptor_complex']).agg(components=('receptor_component',lambda x:';'.join(sorted(set(x))))).reset_index();fn=cc[cc.ligand=='FN1']
for i,q in enumerate(fn.itertuples()):a.text(.01,1-(i+.5)/len(fn),q.ligand+' → '+q.receptor_complex+' : '+q.components,fontsize=9,va='center',transform=a.transAxes)
a.axis('off');lab(a,'E','All FN1 receptor-complex keys and component membership')
panel('S3','E','FN1_complex_keys',fn,'originalcomplexkeys retained; componentmembershipseparate','LRcomplex')
b=fig.add_axes([.20,.16,.71,.30]);gs=['FN1','ITGB1','PLOD2','CCN2','PRSS23'];es=candidate.set_index('gene').loc[gs,['GEP2_frozen_member','n_tested_axes','n_filtered_PRS_axes','Tu_NP_donors_detected','mouse_NP_total_UMI','mouse_AF_total_UMI']].reset_index()
tab=[]
for q in es.itertuples():tab.append([q.gene,'yes' if q.GEP2_frozen_member else 'no',str(q.n_tested_axes),str(q.n_filtered_PRS_axes),str(q.Tu_NP_donors_detected)+'/'+str(int(candidate.set_index('gene').loc[q.gene,'Tu_NP_donors'])),str(int(q.mouse_NP_total_UMI)),str(int(q.mouse_AF_total_UMI))])
b.axis('off');table=b.table(cellText=tab,colLabels=['Gene','GEP2','Target\naxes','Filtered\ntarget axes','Tu NP\ndetected','Mouse NP\nUMI','Mouse AF\nUMI'],loc='center',cellLoc='center',colWidths=[.12,.10,.17,.17,.14,.15,.15]);table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,2)
for (i,j),cell in table.get_celld().items():cell.set_edgecolor('white');cell.set_facecolor('#E9EFF2' if i==0 else ('#f6f7f8' if i%2 else 'white'))
lab(b,'F','Role-specific evidence; no combined score')
panel('S3','F','candidate_evidence_values',es,'target-role counts shown; receptor-node axes separate; original values, not ranked or summed','targetrole distinctfromreceptorrole')
fig.text(.12,.07,f'ITGB1: {int(axes.Receptor.eq("ITGB1").sum()):,} receptor-node axes; receptor and target roles separate.\nCandidate coverage and counterevidence are read from the private input tables.',fontsize=9,linespacing=1.7)
export(fig,'S3_Continued_Complexes_and_Counterevidence',W/'06_supplementary/S3')

# S4: every structure at a fixed anatomy-derived crop, marker and QC completeness.
fig=plt.figure(figsize=(180/25.4,237/25.4));discs=[];crop_rows=[]
for cap in caps:
    for ds in cap['discs']:discs.append((cap,ds))
for i,(cap,ds) in enumerate(discs):
    ax=fig.add_axes([.04+(i%4)*.24,.81-(i//4)*.17,.22,.13]);ax.imshow(images[cap['capture_id']],origin='upper');vv=np.vstack([r['vertices'] for r in ds['regions']]);margin=45;window=[vv[:,0].min()-margin,vv[:,0].max()+margin,vv[:,1].min()-margin,vv[:,1].max()+margin]
    outlines(ax,cap,ds['disc_id'],False);mapstyle(ax,cap,window,bar=200);ax.set_title(ds['disc_id'],fontsize=9,pad=3);crop_rows.append({'capture_id':cap['capture_id'],'disc_id':ds['disc_id'],'crop_x0':window[0],'crop_x1':window[1],'crop_y0':window[2],'crop_y1':window[3],'crop_rule':'all3corepolygonboundingbox+45px; noexpressionselection','brightness_contrast_gamma':'original unchanged','scale_bar_um':200})
fig.text(.035,.965,'A',fontsize=13,fontweight='bold');fig.text(.09,.965,f'All {len(discs)} frozen structures and bilateral AF / central NP cores',fontsize=9)
panel('S4','A','all8_anatomical_crops',pd.DataFrame(crop_rows),'allstructures fixedanatomyboundingbox+45px; sameoriginalimages as4A;200µmcalibratedbar','structure nestedincapture')
b=fig.add_axes([.12,.425,.34,.15]);m=markers.copy();m['disc_region']=m.disc_id+' '+m.histology_core_region.str.replace('_CORE','',regex=False)
dsorder=[(ds['disc_id'],reg) for cp,ds in discs for reg in ['NP_CORE','AF_CORE']]
mm=m.pivot(index='marker',columns='disc_region',values='mean_log1p_per_10000').reindex(index=ms,columns=[di+' '+reg.replace('_CORE','') for di,reg in dsorder]);im=b.imshow(mm,aspect='auto',cmap='magma',vmin=0,vmax=mm.max().max());b.set(xticks=range(16),yticks=range(6));b.set_xticklabels([s.replace('_',' ') for s in mm.columns],rotation=90,fontsize=8);b.set_yticklabels(ms,fontstyle='italic');lab(b,'B','Complete independent markers')
panel('S4','B','all_marker_regions',m,'rawmeanlog1pCP10k beforedisplay; no rowstandardization','6markers×16disc-regions')
mcax=fig.add_axes([.49,.425,.012,.15]);mcb=fig.colorbar(im,cax=mcax);mcb.set_label('Mean log(1 + CP10k)')
c=fig.add_axes([.70,.425,.25,.15]);core=sp[sp.anatomical_region.isin(['NP_CORE','AF_CORE'])].copy();qsum=core.groupby(['genotype','anatomical_region']).agg(n_bins=('spot_id','size'),median_UMI=('n_counts','median'),median_detected_genes=('n_genes','median')).reset_index()
for i,(gt,reg) in enumerate(cols):
    zz=core[(core.genotype==gt)&(core.anatomical_region==reg)]
    c.scatter(zz.n_counts,zz.n_genes,s=9,c=regcol[reg],marker='o' if gt=='CTR' else '^',alpha=.5,rasterized=True,linewidths=0)
c.set_xscale('log');c.set_yscale('log');c.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:format(x,',.0f')));c.yaxis.set_major_formatter(FuncFormatter(lambda x,pos:format(x,',.0f')));c.set(xlabel='Original UMI/bin (log)',ylabel='Detected genes/bin (log)');lab(c,'C','Depth and detected genes')
panel('S4','C','bin_UMI_detected_genes',core[['spot_id','capture_id','genotype','disc_id','anatomical_region','n_counts','n_genes']],'allcorebins; actualtotalUMIversusdetectedgenes; bothaxeslogscale; tealNP orangeAF; circleCTR trianglecKO','bin nestedincapture; noPvalue')
d=fig.add_axes([.13,.12,.34,.16]);ret=retention.copy();ret['label']=ret.disc_id+' '+ret.region.str.replace('_CORE','',regex=False);ret=ret.set_index('label').loc[[di+' '+reg.replace('_CORE','') for di,reg in dsorder]].reset_index()
d.bar(range(16),ret.retention_fraction,color=[regcol[x] for x in ret.region]);d.axhline(1,lw=.5,c='#888888');d.set(xticks=range(16),ylabel='Fraction retained ≥100 UMI',ylim=(0,1.08));d.set_xticklabels(ret.label,rotation=90,fontsize=8);lab(d,'D','Retention in every structure/core')
panel('S4','D','all_disc_region_retention',ret,'eligible>=100originalUMI; noadditionaldraws; adequate≥5bins/core','16disc-regions')
e=fig.add_axes([.70,.12,.25,.16]);pp=allprograms.pivot(index=['capture_id','genotype','disc_id','program'],columns='region',values='mean_program_score').reset_index();pp['NP_minus_AF']=pp.NP_CORE-pp.AF_CORE
pa=pp.pivot(index='disc_id',columns='program',values='NP_minus_AF').reindex(index=[ds['disc_id'] for cap,ds in discs],columns=['GEP1','GEP2','GEP3','GEP4','GEP5']);bound=np.abs(pa.to_numpy()).max();im=e.imshow(pa,aspect='auto',cmap='RdBu_r',vmin=-bound,vmax=bound);e.set(xticks=range(5),yticks=range(8));e.set_xticklabels(pa.columns,rotation=45,ha='right');e.set_yticklabels(pa.index);lab(e,'E','All five frozen programs')
pcax=fig.add_axes([.70,.048,.25,.010]);cb=fig.colorbar(im,cax=pcax,orientation='horizontal');cb.set_label('NP − AF mean score',labelpad=1)
panel('S4','E','all5_program_structure_contrasts',pp,'cachedsparseprojectionofall5fixedprograms; sameorthologmapping; GEP2/GEP5 reproduce cached means at1e-12; nofit',f'{len(discs)} structures × {allprograms.program.nunique()} programs')
export(fig,'S4_All_Structures_Markers_Depth_Programs',W/'06_supplementary/S4')

pd.DataFrame(panels).to_csv(R/'PANEL_SOURCE_DATA_MAP.tsv',sep='\t',index=False)
pd.DataFrame(manifest).drop_duplicates('path').to_csv(R/'INPUT_MANIFEST.tsv',sep='\t',index=False)
(R/'FIGURE_EXPORT_CHECKS.json').write_text(json.dumps(export_checks,indent=2)+'\n')
facts={'patients':int(donors.donor_id.nunique()),'cells':int(donors.Freq.sum()),'AF_NP_LR':len(lr),'AF_NP_LR_p005':int((lr.pval<.05).sum()),'HGNN_tested':len(axes),'PA_active':len(active),'positive_PRS':len(positive),'filtered_PRS':len(filtered),'active_PRS_NA':int(active.PRS.isna().sum()),'spatial_structures':len(depth),'eligible_structures':len(eligible),'depth_positive_structures':int((eligible.depth100_NP_minus_AF>0).sum()),'original_same_eligible_mean_NP_minus_AF':before,'fixed100_mean_NP_minus_AF':after,'attenuation_percent':att,'GEP2_GEP5_previous_summary_max_abs_difference':float(np.max(np.abs(chk.mean_program_score_old-chk.mean_program_score_new))),'no_model_training':True,'no_new_depth_draw':True,'no_ROI_changes':True,'no_genotype_or_animal_P_values':True,'no_GitHub_writes':True}
(R/'RECOMPUTED_FACTS.json').write_text(json.dumps(facts,indent=2)+'\n')
print(json.dumps(facts,indent=2))
