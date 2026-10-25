#!/usr/bin/env python3
"""Rebuild the compact human figure and two static table payloads from external caches.

Set IVDD_REVIEW_ROOT and IVDD_V4_ROOT, or pass --source-root and --v4-root.
No manuscript paragraphs, legends, workbooks or data payloads are bundled.
"""
import argparse,hashlib,json,os,textwrap
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib import colors
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import fitz

NP='NP_matrix_marker_supported';IM='immune_marker_supported'
COMPONENTS={'COL2A1':'Matrix','PRELP':'Matrix','COMP':'Matrix','FMOD':'Matrix',
            'LCN2':'Stress-associated','SELE':'Stress-associated','MMP3':'Stress-associated','NOS2':'Stress-associated'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,
    'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'legend.frameon':False,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
    'xtick.major.size':3,'ytick.major.size':2,'pdf.fonttype':42,'svg.fonttype':'none'})

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_json(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source-root',type=Path,default=os.environ.get('IVDD_REVIEW_ROOT'))
    ap.add_argument('--v4-root',type=Path,default=os.environ.get('IVDD_V4_ROOT'))
    ap.add_argument('--eligibility-table',type=Path,help='Override V4_CACHE/00_inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv')
    ap.add_argument('--output-root',type=Path,required=True)
    a=ap.parse_args()
    if a.source_root is None:ap.error('Provide --source-root or IVDD_REVIEW_ROOT; no frozen result tables are bundled.')
    if a.v4_root is None and a.eligibility_table is None:ap.error('Provide --v4-root/IVDD_V4_ROOT or --eligibility-table.')
    src=a.source_root.resolve();out=a.output_root.resolve()
    for root in [src]+([a.v4_root.resolve()] if a.v4_root else []):
        if out==root or root in out.parents:ap.error('--output-root must be outside the read-only input caches.')
    eligp=a.eligibility_table or a.v4_root/'00_inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv'
    for path in [src,eligp]:
        if not path.exists():ap.error('Required input is absent: '+path.name+' (see README cache contract).')
    out.mkdir(parents=True,exist_ok=True)
    figs=out/'figures';tabs=out/'tables';figs.mkdir(exist_ok=True);tabs.mkdir(exist_ok=True)
    manifest=[]
    def load(rel):
        p=src/rel;manifest.append({'source_project':'frozen V3 review','source_relative_path':rel,'bytes':p.stat().st_size,'sha256':digest(p)})
        return pd.read_csv(p,sep='\t')
    H='01_human_programs/'
    donors=load(H+'COHORT_AND_DONOR_MANIFEST.tsv');counts=load(H+'FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv')
    scores=load(H+'FIGURE2_DONOR_PROGRAM_SOURCE.tsv');members=load(H+'FROZEN_K5_PROGRAMS.tsv')
    coverage=load(H+'PROGRAM_GENE_COVERAGE.tsv');genes=load(H+'TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv')
    brief=load(H+'ALL_PROGRAMS_BRIEF.tsv');support=load(H+'GEP2_DONOR_SUPPORT.tsv')
    tissue=load('05_figures/Figure1_Figure2/Figure1a_donor_tissue_counts.tsv')
    inputs=load('02_sigxtalk_official_new/formal_new/tables/DONOR_GROUP_INPUT_COUNTS.tsv')
    pairs=load('07_public_repo/resources/metadata/SWAHN_VERIFIED_DONOR_PAIRING.tsv')
    elig=pd.read_csv(eligp,sep='\t')
    manifest.append({'source_project':'frozen V4 inventory','source_relative_path':'00_inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv','bytes':eligp.stat().st_size,'sha256':digest(eligp)})
    assert members.groupby('program').size().eq(25).all() and members.program.nunique()==5
    PROGS=sorted(members.program.unique())
    TU=donors.loc[donors.cohort.eq('Tu_GSE165722'),'donor_id'].drop_duplicates().sort_values().tolist()
    DISC=donors.loc[donors.cohort.eq('Swahn_GSE230808'),'donor_id'].drop_duplicates().sort_values().tolist()
    nt,nd=len(TU),len(DISC);n_members=int(members.groupby('program').size().iloc[0]);n_pairs=int(inputs.donor_id.nunique())
    t=donors[donors.cohort.eq('Tu_GSE165722')].set_index('donor_id').reindex(TU).reset_index()
    ct=counts[counts.cohort.eq('Tu_GSE165722')].pivot(index='donor_id',columns='identity',values='n_cells').reindex(TU).fillna(0)
    tr=support[support.cohort.eq('Tu_GSE165722')].set_index('donor_id').reindex(TU)
    assert t.QC_cells.sum()==ct.sum(axis=1).sum() and ct[NP].sum()<=t.QC_cells.sum()
    comparable=tr.context_donor_median.notna()
    positive=comparable & tr.NP_donor_median.gt(tr.context_donor_median)
    n_comparable=int(comparable.sum());n_positive=int(positive.sum())
    absent_comparators=tr.index[~comparable].tolist()
    grade_status=' / '.join(brief.grade_interpretation.dropna().drop_duplicates())
    g5_status=brief.set_index('program').loc['GEP5','expression_interpretation']
    g5_scope=brief.set_index('program').loc['GEP5','interpretation_scope']
    s1a=t[['donor_id','gsm','submitted_cells','QC_cells']].copy();s1a['QC_retention_fraction']=s1a.QC_cells/s1a.submitted_cells
    s1b=pd.DataFrame({'donor_id':TU,'NP_supported_cells':ct[NP].astype(int),'immune_supported_cells':ct[IM].astype(int),
        'comparison_assessable':comparable,'NP_gt_immune':np.where(comparable,np.where(positive,'yes','no'),'NA'),
        'NA_reason':np.where(comparable,'none','no immune-supported cells; no comparator score')}).reset_index(drop=True)
    sel=list(COMPONENTS)
    s1c=genes[genes.program.eq('GEP5')&genes.identity.eq(NP)&genes.original_gene.isin(sel)].copy()
    s1c['component']=s1c.original_gene.map(COMPONENTS)
    s1c=s1c.merge(members[members.program.eq('GEP5')][['gene','rank']].rename(columns={'gene':'original_gene','rank':'frozen_rank'}),on='original_gene',validate='many_to_one')
    assert set(sel).issubset(set(members.loc[members.program.eq('GEP5'),'gene']))
    assert len(s1c)==len(sel)*nt
    for label,df in [('A',s1a),('B',s1b),('C',s1c)]:df.to_csv(figs/f'Supplementary_Figure_S1{label}_source.tsv',sep='\t',index=False,na_rep='NA')
    contract={'figure':'Supplementary Figure S1','core_conclusion':f'Actual donor coverage provides {n_comparable} assessable GEP2 comparisons ({n_positive} NP>immune); fixed GEP5 members provide component context.',
        'archetype':'quantitative grid','width_mm':180,'height_mm':210,'min_font_pt':8,
        'panels':{'A':'Tu QC retention by donor, not another main-figure composition plot','B':f'all{nt}donors with actual NP/immune cell coverage; {n_comparable} assessable; absent comparators remain NA',
                  'C':'exact8previously discussed matrix/stress members of frozenGEP5; not another five-program heatmap'},
        'selection_rule':'Fixed illustrative GEP5 subset before plotting: COL2A1,PRELP,COMP,FMOD(matrix);LCN2,SELE,MMP3,NOS2(stress-associated). All25GEP5members remain inTableS2 overview.',
        'statistical_unit':'donor; cells denote coverage; no inferential tests','bounds':['No grade inference','No new fit or rescoring','GEP5 frozen review: '+g5_status,'Identity support is operational and not histological certification','AUCell/UCell scale magnitudes are not compared'],
        'exports':['PDF','SVG','PNG300dpi'],'original_main_figures_modified':False}
    write_json(out/'FIGURE_S1_CONTRACT.json',contract)
    fig=plt.figure(figsize=(180/25.4,210/25.4));titles=[]
    def panel(box,label,title):
        ax=fig.add_axes(box);fig.text(box[0]-.065,box[1]+box[3]+.018,label,fontsize=12,fontweight='bold')
        tt=fig.text(box[0],box[1]+box[3]+.018,title,fontsize=9,fontweight='bold');titles.append(tt)
        return ax
    ax=panel([.13,.625,.34,.245],'A','Independent-cohort QC retention')
    ax.barh(range(nt),s1a.QC_retention_fraction*100,color='#4D99A5',height=.65)
    for i,row in s1a.iterrows():ax.text(row.QC_retention_fraction*100+1.4,i,f'{int(row.QC_cells):,}',fontsize=8,va='center')
    ax.set_yticks(range(nt),[x.replace('_',' ') for x in TU]);ax.invert_yaxis();ax.set_xlim(0,122)
    ax.set_xticks([0,50,100]);ax.set_xlabel('QC retained / submitted cells (%)')
    ax.text(0,-.22,'Numbers show QC-retained cells.',transform=ax.transAxes,fontsize=8,va='top')
    ax=panel([.63,.625,.31,.245],'B','Within-donor comparison coverage')
    ax.set_xlim(-.5,2.5);ax.set_ylim(float(nt),-1.1);ax.axis('off')
    for j,h in enumerate(['NP n','Immune n','GEP2 pair']):ax.text(j,-.8,h,fontsize=8,fontweight='bold',ha='center',va='center')
    for i,row in s1b.iterrows():
        ax.add_patch(Rectangle((-.5,i-.45),3,.9,facecolor='#ECF4F3' if row.comparison_assessable else '#E5E5E5',edgecolor='white',lw=.5))
        ax.text(-.60,i,row.donor_id.replace('_',' '),fontsize=8,ha='right',va='center')
        for j,v in enumerate([f'{row.NP_supported_cells:,}',f'{row.immune_supported_cells:,}','Yes' if row.comparison_assessable else 'NA']):ax.text(j,i,v,fontsize=8,ha='center',va='center')
    missing_label=', '.join(d.replace('_',' ') for d in absent_comparators)
    ax.text(0,-.22,f'GEP2 NP > immune: {n_positive}/{n_comparable} assessable.\n'+(missing_label+': comparator absent.' if missing_label else 'All donor comparators available.'),transform=ax.transAxes,fontsize=8,va='top')
    ax=panel([.13,.195,.81,.285],'C','GEP5: matrix members (left) and stress members (right)')
    norm=colors.Normalize(vmin=0,vmax=float(s1c.mean_log1p_published_expression.max()))
    x={v:i for i,v in enumerate(sel)};y={v:i for i,v in enumerate(TU)}
    valid=s1c.mean_log1p_published_expression.notna()&s1c.fraction_detected.notna()
    d=s1c[valid&s1c.fraction_detected.gt(0)]
    im=ax.scatter(d.original_gene.map(x),d.donor_id.map(y),s=85*d.fraction_detected,c=d.mean_log1p_published_expression,cmap='viridis',norm=norm,linewidths=0,zorder=3)
    z=s1c[valid&s1c.fraction_detected.eq(0)];ax.scatter(z.original_gene.map(x),z.donor_id.map(y),s=10,facecolors='none',edgecolors='#888888',lw=.6,zorder=3)
    z=s1c[~valid];ax.scatter(z.original_gene.map(x),z.donor_id.map(y),s=22,c='#777777',marker='x',lw=.7,zorder=3)
    ax.set_xticks(range(len(sel)),sel,rotation=45,ha='right');ax.set_yticks(range(nt),[v.replace('_',' ') for v in TU])
    ax.set_xlim(-.6,len(sel)-.4);ax.set_ylim(nt-.4,-.6);ax.spines[:].set_visible(False);ax.tick_params(length=0)
    for i in range(nt):ax.axhline(i,lw=.5,color='#EFEFEF',zorder=0)
    ax.axvline(3.5,color='#AAAAAA',lw=.7,ls=':')
    cax=fig.add_axes([.16,.095,.29,.012]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label('Mean log1p(processed expression)',fontsize=8,labelpad=2);cb.ax.tick_params(labelsize=8,pad=1,length=2)
    handles=[plt.scatter([],[],s=85*v,c='#555555') for v in [.1,.5,1]]
    fig.legend(handles,['10%','50%','100%'],title='Detection (dot area)',loc='center left',bbox_to_anchor=(.52,.103),ncol=3,title_fontsize=8,handletextpad=.2,columnspacing=.8)
    fig.text(.10,.035,textwrap.fill('GEP5: '+g5_scope,95),fontsize=8,ha='left')
    fig.canvas.draw();ren=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height();outside=[];overlaps=[]
    texts=[v for v in fig.findobj(Text) if v.get_visible() and v.get_text().strip()]
    for tt in texts:
        bb=tt.get_window_extent(ren)
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1>w+1 or bb.y1>h+1:outside.append(tt.get_text())
    for title in titles:
        bb=title.get_window_extent(ren)
        for tt in texts:
            if tt is not title and bb.overlaps(tt.get_window_extent(ren)):overlaps.append([title.get_text(),tt.get_text()])
    stem=figs/'Supplementary_Figure_S1_Compact_Human_QC'
    fig.savefig(str(stem)+'.pdf');fig.savefig(str(stem)+'.svg');fig.savefig(str(stem)+'.png',dpi=300);plt.close(fig)
    pdf=fitz.open(str(stem)+'.pdf')
    for dpi,name in [(100,'100pct'),(150,'150pct')]:pdf[0].get_pixmap(matrix=fitz.Matrix(dpi/72,dpi/72),alpha=False).save(str(stem)+'_'+name+'.png')
    # S1 is a single tidy metric table rather than many duplicated small sheets.
    cols1=['cohort_GEO','sample_GSM','study_role','donor_or_scope','tissue_or_identity','measurement','value','unit','status_or_limit']
    r1=[]
    def add1(co,gsm,role,donor,identity,measure,value,unit,limit):r1.append(dict(zip(cols1,[co,gsm,role,donor,identity,measure,value,unit,limit])))
    for row in t.itertuples():
        common=('GSE165722',row.gsm,'Independent adult NP expression',row.donor_id)
        for identity,measure,value,unit,limit in [
          ('NP specimen','Submitted cells',int(row.submitted_cells),'cells','BD Rhapsody; published processed expression, not certified raw UMI'),
          ('All retained identities','QC-retained cells',int(row.QC_cells),'cells','≥200detected genes; positive sum; mitochondrial fraction≤0.20'),
          ('NP marker support','Operational NP-supported cells',int(ct.loc[row.donor_id,NP]),'cells','Marker co-detection; not histological purity'),
          ('Immune marker support','Immune-supported cells',int(ct.loc[row.donor_id,IM]),'cells','Zero cell count is a real absent comparator, not a zero score'),
          ('All retained identities','QC retention',float(row.QC_cells/row.submitted_cells),'fraction of submitted cells','Descriptive coverage; donor is observation unit'),
          ('NP versus immune','GEP2 comparison assessable',1 if comparable.loc[row.donor_id] else None,'assessable donor indicator',('NP > immune' if positive.loc[row.donor_id] else 'NP <= immune') if comparable.loc[row.donor_id] else 'NA: no immune-supported cells; score comparison unavailable'),
          ('NP specimen','Clinical grade crosswalk',None,'not analyzed',f'Reported GEO {row.reported_grade}; putative paper grade {row.putative_paper_grade}; crosswalk UNVERIFIED, no grade test')]:
            add1(*common,identity,measure,value,unit,limit)
    dp=donors[donors.cohort.eq('Swahn_GSE230808')].set_index('donor_id')
    dc=counts[counts.cohort.eq('Swahn_GSE230808')].pivot(index='donor_id',columns='identity',values='n_cells')
    at=tissue[tissue.cohort.eq('Swahn_GSE230808')].set_index('donor_id');pm=pairs.set_index('donor_id')
    for donor in DISC:
        for identity,measure,n,gsm in [('NP-source broad stromal','Retained NP-source stromal cells',int(dp.loc[donor,'QC_cells']),dp.loc[donor,'gsm']),
              ('NP-like','High-confidence NP-like cells',int(dc.loc[donor,'NP_chondrocyte_like']),dp.loc[donor,'gsm']),
              ('AF stromal','Retained paired AF stromal cells',int(at.loc[donor,'AF_retained_cells']),pm.loc[donor,'AF_GSM'])]:
            add1('GSE230808',gsm,'Adult disease discovery / paired AF-NP source',donor,identity,measure,n,'cells','10x 3prime v3.1; frozen identity/QC; cells are not independent patients')
        add1('GSE230808',dp.loc[donor,'gsm'],'Adult disease discovery / paired AF-NP source',donor,'Paired donor','Verified donor age',int(pm.loc[donor,'age']),'years',
             'Disease, age and series are confounded; no cross-series disease effect estimated')
    for row in inputs.itertuples():
        gsm=pm.loc[row.donor_id,'AF_GSM'] if row.group=='AF_stromal' else pm.loc[row.donor_id,'NP_GSM']
        add1('GSE230808',gsm,'Pooled new official communication input',row.donor_id,row.group,'Actually selected input cells',int(row.Freq),'cells',f'Within-donor quartile rule; cap 150/group; model input pooled across {n_pairs} donors')
    for row in elig[~elig.eligible].itertuples():
        for name,value in [('Lower quartile tail cells',row.n_low),('Upper quartile tail cells',row.n_high)]:
            add1('GSE230808',dp.loc[row.donor_id,'gsm'],'Excluded from communication input',row.donor_id,'High-confidence NP-like',name,int(value),'cells','Both tails below 20; retained NP count is not adequate quartile support')
    group_totals=inputs.groupby('group').Freq.sum()
    summary=[('GSE165722','Independent adult NP','Donor count',nt,'donors',f'{int(ct[NP].gt(0).sum())} donors have NP support'),
      ('GSE165722','NP versus immune','Assessable positive GEP2 comparisons',n_positive,f'of {n_comparable} assessable; {nt} donors total','Absent comparator(s): '+(', '.join(absent_comparators) or 'none')),
      ('GSE165722','Clinical grade','Grade analysis available',None,'not tested',grade_status+': GEO I–IV versus paper II–V crosswalk unresolved'),
      ('GSE230808','Source cohort','Disease NP donors with paired AF availability',nd,'paired donors','Linked GSE229711 is the same study system; no new healthy/disease test'),
      ('GSE230808','New official communication','Actually eligible paired donors',n_pairs,'donors','Pooled input does not imply independently trained/validated patient models'),
      ('GSE230808','New official communication','Selected pooled input',int(inputs.Freq.sum()),'cells','; '.join(f'{int(v)} {k}' for k,v in group_totals.items())),
      ('Human cohorts','Cohort design','Cross-series disease effect identifiable',None,'not estimated','Disease/age/series/platform confounding retained; no Harmony cure asserted')]
    for co,ident,measure,value,unit,limit in summary:add1(co,'series-level summary','Study scope','cohort',ident,measure,value,unit,limit)
    d1=pd.DataFrame(r1,columns=cols1);assert len(d1)<150
    d1.to_csv(tabs/'S1_Cohorts_QC.tsv',sep='\t',index=False,na_rep='NA')

    # S2 retains all important results in 135 tidy rows, not repeated 5program heatmaps.
    cols2=['entry_type','cohort_or_scope','donor','program','gene_or_context','frozen_rank','measurement','value','unit','availability','interpretation_or_NA_reason']
    r2=[]
    def add2(typ,co,donor,prog,gene,rank,measurement,value,unit,avail,note):r2.append(dict(zip(cols2,[typ,co,donor,prog,gene,rank,measurement,value,unit,avail,note])))
    for prog in PROGS:
        br=brief[brief.program.eq(prog)].iloc[0];lis=members[members.program.eq(prog)].sort_values('rank').gene.tolist()
        add2('Program overview','Frozen discovery K=5','all',prog,'; '.join(lis),None,'Frozen ranked membership',len(lis),'genes',f'Tu mapping {br.mapped_genes_min}–{br.mapped_genes_max}/{len(lis)}',
             br.expression_interpretation+': '+br.interpretation_scope+'; no refit or clinical-grade inference')
    g2=genes[genes.program.eq('GEP2')&genes.identity.eq(NP)]
    for row in members[members.program.eq('GEP2')].sort_values('rank').itertuples():
        df=g2[g2.original_gene.eq(row.gene)];n=df.mean_log1p_published_expression.notna().sum();miss=df.loc[df.mean_log1p_published_expression.isna(),'donor_id'].tolist()
        for metric,field,unit in [('Mean donor NP detection','fraction_detected','cell detection fraction; equal donor weight'),('Mean donor NP expression','mean_log1p_published_expression','log1p published processed expression; equal donor weight')]:
            v=df[field].dropna();note=(f'Available-donor range {v.min():.6g}–{v.max():.6g}; no zero imputation. ' if len(v) else 'No mapped measurement. ')
            note+=('Missing: '+','.join(miss)) if miss else f'All {nt} donors mapped.'
            released=sorted(set(df.mapped_gene.dropna().astype(str)))
            if released!=[row.gene]:note+=' Released feature(s): '+', '.join(released)+'. CCN2/CTGF is one gene.'
            add2('GEP2 frozen member','GSE165722','all available',row.program,row.gene,int(row.rank),metric,float(v.mean()) if len(v) else None,unit,f'{n}/{nt} mapped donors',note)
    for row in scores[scores.cohort.eq('Tu_GSE165722')&scores.identity.eq(NP)].itertuples():
        add2('NP donor score','GSE165722',row.donor_id,row.program,'NP marker support',None,'Donor median',float(row.score_median),'AUCell score; decoupler 2.2.0','observed',
            'UCell/AUCell magnitudes not comparable; grade analysis '+grade_status+'; NP identity operational')
    for donor in TU:
        row=tr.loc[donor];value=float(row.context_donor_median) if pd.notna(row.context_donor_median) else None
        add2('GEP2 comparator','GSE165722',donor,'GEP2','Immune marker support',None,'Donor median',value,'AUCell score; decoupler 2.2.0','observed' if value is not None else 'NA',
            (('NP > immune' if positive.loc[donor] else 'NP <= immune')+f'; one of {n_comparable} assessable pairs') if value is not None else 'No immune-supported cells; cannot compare; not a negative or zero score')
    discovery_support=support[support.cohort.eq('Swahn_GSE230808')]
    discovery_positive=int(discovery_support.NP_donor_median.gt(discovery_support.context_donor_median).sum())
    discovery_assessable=int(discovery_support.context_donor_median.notna().sum())
    for row in scores[scores.cohort.eq('Swahn_GSE230808')&scores.program.eq('GEP2')].itertuples():
        add2('Discovery GEP2 context','GSE230808',row.donor_id,'GEP2',row.identity,None,'Donor median',float(row.score_median),'UCell score; R UCell 2.4.0','verified cached score',
            f'Within-cohort context only; {discovery_positive}/{discovery_assessable} NP > other stromal and {discovery_assessable-discovery_positive} nonpositive; do not compare numeric magnitude with AUCell')
    for row in coverage[coverage.program.eq('GEP2')].itertuples():
        unavailable=g2.loc[g2.donor_id.eq(row.donor_id)&g2.mean_log1p_published_expression.isna(),'original_gene'].tolist()
        add2('GEP2 coverage','GSE165722',row.donor_id,'GEP2','Unchanged top25',None,'Mapped member count',int(row.n_mapped_genes),f'of {n_members} frozen genes','assessable',
            'Unavailable in this donor: '+(', '.join(unavailable) or 'none')+'; effective denominators documented in retained analysis')
    for gene,component in COMPONENTS.items():
        sub=s1c[s1c.original_gene.eq(gene)];v=sub.fraction_detected.dropna();fr=int(sub.frozen_rank.iloc[0]);n=len(v)
        add2('GEP5 component context','GSE165722','all available','GEP5',gene,fr,'Mean donor NP detection',float(v.mean()) if n else None,'cell detection fraction; equal donor weight',f'{n}/{nt} available donors',
            f'{component}; frozen illustrative component subset, not enrichment test; GEP5 '+g5_status+'; missing mapping differs from true zero')
    d2=pd.DataFrame(r2,columns=cols2);assert len(d2)<=150
    d2.to_csv(tabs/'S2_Fixed_NP_Programs.tsv',sep='\t',index=False,na_rep='NA')
    def sanitized(df):
        return df.astype(object).where(pd.notna(df),None).to_dict('records')
    payload={'source_status':'CACHED_INPUT_NOT_INDEPENDENTLY_CERTIFIED','only_compact_human_scope':True,
       'sheets':[{'name':'S1_Cohorts_QC','title':'Human cohort coverage and key QC','headers':cols1,'rows':[[r[c] for c in cols1] for r in sanitized(d1)],
            'types':['number' if c=='value' else 'text' for c in cols1],'source':'human/tables/S1_Cohorts_QC.tsv','source_sha256':digest(tabs/'S1_Cohorts_QC.tsv'),
            'description':'One tidy metric table; donor/tissue/sample identity retained; cells quantify coverage, not patient replication. NA is unavailable/unassessable, never a zero score.'},
           {'name':'S2_Fixed_NP_Programs','title':'Frozen programs and focused donor evidence','headers':cols2,'rows':[[r[c] for c in cols2] for r in sanitized(d2)],
            'types':['number' if c=='value' else 'integer' if c=='frozen_rank' else 'text' for c in cols2],'source':'human/tables/S2_Fixed_NP_Programs.tsv','source_sha256':digest(tabs/'S2_Fixed_NP_Programs.tsv'),
            'description':'Five frozen25member lists, full GEP2top25expression/detection, donor scores, coverage and mixedGEP5 context; UCell/AUCell separate scales.'}],
       'source_manifest':manifest,'no_formula_cells':True,'source_data_rows_are_static':True}
    write_json(out/'HUMAN_WORKBOOK_PAYLOAD.json',payload)
    dictionary=[]
    for sh,cols in [('S1_Cohorts_QC',cols1),('S2_Fixed_NP_Programs',cols2)]:
        for c in cols:dictionary.append({'sheet':sh,'field':c,'definition':{
            'value':'Static source-derived numeric count/score/fraction. Null means unavailable, not zero; lower and upper tail counts occupy distinct rows.',
            'sample_GSM':'Actual GEO sample accession; series-level or cross-cohort rows explicitly labelled.',
            'gene_or_context':'Original frozen member, identity comparator or complete ranked25member string in program-overview rows.',
            'frozen_rank':'Original rank1–25 when applicable; null for summaries.',
            'availability':'Observed, mapped donor count, verified cached score or explicitNA.',
            'unit':'Metric-specific denominator or score implementation; AUCell andUCell magnitudes not comparable.',
            'status_or_limit':'QC criterion, design boundary or precise unassessable reason; not a significance score.',
            'interpretation_or_NA_reason':'Retained positive, mixed, reversed or missing evidence; no composite score.'}.get(c,'Identifier, cohort role, measurement name or declared evidence context.'),
            'missing_rule':'NA/null is not measured zero. An observed zero immune cell count gives an unavailable comparator score.'})
    pd.DataFrame(dictionary).to_csv(out/'HUMAN_DATA_DICTIONARY.tsv',sep='\t',index=False)
    pd.DataFrame(manifest).to_csv(out/'HUMAN_INPUT_MANIFEST.tsv',sep='\t',index=False)
    g2cov=coverage[coverage.program.eq('GEP2')]
    stats={'Tu_donors':nt,'Tu_submitted_cells':int(t.submitted_cells.sum()),'Tu_QC_cells':int(t.QC_cells.sum()),'Tu_NP_supported_cells':int(ct[NP].sum()),'Tu_GEP2_comparable_donors':n_comparable,'Tu_GEP2_positive_comparisons':n_positive,'Tu_missing_comparator':absent_comparators,
        'frozen_K':len(PROGS),'frozen_members_per_program':n_members,'GEP2_mapped_range':[int(g2cov.n_mapped_genes.min()),int(g2cov.n_mapped_genes.max())],'GEP5_expression_status':g5_status,'grade_analysis':grade_status,
        'S1_rows':len(d1),'S2_rows':len(d2),'fig_panel_source_rows':{'A':len(s1a),'B':len(s1b),'C':len(s1c)}}
    write_json(out/'HUMAN_COMPACT_FACTS.json',stats)
    assert not d1.duplicated(['cohort_GEO','sample_GSM','study_role','donor_or_scope','tissue_or_identity','measurement']).any()
    assert not d2.duplicated(['entry_type','cohort_or_scope','donor','program','gene_or_context','measurement']).any()
    for sheet in payload['sheets']:
        assert len(sheet['headers'])==len(sheet['types'])
        assert all(len(row)==len(sheet['headers']) for row in sheet['rows'])
        for row in sheet['rows']:
            for val,typ in zip(row,sheet['types']):
                if val is not None and typ in ['integer','number']:assert isinstance(val,(int,float)),(sheet['name'],val,typ)
    qa={'status':'PASS_AUTOMATED_AWAITING_VISUAL_QA' if not outside and not overlaps else 'REVISE','canvas_text_overflow':outside,'panel_title_text_intersections':overlaps,
        'width_mm':180,'height_mm':210,'font_min_pt':8,'PDF_sha256':digest(str(stem)+'.pdf'),'export_paths':[str((figs/(stem.name+ext)).relative_to(out)).replace('\\','/') for ext in ['.pdf','.svg','.png']],
        'actual_visual_review':'PENDING'}
    write_json(out/'HUMAN_FIGURE_RENDER_QA.json',qa)
    write_json(out/'HUMAN_COMPACT_COMPLETION.json',{'status':'RENDERED_AWAITING_VISUAL_QA','facts':stats,'old_outputs_modified':False,'final_XLSX_created_by_this_agent':False})
    print(json.dumps({'stats':stats,'QA':qa},indent=2))

if __name__=='__main__':main()
