#!/usr/bin/env python3
"""Rebuild manuscript human figures from frozen, donor-level V3 tables.

Usage: python build_human_figures.py --source-root REVIEW --v4-root V4_CACHE --output-root OUTPUT
IVDD_REVIEW_ROOT and IVDD_V4_ROOT may replace the two input-root arguments.
Only frozen summary tables are read. Outputs must be outside input caches.
No cells, models, thresholds or programs are fitted. No manuscript text is emitted.
"""
import argparse
import hashlib
import json
import os
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
import pandas as pd
from PIL import Image

NP='NP_matrix_marker_supported'
ST='disc_stromal_marker_supported'
IM='immune_marker_supported'
EN='endothelial_marker_supported'
MI='mixed_lineage_marker_support'
UN='unresolved_low_marker_support'
IDS=[NP,ST,IM,EN,MI,UN]
IDLAB={NP:'NP support',ST:'Disc stromal',IM:'Immune',EN:'Endothelial',MI:'Mixed',UN:'Unresolved',
       'NP_chondrocyte_like':'NP-like','other_disc_stromal':'Other stromal'}
PAL={NP:'#168D9C',ST:'#DA9150',IM:'#8A6FA6',EN:'#6E98BD',MI:'#9A9A9A',UN:'#D9D9D9',
     'NP_chondrocyte_like':'#168D9C','other_disc_stromal':'#DA9150'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
    'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
    'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':2.5,'ytick.major.size':2.5,
    'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white','legend.frameon':False})


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source-root',type=Path,default=os.environ.get('IVDD_REVIEW_ROOT'))
    ap.add_argument('--v4-root',type=Path,default=os.environ.get('IVDD_V4_ROOT'),help='Read-only authoring cache with eligibility and Figure1A design source')
    ap.add_argument('--output-root',type=Path,required=True)
    ap.add_argument('--eligibility-table',type=Path,help='Override V4_CACHE/00_inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv')
    ap.add_argument('--design-table',type=Path,help='Override V4_CACHE/03_tables/human/Figure1A_source.tsv (layer, label, scope)')
    a=ap.parse_args()
    if a.source_root is None:ap.error('Provide --source-root or IVDD_REVIEW_ROOT; no frozen result tables are bundled.')
    if a.v4_root is None and (a.eligibility_table is None or a.design_table is None):
        ap.error('Provide --v4-root/IVDD_V4_ROOT or both --eligibility-table and --design-table.')
    src=a.source_root.resolve();out=a.output_root.resolve()
    for root in [src]+([a.v4_root.resolve()] if a.v4_root else []):
        if out==root or root in out.parents:ap.error('--output-root must be outside the read-only input caches.')
    eligpath=a.eligibility_table or a.v4_root/'00_inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv'
    designpath=a.design_table or a.v4_root/'03_tables/human/Figure1A_source.tsv'
    for path in [src,eligpath,designpath]:
        if not path.exists():ap.error('Required input is absent: '+path.name+' (see README cache contract).')
    tab=out/'03_tables/human';qa=out/'08_review/human'
    for d in [tab,qa,out/'04_figures/Figure1',out/'04_figures/Figure2',out/'06_supplementary/figures']:
        d.mkdir(parents=True,exist_ok=True)
    maps=[];facts=[];manifest=[];qarecs=[]
    def read(rel):
        p=src/rel;f=pd.read_csv(p,sep='\t')
        manifest.append({'relative_source':rel,'sha256':sha(p),'bytes':p.stat().st_size,
                         'role':'frozen donor-level computational source','submission_status':'CACHED_INPUT_NOT_INDEPENDENTLY_CERTIFIED'})
        return f
    h='01_human_programs/'
    cohort=read(h+'COHORT_AND_DONOR_MANIFEST.tsv')
    cnt=read(h+'FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv')
    score=read(h+'FIGURE2_DONOR_PROGRAM_SOURCE.tsv')
    frozen=read(h+'FROZEN_K5_PROGRAMS.tsv')
    coverage=read(h+'PROGRAM_GENE_COVERAGE.tsv')
    mapping=read(h+'PROGRAM_GENE_MAPPING.tsv')
    effective=read(h+'AUCELL_EFFECTIVE_PROGRAM_COVERAGE.tsv')
    m=read(h+'TU_MARKER_DOTPLOT_SOURCE.tsv')
    g=read(h+'TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv')
    cand=read(h+'TU_CANDIDATE_GENE_BY_DONOR_IDENTITY.tsv')
    brief=read(h+'ALL_PROGRAMS_BRIEF.tsv')
    support=read(h+'GEP2_DONOR_SUPPORT.tsv')
    qc=read(h+'TU_MATRIX_AND_QC_MANIFEST.tsv')
    tissue=read('05_figures/Figure1_Figure2/Figure1a_donor_tissue_counts.tsv')
    comm=read('02_sigxtalk_official_new/formal_new/tables/DONOR_GROUP_INPUT_COUNTS.tsv')
    elig=pd.read_csv(eligpath,sep='\t')
    design=pd.read_csv(designpath,sep='\t').sort_values('layer')
    assert list(design.columns)==['layer','label','scope'] and len(design)==4,'Expected frozen four-layer study-design source.'
    manifest.append({'relative_source':'V4/03_tables/human/'+designpath.name,'sha256':sha(designpath),'bytes':designpath.stat().st_size,
        'role':'frozen schematic labels and counts; no experiment values','submission_status':'CACHED_INPUT_NOT_INDEPENDENTLY_CERTIFIED'})
    manifest.append({'relative_source':'V4 inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv','sha256':sha(eligpath),
        'bytes':eligpath.stat().st_size,'role':'actual frozen donor eligibility recovered read-only from V3 runtime','submission_status':'CACHED_INPUT_NOT_INDEPENDENTLY_CERTIFIED'})
    elig.to_csv(tab/'S2_COMMUNICATION_DONOR_ELIGIBILITY.tsv',sep='\t',index=False)
    for rel in [h+'AUCELL_RUN_METADATA.json',h+'PROVENANCE_AND_SCORE_SCALE.md',
                '07_public_repo/resources/frozen_k5/K_selection_and_parameters.tsv',
                '07_public_repo/resources/frozen_k5/PROGRAM_LOCK.json']:
        p=src/rel
        manifest.append({'relative_source':rel,'sha256':sha(p),'bytes':p.stat().st_size,
            'role':'executed method parameters/provenance; no re-estimation','submission_status':'CACHED_INPUT_NOT_INDEPENDENTLY_CERTIFIED'})
    kparam=read('07_public_repo/resources/frozen_k5/K_selection_and_parameters.tsv')
    kparam.to_csv(tab/'S3_CNMf_FROZEN_K_SELECTION_PARAMETERS.tsv',sep='\t',index=False)
    P=sorted(frozen.program.unique());n_members=frozen.groupby('program').size()
    assert len(P)==5 and n_members.eq(25).all(),'This builder implements the locked K=5/top25 analysis contract.'
    TUD=cohort.loc[cohort.cohort.eq('Tu_GSE165722'),'donor_id'].drop_duplicates().sort_values().tolist()
    DONS=cohort.loc[cohort.cohort.eq('Swahn_GSE230808'),'donor_id'].drop_duplicates().sort_values().tolist()
    nt,nd=len(TUD),len(DONS);nc=int(comm.donor_id.nunique())
    assert nt and nd and set(comm.donor_id).issubset(DONS)
    assert set(score.score_method)=={'UCell_2.4.0','AUCell_decoupler_2.2.0'}
    assert cnt.loc[cnt.cohort.eq('Tu_GSE165722'),'n_cells'].sum()==qc.n_qc_pass.sum()
    assert cnt.loc[cnt.identity.eq(NP),'n_cells'].sum()<=qc.n_qc_pass.sum()
    assert not cand.duplicated(['gsm','identity','gene']).any()
    for df,name in [(cohort,'S1_PUBLIC_COHORTS'),(cnt,'S2_DONOR_IDENTITY_COUNTS'),(qc,'S2_TU_QC'),
            (frozen,'S3_ALL125_FROZEN_PROGRAM_MEMBERS'),(coverage,'S3_PROGRAM_COVERAGE'),
            (mapping,'S3_GENE_MAPPING'),(effective,'S3_AUCELL_EFFECTIVE_DENOMINATORS'),
            (score,'S4_DONOR_PROGRAM_SCORES'),(support,'S4_GEP2_DONOR_SUPPORT'),
            (brief,'S4_PROGRAM_INTERPRETATION'),(g,'S4_ALL_PROGRAM_GENE_EXPRESSION'),
            (m,'S2_COMPLETE_MARKER_EXPRESSION'),(cand,'S4_CANDIDATE_EXPRESSION')]:
        df.to_csv(tab/(name+'.tsv'),sep='\t',index=False,na_rep='NA')
    def save_source(df,panel,source,rule):
        p=tab/(panel+'_source.tsv');df.to_csv(p,sep='\t',index=False,na_rep='NA')
        maps.append({'panel_id':panel,'input_file':source,'variables':'see source columns',
            'transformation':rule,'statistical_unit':'donor unless separately stated',
            'independent_n':('not applicable; analysis schematic' if panel=='Figure1A' else f'{nd} discovery and {nt} independent donors; communication input from {nc} paired donors'),
            'missing_reason':('Absent context cells yield an unavailable comparator median, not a zero score' if panel=='Figure2E' else
                              'Unavailable mappings recorded in the source table; no zero imputation' if panel in ['Figure2B','Figure2C'] else
                              'Unmapped genes/absent identities remain NA; actual zero detection retained separately'),
            'statistical_test':'none; descriptive','source_data':p.relative_to(out).as_posix(),
            'source_rows':len(df),'source_sha256':sha(p)})
        return df
    def fact(fid,value,unit,source,rule,fig):
        facts.append({'fact_id':fid,'value':value,'unit':unit,'source_file':source,
            'row_or_rule':rule,'calculation_script':'scripts/authoring/v4/human/build_human_figures.py',
            'independent_n':f'{nd} discovery and {nt} independent donors; communication subset {nc} donors',
            'figure':fig,'table':'Supplementary Tables S1-S4','scope':'expression and identity context; no clinical-grade inference'})
    fact('HUMAN_TU_SUBMITTED',int(qc.n_cells.sum()),'cells',h+'TU_MATRIX_AND_QC_MANIFEST.tsv','sum(n_cells)','FigureS1A')
    fact('HUMAN_TU_QC',int(qc.n_qc_pass.sum()),'cells',h+'TU_MATRIX_AND_QC_MANIFEST.tsv','sum(n_qc_pass)','Figure1C')
    fact('HUMAN_TU_NP',int(cnt.loc[cnt.identity.eq(NP),'n_cells'].sum()),'cells',h+'FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv','identity=NP_matrix_marker_supported; sum(n_cells)','Figure1C')
    fact('HUMAN_DISCOVERY_STROMAL',int(cnt.loc[cnt.cohort.eq('Swahn_GSE230808'),'n_cells'].sum()),'cells',h+'FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv','cohort=Swahn_GSE230808; sum(n_cells)','Figure1E')
    fact('HUMAN_DISCOVERY_NP',int(cnt.loc[cnt.identity.eq('NP_chondrocyte_like'),'n_cells'].sum()),'cells',h+'FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv','identity=NP_chondrocyte_like; sum(n_cells)','Figure1B')
    tu_support=support[support.cohort.eq('Tu_GSE165722')].copy()
    paired=tu_support[tu_support.context_donor_median.notna()]
    n_positive=int((paired.NP_donor_median>paired.context_donor_median).sum())
    missing_comparators=tu_support.loc[tu_support.context_donor_median.isna(),'donor_id'].tolist()
    fact('GEP2_TU_NP_GT_IMMUNE',n_positive,f'of {len(paired)} comparable donors; {nt} total',h+'GEP2_DONOR_SUPPORT.tsv','finite context median and NP>context; missing contexts remain NA','Figure2E')
    fact('GEP2_TU_COVERAGE_MIN',int(coverage[coverage.program.eq('GEP2')].n_mapped_genes.min()),'of 25 frozen genes',h+'PROGRAM_GENE_COVERAGE.tsv','program=GEP2; min(n_mapped_genes)','Figure2B')
    fact('GEP2_TU_COVERAGE_MAX',int(coverage[coverage.program.eq('GEP2')].n_mapped_genes.max()),'of 25 frozen genes',h+'PROGRAM_GENE_COVERAGE.tsv','program=GEP2; max(n_mapped_genes)','Figure2B')
    for co,sub in cohort.groupby('cohort'):
        fact('HUMAN_DONORS_'+co,int(sub.donor_id.nunique()),'donors',h+'COHORT_AND_DONOR_MANIFEST.tsv','nunique(donor_id) for '+co,'Figure1')
    fact('COMMUNICATION_ELIGIBLE_DONORS',int(comm.donor_id.nunique()),'paired donors','DONOR_GROUP_INPUT_COUNTS.tsv','nunique(donor_id)','Figure1B')
    for group in ['AF_stromal','NP_GEP2_high','NP_GEP2_low']:
        fact('COMMUNICATION_'+group.upper(),int(comm.loc[comm.group.eq(group),'Freq'].sum()),'selected cells','DONOR_GROUP_INPUT_COUNTS.tsv','sum(Freq) for '+group,'Figure1B')
    for d in elig.loc[~elig.eligible,'donor_id']:
        er=elig[elig.donor_id.eq(d)].iloc[0]
        for field in ['n_NP','n_low','n_high']:
            fact('EXCLUDED_'+d+'_'+field,int(er[field]),'cells','RECOMPUTED_DONOR_ELIGIBILITY.tsv','donor_id='+d+'; '+field,'Figure1B')
    fact('FROZEN_PROGRAM_NUMBER',int(frozen.program.nunique()),'programs',h+'FROZEN_K5_PROGRAMS.tsv','nunique(program)','Figure2A')
    fact('FROZEN_MEMBERS_PER_PROGRAM',int(n_members.iloc[0]),'members per program',h+'FROZEN_K5_PROGRAMS.tsv','groupby(program).size()','Figure2C')
    dsupp=support[support.cohort.eq('Swahn_GSE230808')]
    fact('GEP2_DISCOVERY_NP_GT_STROMAL',int((dsupp.NP_donor_median>dsupp.context_donor_median).sum()),f'of {nd} donors',h+'GEP2_DONOR_SUPPORT.tsv','Swahn; NPmedian>contextmedian','Figure2D')
    for gene in ['PLOD2','CCN2','PRSS23']:
        vals=g[g.program.eq('GEP2')&g.identity.eq(NP)&g.original_gene.eq(gene)].fraction_detected.dropna()
        fact(gene+'_TU_NP_DETECTION_MIN',float(vals.min()),'detection fraction',h+'TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv','GEP2; NP support; original_gene='+gene+'; minimum acrossdonors','Figure2C')
        fact(gene+'_TU_NP_DETECTION_MAX',float(vals.max()),'detection fraction',h+'TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv','GEP2; NP support; original_gene='+gene+'; maximum acrossdonors','Figure2C')
    for pid,sub in score[score.identity.eq(NP)].groupby('program'):
        fact(pid+'_TU_NP_SCORE_MIN',sub.score_median.min(),'AUCell',h+'FIGURE2_DONOR_PROGRAM_SOURCE.tsv','NP_matrix_marker_supported; min donor median','Figure2A')
        fact(pid+'_TU_NP_SCORE_MAX',sub.score_median.max(),'AUCell',h+'FIGURE2_DONOR_PROGRAM_SOURCE.tsv','NP_matrix_marker_supported; max donor median','Figure2A')

    # Figure contracts are saved before drawing.
    contracts={
      'Figure1':{'claim':'Adult disc cohorts provide distinct donor-resolved tissue and marker-support contexts for fixed-program analysis.',
        'archetype':'schematic-led composite','hero':'D independent identity-marker dotplot',
        'panels':{'A':'four-layer data-to-question schematic','B':'actual paired discovery availability and communication inclusion',
          'C':'Tu donor identity composition','D':'independent canonical markers','E':'Swahn donor identity counts; permitted matrix substitute for unavailable cached embeddings'},
        'boundaries':['no invented UMAP','cohort and disease/platform effects are confounded','marker support is operational identity','no inferred causal transfer arrows']},
      'Figure2':{'claim':f'The unchanged GEP2 gene set has {n_positive} of {len(paired)} assessable NP-versus-immune donor comparisons in the positive direction. All program contexts remain visible.',
        'archetype':'quantitative grid','hero':f'C all{n_members.GEP2} frozen GEP2 genes in all{nt} Tu donors',
        'panels':{'A':'all5programs, separately scaled UCell and AUCell','B':'actual GEP2 effective coverage and missing members','C':'25member gene expression/detection bubble matrix',
          'D':'discovery donor NP/other-stromal medians','E':'Tu pairedNP/immune with absent comparators explicitly NA','F':'all5program support and discordance'},
        'boundaries':['no cross-method score-magnitude comparison','no grade inference','no refit or gene substitution','missing mappings not zero expression']},
      'FigureS1':{'claim':'The full donor QC, composition and independent marker panels expose coverage and identity limitations.',
        'archetype':'quantitative grid','hero':'C-F complete pre-existing canonical marker sets',
        'panels':{'A':'allTu submitted/QC/NP counts','B':'allSwahn tissue count coverage','C':'NPmatrix+discECM markers','D':'immune+myeloid markers','E':'endothelial markers','F':'contractile markers'},
        'boundaries':['processed expression is not certified raw UMI','marker availability differs','small identities and unresolved labels retained']},
      'FigureS2':{'claim':'All five frozen 25-gene programs show uneven identity context and mapping coverage; mixed programs remain visible.',
        'archetype':'quantitative grid','hero':'B-F complete125members by all6 Tu identities',
        'panels':{'A':f'{len(P)} programs × {nt} donor gene coverage','B':'GEP1 all25members','C':'GEP2 all25members','D':'GEP3 all25members','E':'GEP4 all25members','F':'GEP5 all25members'},
        'boundaries':['equal donor weighting','missing mappings explicit','no specificity claim for overlapping annotation markers','Frozen GEP5 expression review: '+brief.set_index('program').loc['GEP5','expression_interpretation']]}}
    for k,v in contracts.items():
        v.update({'backend':'Python matplotlib','width_mm':180,'height_mm':240 if k=='FigureS2' else 230,
           'min_text_pt':8,'exports':['PDF','SVG','PNG300dpi','TIFF600dpi'],
           'statistics':'Donor summaries only. No inferential P values, no cells-as-patients tests.'})
    (qa/'FIGURE_CONTRACTS.json').write_text(json.dumps(contracts,indent=2),encoding='utf8')

    def fig_new(height=230):
        return plt.figure(figsize=(180/25.4,height/25.4))
    def panel(fig,box,letter,title):
        ax=fig.add_axes(box)
        fig.text(max(.012,box[0]-.055),box[1]+box[3]+.018,letter,fontsize=12,fontweight='bold',va='bottom')
        artist=fig.text(box[0],box[1]+box[3]+.018,title,fontsize=9,fontweight='bold',va='bottom')
        artist._qa_role='panel_title'
        return ax
    def heat(ax,df,xcol,ycol,value,xs,ys,cmap='viridis',vmin=0,vmax=None,annotate=False,fmt='.2f'):
        ar=df.pivot(index=ycol,columns=xcol,values=value).reindex(index=ys,columns=xs).to_numpy(float)
        cm=plt.get_cmap(cmap).copy();cm.set_bad('#E0E0E0')
        im=ax.imshow(ar,aspect='auto',cmap=cm,vmin=vmin,vmax=vmax,interpolation='nearest')
        ax.set_xticks(range(len(xs)),xs);ax.set_yticks(range(len(ys)),ys)
        ax.tick_params(length=0,pad=3);ax.spines[:].set_visible(False)
        if annotate:
            for i in range(len(ys)):
                for j in range(len(xs)):
                    txt='NA' if not np.isfinite(ar[i,j]) else format(ar[i,j],fmt)
                    ax.text(j,i,txt,ha='center',va='center',fontsize=8,
                        color='white' if np.isfinite(ar[i,j]) and ar[i,j]>(vmax or np.nanmax(ar))*.58 else '#222222')
        return im
    def cbar(fig,im,box,label,ticks=None):
        cax=fig.add_axes(box);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=ticks)
        cb.set_label(label,fontsize=8,labelpad=2);cb.ax.tick_params(labelsize=8,pad=1,length=2)
        return cb
    def bubble(ax,df,xcol,ycol,xs,ys,vmax=None,fs=8,area=42):
        df=df.copy();xx={x:i for i,x in enumerate(xs)};yy={y:i for i,y in enumerate(ys)}
        df=df[df[xcol].isin(xs)&df[ycol].isin(ys)]
        present=df.mean_log1p_published_expression.notna() & df.fraction_detected.notna()
        norm=colors.Normalize(0,vmax or df.mean_log1p_published_expression.max())
        ok=df[present & df.fraction_detected.gt(0)]
        im=ax.scatter(ok[xcol].map(xx),ok[ycol].map(yy),s=area*ok.fraction_detected,
            c=ok.mean_log1p_published_expression,cmap='viridis',norm=norm,linewidths=0,zorder=3)
        zero=df[present & df.fraction_detected.eq(0)]
        ax.scatter(zero[xcol].map(xx),zero[ycol].map(yy),s=5,facecolors='none',edgecolors='#888888',linewidths=.5,zorder=3)
        missing=df[~present]
        ax.scatter(missing[xcol].map(xx),missing[ycol].map(yy),s=14,c='#888888',marker='x',linewidths=.6,zorder=3)
        ax.set_xlim(-.6,len(xs)-.4);ax.set_ylim(len(ys)-.5,-.5)
        ax.set_xticks(range(len(xs)),xs,rotation=65,ha='right',fontsize=fs)
        ax.set_yticks(range(len(ys)),[IDLAB.get(y,y.replace('Tu_','Tu ')) for y in ys],fontsize=fs)
        ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)
        for y in range(len(ys)):
            ax.axhline(y,color='#F0F0F0',lw=.5,zorder=0)
        return im
    def dot_legend(fig,x,y,area=42):
        handles=[plt.scatter([],[],s=area*p,c='#555555') for p in [.1,.5,1]]
        fig.legend(handles,['10%','50%','100%'],loc='center left',bbox_to_anchor=(x,y),ncol=3,
            frameon=False,handletextpad=.2,columnspacing=.8,title='Detection (dot area)',title_fontsize=8)
    def save(fig,name,directory):
        fig.canvas.draw();renderer=fig.canvas.get_renderer();W,H=fig.canvas.get_width_height()
        outside=[];container_overflow=[];title_overlap=[];texts=[]
        for artist in fig.findobj(Text):
            if not artist.get_visible() or not artist.get_text().strip():continue
            bb=artist.get_window_extent(renderer)
            texts.append((artist,bb))
            if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W+1 or bb.y1 > H+1:
                outside.append({'text':artist.get_text(),'bbox':list(map(float,bb.extents))})
            if hasattr(artist,'_qa_container'):
                bb0=artist.axes.transAxes.transform_bbox(matplotlib.transforms.Bbox.from_extents(*artist._qa_container))
                if bb.x0<bb0.x0-1 or bb.x1>bb0.x1+1 or bb.y0<bb0.y0-1 or bb.y1>bb0.y1+1:
                    container_overflow.append({'text':artist.get_text(),'bbox':list(map(float,bb.extents)),'container':list(map(float,bb0.extents))})
        for artist,bb in texts:
            if getattr(artist,'_qa_role','')!='panel_title':continue
            for other,bb2 in texts:
                if artist is other:continue
                if bb.overlaps(bb2):
                    title_overlap.append({'title':artist.get_text(),'overlapping_text':other.get_text()})
        p=directory/name
        fig.savefig(str(p)+'.pdf');fig.savefig(str(p)+'.svg');fig.savefig(str(p)+'.png',dpi=300)
        fig.savefig(str(p)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        # Render the actual exported PDF in Python; keep independent TIFF check.
        import fitz
        doc=fitz.open(str(p)+'.pdf');page=doc[0]
        for scale in [1,2]:
            px=page.get_pixmap(matrix=fitz.Matrix(100/72*scale,100/72*scale),alpha=False)
            px.save(str(p)+('_PDF_100pct.png' if scale==1 else '_PDF_200pct.png'))
        text=page.get_text()
        ti=Image.open(str(p)+'.tiff')
        ti.thumbnail((1200,1600));ti.convert('RGB').save(str(p)+'_TIFF_check.png')
        qarecs.append({'figure':name,'status':'PASS_AUTOMATED' if not (outside or container_overflow or title_overlap) else 'REVISE',
            'width_mm':round(page.rect.width/72*25.4,3),'height_mm':round(page.rect.height/72*25.4,3),
            'outside_text':outside,'editable_pdf_text_characters':len(text),
            'panel_container_overflow':container_overflow,'panel_title_text_overlaps':title_overlap,
            'png_size':list(Image.open(str(p)+'.png').size),'tiff_dpi':list(map(float,Image.open(str(p)+'.tiff').info.get('dpi',()))),
            'min_font_pt':8,'visual_status':'PENDING_ACTUAL_INSPECTION'})
        plt.close(fig)

    # ---- Figure 1 ------------------------------------------------------
    fig=fig_new()
    ax=panel(fig,[.075,.805,.88,.14],'A','Study design and evidence layers');ax.axis('off')
    stages=[(r.label,r.scope.replace('; ','\n')) for r in design.itertuples()]
    for j,(title,body) in enumerate(stages):
        x=j*.253
        ax.add_patch(Rectangle((x,.14),.232,.73,facecolor=['#E9F3F5','#F1EFF7','#F7F0E7','#F1F1F1'][j],edgecolor='#CCCCCC',lw=.6))
        title_art=ax.text(x+.012,.82,textwrap.fill(title,18),fontsize=8,fontweight='bold',va='top',linespacing=1.15)
        title_art._qa_container=(x,.14,x+.232,.87)
        wrapped='\n'.join(textwrap.fill(line,25) for line in body.split('\n'))
        body_art=ax.text(x+.012,.53,wrapped,fontsize=8,va='top',linespacing=1.2)
        body_art._qa_container=(x,.14,x+.232,.87)
        if j<3:ax.annotate('',xy=(x+.25,.50),xytext=(x+.233,.50),arrowprops={'arrowstyle':'->','lw':.8,'color':'#555555'})
    save_source(design,'Figure1A','study design and executed source metadata','conceptual analysis relationships, not causal molecular arrows')
    sel=comm.pivot(index='donor_id',columns='group',values='Freq').reindex(DONS)
    ds=cnt[cnt.cohort.eq('Swahn_GSE230808')].pivot(index='donor_id',columns='identity',values='n_cells').reindex(DONS)
    avail=tissue[tissue.cohort.eq('Swahn_GSE230808')].set_index('donor_id').reindex(DONS)
    bdata=pd.DataFrame({'donor_id':DONS,'AF_retained_cells':avail.AF_retained_cells.astype(int).to_numpy(),
        'NP_like_cells':ds.NP_chondrocyte_like.astype(int).to_numpy(),
        'NP_high_input':sel.NP_GEP2_high.fillna(0).astype(int).to_numpy(),
        'NP_low_input':sel.NP_GEP2_low.fillna(0).astype(int).to_numpy(),
        'communication_included':sel.AF_stromal.notna().to_numpy()})
    bdata=bdata.merge(elig[['donor_id','n_low','n_high','q25','q75','eligible']],on='donor_id',validate='one_to_one')
    bdata['exclusion_reason']=np.where(bdata.communication_included,'included',
        'both NP quartile tails below min20; actual n_low/n_high in this table')
    save_source(bdata,'Figure1B','Figure1a_donor_tissue_counts.tsv; DONOR_GROUP_INPUT_COUNTS.tsv; donor identity counts',
        'actual retained AF/NP-like counts; high/low input displayed as dash when not included; original fixed selection rule')
    ax=panel(fig,[.16,.555,.33,.205],'B','Paired discovery donors')
    for i,row in bdata.iterrows():
        ax.add_patch(Rectangle((-.5,i-.45),4,.9,facecolor='#E9F3F5' if row.communication_included else '#F1F1F1',edgecolor='white'))
        vals=[str(row.AF_retained_cells),str(row.NP_like_cells),str(row.NP_high_input) if row.communication_included else '—',str(row.NP_low_input) if row.communication_included else '—']
        for j,v in enumerate(vals):ax.text(j,i,v,ha='center',va='center',fontsize=8)
    ax.set_xlim(-.5,3.5);ax.set_ylim(nd-.5,-.5);ax.set_xticks(range(4),['AF','NP-like','High¹','Low¹'])
    ax.set_yticks(range(nd),DONS);ax.tick_params(length=0);ax.spines[:].set_visible(False)
    excluded_tails=' or '.join(f'{int(r.n_low)}/{int(r.n_high)}' for r in elig[~elig.eligible].itertuples())
    ax.text(0,-.16,f'¹ Selected cells; blue: {nc} included pairs.\nGrey: NP tails {excluded_tails}; minimum 20 each.',transform=ax.transAxes,fontsize=8,va='top')
    tu=cnt[cnt.cohort.eq('Tu_GSE165722')].pivot(index='donor_id',columns='identity',values='n_cells').reindex(index=TUD,columns=IDS).fillna(0)
    cdata=tu.reset_index().melt(id_vars='donor_id',var_name='identity',value_name='n_cells')
    cdata['donor_total']=cdata.donor_id.map(tu.sum(axis=1));cdata['fraction']=cdata.n_cells/cdata.donor_total
    save_source(cdata,'Figure1C',h+'FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv','within-donor fractions of all QC cells; no pooling across donors')
    ax=panel(fig,[.62,.555,.335,.205],'C','Independent NP-source cohort')
    left=np.zeros(nt)
    for identity in IDS:
        frac=tu[identity]/tu.sum(axis=1);ax.barh(range(nt),frac,left=left,color=PAL[identity],height=.7,label=IDLAB[identity]);left+=frac.to_numpy()
    ax.set_yticks(range(nt),[x.replace('_',' ') for x in TUD]);ax.set_ylim(nt-.4,-.6);ax.set_xlim(0,1)
    ax.set_xticks([0,.5,1],['0','50','100']);ax.set_xlabel('QC cell composition (%)')
    ax.legend(loc='upper left',bbox_to_anchor=(-.10,-.30),ncol=2,handlelength=.9,columnspacing=.7,labelspacing=.35,borderaxespad=0)
    md=m.groupby(['identity','gene'],as_index=False).agg(fraction_detected=('fraction_detected','mean'),
        mean_log1p_published_expression=('mean_log1p_published_expression','mean'),n_contributing_donors=('donor_id','nunique'))
    markers=['ACAN','COL2A1','SOX9','KRT19','COL1A1','DCN','PTPRC','LST1','CD3D','PECAM1','VWF']
    ddata=md[md.gene.isin(markers)].copy()
    save_source(ddata,'Figure1D',h+'TU_MARKER_DOTPLOT_SOURCE.tsv','equal mean across available donor/identity summaries; pre-existing marker panels, independent of GEP2 selection')
    ax=panel(fig,[.20,.245,.72,.145],'D','Canonical identity-marker context')
    im=bubble(ax,ddata,'gene','identity',markers,IDS,area=65)
    cbar(fig,im,[.23,.150,.28,.012],'Mean donor log1p(expression)')
    dot_legend(fig,.56,.157,area=65)
    edata=ds.reset_index().melt(id_vars='donor_id',var_name='identity',value_name='n_cells')
    save_source(edata,'Figure1E',h+'FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv',
        'all discovery stromal identities retained; matrix substitute for unavailable cached cell-level embedding')
    # Full source-count matrix is given directly as a small final panel.
    ax=panel(fig,[.17,.035,.76,.040],'E','Discovery identity composition (donor order as B)')
    vals=ds[['NP_chondrocyte_like','other_disc_stromal']].T
    im=ax.imshow(vals.to_numpy(),aspect='auto',cmap='Blues',norm=colors.LogNorm(vmin=1,vmax=float(vals.max().max())))
    ax.set_xticks(range(nd),[d.replace('SP21.','').replace('SP','') for d in DONS])
    ax.set_yticks(range(2),['NP-like','Other stromal']);ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)
    for i in range(2):
        for j in range(nd):
            v=int(vals.iloc[i,j]);ax.text(j,i,f'{v:,}',ha='center',va='center',fontsize=8,color='white' if v>1200 else '#222222')
    save(fig,'Figure1_Cohort_and_Cellular_Framework',out/'04_figures/Figure1')

    # ---- Figure 2 ------------------------------------------------------
    fig=fig_new()
    a_data=score[score.identity.isin([NP,'NP_chondrocyte_like'])].copy()
    save_source(a_data,'Figure2A',h+'FIGURE2_DONOR_PROGRAM_SOURCE.tsv','donor median; only NP-support identity; independent method-specific colour scales')
    for idx,(co,donors,title,box,vmax) in enumerate([
      ('Swahn_GSE230808',DONS,'Discovery NP-like · UCell',[.13,.835,.34,.095],.8),
      ('Tu_GSE165722',TUD,'Independent NP support · AUCell',[.61,.835,.34,.095],.45)]):
        ax=panel(fig,box,'A' if idx==0 else '',title)
        sub=a_data[a_data.cohort.eq(co)];im=heat(ax,sub,'donor_id','program','score_median',donors,P,vmax=vmax)
        ax.set_xticklabels([d.split('.')[-1] if co.startswith('Swahn') else d.replace('Tu_','') for d in donors])
        cbar(fig,im,[box[0]+.035,.790,.255,.010],title.split('·')[-1].strip()+' score')
    bdata=coverage[coverage.program.eq('GEP2')].copy()
    mapping['donor_id']=mapping.gsm.map(cohort.set_index('gsm').donor_id)
    missing=mapping[(mapping.program.eq('GEP2')) & mapping.mapped_gene.isna()]
    missing_by_d=missing.groupby('donor_id').original_gene.agg(';'.join)
    bdata['unavailable_frozen_genes']=bdata.donor_id.map(missing_by_d).fillna('none')
    bdata=bdata.set_index('donor_id').reindex(TUD).reset_index()
    save_source(bdata,'Figure2B',h+'PROGRAM_GENE_COVERAGE.tsv; '+h+'PROGRAM_GENE_MAPPING.tsv','mapped members/25; unavailable genes are not zeros')
    ax=panel(fig,[.14,.675,.80,.040],'B','GEP2 mapped coverage (25 fixed members)')
    ax.bar(range(nt),bdata.n_mapped_genes,color='#168D9C',width=.62)
    ax.set_ylim(0,29);ax.set_yticks([0,25]);ax.set_xticks(range(nt),[d.replace('_',' ') for d in TUD]);ax.set_ylabel('Genes')
    for j,v in enumerate(bdata.n_mapped_genes):ax.text(j,v+.5,str(v),ha='center',va='bottom',fontsize=8)
    genes=frozen[frozen.program.eq('GEP2')].sort_values('rank').gene.tolist()
    cdata=g[g.program.eq('GEP2') & g.identity.eq(NP)].copy()
    cdata['frozen_rank']=cdata.original_gene.map(dict(zip(genes,range(1,len(genes)+1))))
    cdata['mapping_state']=np.where(cdata.gene_available.astype(str).str.lower().eq('true'),'mapped','unavailable')
    assert len(cdata)==len(genes)*nt
    save_source(cdata,'Figure2C',h+'TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv','all25genes in frozen rank order; actual donor NP-support means and detection; grey x for unavailable mapping')
    ax=panel(fig,[.14,.465,.80,.145],'C','All 25 GEP2 members in independent NP-supported cells')
    im=bubble(ax,cdata,'original_gene','donor_id',genes,TUD,area=43)
    for tick in ax.get_xticklabels():
        if tick.get_text() in ['PRSS23','CCN2','PLOD2']:tick.set_fontweight('bold')
    cbar(fig,im,[.17,.365,.27,.010],'Mean log1p(expression)')
    dot_legend(fig,.49,.373,area=43)
    def paired_panel(ax,rows,primary,context,xlab,donorlabels):
        for i,row in rows.iterrows():
            y=donorlabels.index(row.donor_id)
            if pd.notna(row.context_donor_median):
                ax.plot([row.context_donor_median,row.NP_donor_median],[y,y],color='#A7A7A7',lw=.8,zorder=1)
                ax.scatter(row.context_donor_median,y,s=22,c=PAL.get(context,'#8A6FA6'),zorder=3)
            ax.scatter(row.NP_donor_median,y,s=24,c=PAL.get(primary,'#168D9C'),zorder=3)
        ax.set_yticks(range(len(donorlabels)),[x.replace('SP','').replace('_',' ') for x in donorlabels]);ax.set_ylim(len(donorlabels)-.4,-.6)
        ax.set_xlabel(xlab);ax.grid(axis='x',color='#EEEEEE',lw=.5);ax.set_axisbelow(True)
    ddata=support[support.cohort.eq('Swahn_GSE230808')].copy()
    save_source(ddata,'Figure2D',h+'GEP2_DONOR_SUPPORT.tsv','paired actual NP-like and other-stromal donor medians; no refit/no new verification')
    ax=panel(fig,[.13,.180,.33,.120],'D','Discovery donor context')
    paired_panel(ax,ddata,'NP_chondrocyte_like','other_disc_stromal','GEP2 UCell',DONS);ax.set_xlim(.3,.8)
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=PAL['NP_chondrocyte_like'],label='NP-like'),
        Line2D([],[],marker='o',ls='',color=PAL['other_disc_stromal'],label='Other stromal')],loc='upper left',bbox_to_anchor=(-.10,-.30),ncol=2,handletextpad=.25,columnspacing=.7)
    save_source(tu_support,'Figure2E',h+'GEP2_DONOR_SUPPORT.tsv',f'{len(paired)} actual NP/immune pairs; absent immune comparators remain unavailable and are not imputed')
    ax=panel(fig,[.62,.180,.33,.120],'E','Independent donor context')
    paired_panel(ax,tu_support,NP,IM,'GEP2 AUCell',TUD);ax.set_xlim(.08,.45)
    ax.set_xticks([.1,.2,.3,.4])
    for missing_donor in missing_comparators:
        ax.text(.447,TUD.index(missing_donor),'NA¹',fontsize=8,ha='right',va='center',color='#777777')
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=PAL[NP],label='NP support'),
        Line2D([],[],marker='o',ls='',color=PAL[IM],label='Immune')],loc='upper left',bbox_to_anchor=(-.10,-.30),ncol=2,handletextpad=.25,columnspacing=.7)
    fdata=brief[['program','mapped_genes_min','mapped_genes_max','expression_interpretation','interpretation_scope','grade_interpretation']].copy()
    save_source(fdata,'Figure2F',h+'ALL_PROGRAMS_BRIEF.tsv','pre-existing hash-bound qualitative review; no new significance score')
    ax=panel(fig,[.09,.018,.87,.052],'F','All-program evidence summary');ax.axis('off')
    labels=['GEP1\nNP context','GEP2\nECM context','GEP3\nMixed lineage','GEP4\nSparse cycling','GEP5\nMixed stress']
    for i,row in enumerate(fdata.itertuples()):
        x=i/5
        ax.add_patch(Rectangle((x,0),.196,1,facecolor='#E6F0EE' if row.expression_interpretation=='CONSISTENT' else '#F3EEE5',edgecolor='white'))
        ax.text(x+.098,.90,labels[i],fontsize=8,ha='center',va='top',linespacing=1.15)
        ax.text(x+.098,.06,'Consistent' if row.expression_interpretation=='CONSISTENT' else 'Mixed',fontsize=8,ha='center',va='bottom')
    # Bottom context detail sits in legend rather than in art.
    save(fig,'Figure2_Frozen_Programs_CrossCohort',out/'04_figures/Figure2')

    # ---- Supplementary Figure S1 ---------------------------------------
    fig=fig_new()
    s1a=cohort[cohort.cohort.eq('Tu_GSE165722')].copy().set_index('donor_id').reindex(TUD).reset_index()
    s1a['NP_supported_cells']=s1a.donor_id.map(tu[NP]);s1a['qc_retention_fraction']=s1a.QC_cells/s1a.submitted_cells
    save_source(s1a,'FigureS1A',h+'COHORT_AND_DONOR_MANIFEST.tsv; donor identity counts','submitted→QC→NP support as nested coverage, not inferential biological effects')
    ax=panel(fig,[.14,.800,.33,.14],'A','Tu donor coverage')
    y=np.arange(nt)
    ax.barh(y,s1a.submitted_cells,color='#E4E4E4',height=.72,label='Submitted')
    ax.barh(y,s1a.QC_cells,color='#6FA6B3',height=.49,label='QC retained')
    ax.barh(y,s1a.NP_supported_cells,color='#168D9C',height=.25,label='NP support')
    ax.set_yticks(y,[d.replace('_',' ') for d in TUD]);ax.invert_yaxis();ax.set_xlabel('Cells');ax.set_xlim(0,7600)
    ax.set_xticks([0,3000,6000],['0','3,000','6,000']);ax.legend(loc='upper left',bbox_to_anchor=(-.1,-.26),ncol=1,labelspacing=.25,handlelength=1)
    s1b=avail.reset_index()[['donor_id','AF_retained_cells','NP_retained_cells']]
    save_source(s1b,'FigureS1B','Figure1a_donor_tissue_counts.tsv','retained source coverage; NP-source counts include broad stromal identities')
    ax=panel(fig,[.63,.800,.32,.14],'B','Swahn paired tissue coverage')
    for i,(col,c,label) in enumerate([('AF_retained_cells','#DA9150','AF stromal'),('NP_retained_cells','#168D9C','NP-source stromal')]):
        ax.scatter(s1b[col],np.arange(nd)+(i-.5)*.16,s=22,c=c,label=label)
    ax.set_yticks(np.arange(nd),[d.replace('SP','') for d in DONS]);ax.invert_yaxis();ax.set_xlabel('Retained cells');ax.set_xlim(0,16000)
    ax.set_xticks([0,7500,15000],['0','7,500','15,000']);ax.legend(loc='upper left',bbox_to_anchor=(-.1,-.26),labelspacing=.25,handlelength=1)
    groups=[('C',['NP_matrix','disc_ECM'],'NP-matrix and stromal markers',.535,.105),
            ('D',['immune','myeloid'],'Immune and myeloid markers',.325,.100),
            ('E',['endothelial'],'Endothelial markers',.12,.105),
            ('F',['contractile'],'Contractile markers',.12,.105)]
    full_max=md.mean_log1p_published_expression.max()
    for letter,panels,title,bottom,height in groups:
        names=list(dict.fromkeys(m[m.marker_panel.isin(panels)].gene))
        df=md[md.gene.isin(names)].copy()
        save_source(df,'FigureS1'+letter,h+'TU_MARKER_DOTPLOT_SOURCE.tsv','equal donor weighting within identity; '+','.join(panels)+' pre-existing marker set')
        if letter in ['E','F']:
            box=[.19 if letter=='E' else .62,bottom,.30,height]
        else:box=[.19,bottom,.75,height]
        ax=panel(fig,box,letter,title);im=bubble(ax,df,'gene','identity',names,IDS,vmax=full_max,area=45)
        if letter=='F':ax.set_yticklabels([])
    cbar(fig,im,[.24,.040,.27,.010],'Mean donor log1p(expression)')
    dot_legend(fig,.56,.048,area=45)
    save(fig,'FigureS1_Donor_QC_and_All_Identity_Markers',out/'06_supplementary/figures')

    # ---- Supplementary Figure S2 ---------------------------------------
    fig=fig_new(240)
    save_source(coverage,'FigureS2A',h+'PROGRAM_GENE_COVERAGE.tsv',f'all{len(P)}programs×{nt}donors; mapped fixed members/25')
    ax=panel(fig,[.17,.895,.73,.060],'A','All-program mapped gene coverage (of 25)')
    im=heat(ax,coverage,'donor_id','program','n_mapped_genes',TUD,P,cmap='Blues',vmin=0,vmax=25,annotate=True,fmt='.0f')
    ax.set_xticklabels([d.replace('_',' ') for d in TUD]);ax.set_yticklabels(P)
    # Raw donor values, equally averaged only among actual available entries.
    ga=g.groupby(['program','original_gene','identity'],as_index=False).agg(
        mean_log1p_published_expression=('mean_log1p_published_expression','mean'),
        fraction_detected=('fraction_detected','mean'),n_contributing_donors=('donor_id','nunique'),
        n_available_donors=('mean_log1p_published_expression','count'))
    ga['aggregation']='equal donor weight among available measurements; unavailable donors not zero-filled'
    ga['expression_interpretation']=ga.program.map(brief.set_index('program').expression_interpretation)
    gmax=ga.mean_log1p_published_expression.max()
    for k,program in enumerate(P):
        letter=chr(ord('B')+k);names=frozen[frozen.program.eq(program)].sort_values('rank').gene.tolist()
        sdf=ga[ga.program.eq(program)].copy()
        save_source(sdf,'FigureS2'+letter,h+'TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv','all25frozen members; all6identities; equal-donor available means; no member substitution')
        bottom=.505 if k<3 else .110
        left=[.18,.50,.82][k%3]
        ax=panel(fig,[left,bottom,.145,.300],letter,program+' · '+brief.set_index('program').loc[program,'expression_interpretation'].capitalize())
        im=bubble(ax,sdf,'identity','original_gene',IDS,names,vmax=gmax,area=25)
        ax.set_xticklabels(['NP','Stroma','Immune','Endoth.','Mixed','Unres.'],rotation=65,ha='right')
    ax=fig.add_axes([.76,.14,.21,.23]);ax.axis('off')
    ax.text(0,1,'Interpretation boundary',fontweight='bold',fontsize=8,va='top')
    boundary='\n\n'.join(p+': '+brief.set_index('program').loc[p,'expression_interpretation'] for p in P[2:])
    grade=' / '.join(brief.grade_interpretation.dropna().drop_duplicates())
    ax.text(0,.88,boundary+'\n\nClinical grade:\n'+grade,fontsize=8,va='top',linespacing=1.25)
    cbar(fig,im,[.22,.040,.27,.009],'Mean donor log1p(expression)')
    dot_legend(fig,.53,.049,area=25)
    save(fig,'FigureS2_All_Five_Frozen_Programs',out/'06_supplementary/figures')

    pd.DataFrame(manifest).drop_duplicates('relative_source').to_csv(qa/'HUMAN_INPUT_MANIFEST.tsv',sep='\t',index=False)
    pd.DataFrame(maps).to_csv(tab/'HUMAN_DATA_TO_FIGURE_MAP.tsv',sep='\t',index=False)
    pd.DataFrame(facts).to_csv(tab/'HUMAN_FACT_REGISTRY.tsv',sep='\t',index=False)
    (qa/'FIGURE_RENDER_QA.json').write_text(json.dumps(qarecs,indent=2),encoding='utf8')
    metadata={'status':'FIGURES_EXPORTED_AWAITING_VISUAL_QA','source_tables':len({r['relative_source'] for r in manifest}),
        'figures':[r['figure'] for r in qarecs],'n_Tu_submitted':int(qc.n_cells.sum()),'n_Tu_QC':int(qc.n_qc_pass.sum()),
        'n_Tu_NP':int(cnt.loc[cnt.identity.eq(NP),'n_cells'].sum()),'n_discovery_stromal':int(cnt.loc[cnt.cohort.eq('Swahn_GSE230808'),'n_cells'].sum()),
        'n_discovery_NP':int(cnt.loc[cnt.identity.eq('NP_chondrocyte_like'),'n_cells'].sum()),'GEP2_Tu_comparisons':{'direction_positive':n_positive,'assessable':len(paired),'donors_total':nt,'missing_comparator':missing_comparators},
        'gene_programs':len(P),'frozen_members_per_program':int(n_members.iloc[0]),'refitting_performed':False,'clinical_grade_test':' / '.join(brief.grade_interpretation.dropna().drop_duplicates()),
        'source_data_panels':len(maps),'automated_outside_text_count':sum(len(r['outside_text']) for r in qarecs),
        'missing_inputs':['No cached per-cell embedding; donor-identity matrix used as explicitly allowed.'],
        'code_publication_performed':False}
    (qa/'HUMAN_COMPLETION.json').write_text(json.dumps(metadata,indent=2),encoding='utf8')
    print(json.dumps(metadata,indent=2))

if __name__=='__main__':main()
