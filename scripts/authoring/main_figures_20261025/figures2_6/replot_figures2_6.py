#!/usr/bin/env python3
"""Recompose four human and six experimental data panels from frozen V4 tables.

No old artwork is cropped, masked or raster-pasted. No supplementary material,
source-data workbook, model output, or manuscript is written by this command.
"""
import argparse
import hashlib
import json
from pathlib import Path

from _io_guard import independent_output, reject_output_links

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
import pandas as pd
from scipy import stats

NP='NP_matrix_marker_supported'
ST='disc_stromal_marker_supported'
IM='immune_marker_supported'
EN='endothelial_marker_supported'
MI='mixed_lineage_marker_support'
UN='unresolved_low_marker_support'
IDS=[NP,ST,IM,EN,MI,UN]
IDLAB={NP:'NP support',ST:'Disc stromal',IM:'Immune',EN:'Endothelial',MI:'Mixed',UN:'Unresolved'}
PAL={NP:'#168D9C',ST:'#DA9150',IM:'#8A6FA6',EN:'#6E98BD',MI:'#9A9A9A',UN:'#D9D9D9'}
DONS=['SP21.007','SP21.011','SP21.013','SP21.014','SP21.016','SP21.017','SP22.002','SP22.003']
TUD=['Tu_'+str(i) for i in range(1,9)]
MARKERS=['ACAN','COL2A1','SOX9','KRT19','COL1A1','DCN','PTPRC','LST1','CD3D','PECAM1','VWF']
COLORS=['#526578','#1C7F85','#BE6A32']
BATCH_MARKERS=['o','s','^','D']
GROUPS=['G1','G2','G3','G4','G5','G6']
GROUP_LABELS=['NC−','NC+','si1−','si1+','si2−','si2+']
SECONDARY_GENES=['PRSS23','CCN2','COL1A1','ACAN']
SECONDARY_CONTRASTS=['FN_siNC','FN_si1','FN_si2','interaction_si1','interaction_si2']
SOURCE_STATUS='INTERNAL_REVIEW_SOURCE_UNRESOLVED'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_value(obj):
    if isinstance(obj, np.ndarray): return [json_value(v) for v in obj.tolist()]
    if isinstance(obj, (np.integer,np.floating)): return json_value(obj.item())
    if isinstance(obj, dict): return {str(k):json_value(v) for k,v in obj.items()}
    if isinstance(obj, (list,tuple)): return [json_value(v) for v in obj]
    if isinstance(obj, float) and not np.isfinite(obj): return None
    return obj


def signature(obj):
    return hashlib.sha256(json.dumps(json_value(obj),sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def setup():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
        'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
        'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
        'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':2.5,'ytick.major.size':2.5,
        'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})


def label(fig, ax, letter, title, xoffset=.055, title_offset=.015):
    box=ax.get_position()
    fig.text(max(.012,box.x0-xoffset),box.y1+title_offset,letter,fontsize=12,fontweight='bold',va='bottom')
    title_artist=fig.text(box.x0,box.y1+title_offset,title,fontsize=9,fontweight='bold',va='bottom')
    title_artist._qa_role='panel_title'


def add_panel(fig,box,letter,title):
    ax=fig.add_axes(box);label(fig,ax,letter,title);return ax


def human_payloads(source):
    b=source['Figure1B']; c=source['Figure1C']; d=source['Figure1D']; e=source['Figure1E']
    # These are the unchanged original B/C/D/E drawing inputs and transformations.
    tu=c.pivot(index='donor_id',columns='identity',values='n_cells').reindex(index=TUD,columns=IDS).fillna(0)
    return {'Figure2A': {'donors':b.donor_id.tolist(),'counts':b[['AF_retained_cells','NP_like_cells','NP_high_input','NP_low_input']].to_numpy(),
                'included':b.communication_included.tolist()},
        'Figure2B': {'donors':TUD,'identities':IDS,'n_cells':tu.to_numpy(),'fractions':tu.div(tu.sum(axis=1),axis=0).to_numpy()},
        'Figure2C':d.sort_values(['identity','gene']).to_dict(orient='list'),
        'Figure2D':e.pivot(index='identity',columns='donor_id',values='n_cells').reindex(index=['NP_chondrocyte_like','other_disc_stromal'],columns=DONS).to_numpy()}


def draw_human(source):
    b=source['Figure1B']; c=source['Figure1C']; d=source['Figure1D']; e=source['Figure1E']
    fig=plt.figure(figsize=(180/25.4,185/25.4))
    ax=add_panel(fig,[.16,.655,.33,.265],'A','Paired discovery donors')
    for i,row in b.iterrows():
        ax.add_patch(Rectangle((-.5,i-.45),4,.9,facecolor='#E9F3F5' if row.communication_included else '#F1F1F1',edgecolor='white'))
        vals=[str(row.AF_retained_cells),str(row.NP_like_cells),str(row.NP_high_input) if row.communication_included else '—',str(row.NP_low_input) if row.communication_included else '—']
        for j,v in enumerate(vals): ax.text(j,i,v,ha='center',va='center',fontsize=8)
    ax.set_xlim(-.5,3.5);ax.set_ylim(7.5,-.5);ax.set_xticks(range(4),['AF','NP-like','High¹','Low¹'])
    ax.set_yticks(range(8),DONS);ax.tick_params(length=0);ax.spines[:].set_visible(False)
    ax.text(0,-.15,'¹ Selected cells; blue: 6 included pairs.\nGrey: NP tails 12/12 or 4/4; minimum 20 each.',transform=ax.transAxes,fontsize=8,va='top')
    tu=c.pivot(index='donor_id',columns='identity',values='n_cells').reindex(index=TUD,columns=IDS).fillna(0)
    ax=add_panel(fig,[.62,.655,.335,.265],'B','Independent NP-source cohort')
    left=np.zeros(8)
    for identity in IDS:
        frac=tu[identity]/tu.sum(axis=1);ax.barh(range(8),frac,left=left,color=PAL[identity],height=.7,label=IDLAB[identity]);left+=frac.to_numpy()
    ax.set_yticks(range(8),[x.replace('_',' ') for x in TUD]);ax.set_ylim(7.6,-.6);ax.set_xlim(0,1)
    ax.set_xticks([0,.5,1],['0','50','100']);ax.set_xlabel('QC cell composition (%)')
    ax.legend(loc='upper left',bbox_to_anchor=(-.10,-.26),ncol=2,handlelength=.9,columnspacing=.7,labelspacing=.35,borderaxespad=0,frameon=False)
    ax=add_panel(fig,[.20,.305,.72,.175],'C','Canonical identity-marker context')
    xx={x:i for i,x in enumerate(MARKERS)}; yy={y:i for i,y in enumerate(IDS)}
    present=d.mean_log1p_published_expression.notna() & d.fraction_detected.notna()
    norm=colors.Normalize(0,d.mean_log1p_published_expression.max());ok=d[present & d.fraction_detected.gt(0)]
    im=ax.scatter(ok.gene.map(xx),ok.identity.map(yy),s=65*ok.fraction_detected,c=ok.mean_log1p_published_expression,cmap='viridis',norm=norm,linewidths=0,zorder=3)
    zero=d[present & d.fraction_detected.eq(0)]
    ax.scatter(zero.gene.map(xx),zero.identity.map(yy),s=5,facecolors='none',edgecolors='#888888',linewidths=.5,zorder=3)
    missing=d[~present];ax.scatter(missing.gene.map(xx),missing.identity.map(yy),s=14,c='#888888',marker='x',linewidths=.6,zorder=3)
    ax.set_xlim(-.6,len(MARKERS)-.4);ax.set_ylim(len(IDS)-.5,-.5)
    ax.set_xticks(range(len(MARKERS)),MARKERS,rotation=65,ha='right',fontsize=8)
    ax.set_yticks(range(len(IDS)),[IDLAB[y] for y in IDS],fontsize=8);ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)
    for y in range(len(IDS)): ax.axhline(y,color='#F0F0F0',lw=.5,zorder=0)
    cax=fig.add_axes([.23,.193,.28,.014]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label('Mean donor log1p(expression)',fontsize=8,labelpad=2);cb.ax.tick_params(labelsize=8,pad=1,length=2)
    handles=[plt.scatter([],[],s=65*p,c='#555555') for p in [.1,.5,1]]
    fig.legend(handles,['10%','50%','100%'],loc='center left',bbox_to_anchor=(.56,.197),ncol=3,frameon=False,handletextpad=.2,columnspacing=.8,title='Detection (dot area)',title_fontsize=8)
    ax=add_panel(fig,[.17,.040,.76,.060],'D','Discovery identity composition (donor order as A)')
    vals=e.pivot(index='identity',columns='donor_id',values='n_cells').reindex(index=['NP_chondrocyte_like','other_disc_stromal'],columns=DONS)
    ax.imshow(vals.to_numpy(),aspect='auto',cmap='Blues',norm=colors.LogNorm(vmin=1,vmax=float(vals.max().max())))
    ax.set_xticks(range(8),['007','011','013','014','016','017','22.002','22.003'])
    ax.set_yticks(range(2),['NP-like','Other stromal']);ax.tick_params(length=0,pad=2);ax.spines[:].set_visible(False)
    for i in range(2):
        for j in range(8):
            v=int(vals.iloc[i,j]);ax.text(j,i,f'{v:,}',ha='center',va='center',fontsize=8,color='white' if v>1200 else '#222222')
    return fig


def experiment_payloads(source):
    d=source['QPCR_BATCH_DELTACt'];c=source['QPCR_CONTRASTS_AND_CI'];bc=source['QPCR_BATCH_CONTRASTS'];cc=source['CCK8_BATCH_SUMMARY']
    a=[]
    for i,g in enumerate(GROUPS):
        v=np.log2(d[(d.Gene=='ITGB1')&(d.Group_ID==g)].sort_values('Batch_ID').fold_vs_same_batch_G1.values)
        a.append({'group':g,'x':i+np.linspace(-.1,.1,4),'log2_fold':v,'mean':v.mean(),'ci_half_width':stats.t.ppf(.975,3)*stats.sem(v)})
    pb=d[d.Gene=='PLOD2'].copy()
    batches=[{'batch':b,'log2_fold':np.log2(sub.set_index('Group_ID').loc[GROUPS,'fold_vs_same_batch_G1'].to_numpy())} for b,sub in pb.groupby('Batch_ID')]
    fn=[{'batch':b,'fold_ratio':sub.set_index('contrast').loc[['FN_siNC','FN_si1','FN_si2'],'fold_ratio'].to_numpy()} for b,sub in bc[(bc.gene=='PLOD2')&bc.contrast.isin(['FN_siNC','FN_si1','FN_si2'])].groupby('Batch_ID')]
    pc=c[(c.gene=='PLOD2')&c.contrast.str.startswith('interaction')].sort_values('contrast')
    mat=c[c.gene.isin(SECONDARY_GENES)].pivot(index='gene',columns='contrast',values='fold_ratio').loc[SECONDARY_GENES,SECONDARY_CONTRASTS]
    cck=[]
    for g in GROUPS:
        sub=cc[cc.Group_ID==g];summ=sub.groupby('actual_read_hours_after_reseeding').ratio_vs_same_batch_time_G1.agg(['mean','sem'])
        cck.append({'group':g,'summary':summ.reset_index().to_dict(orient='list'),'batch_points':sub[['Batch_ID','actual_read_hours_after_reseeding','ratio_vs_same_batch_time_G1']].to_dict(orient='list')})
    return {'Figure6A':a,'Figure6B':batches,'Figure6C':fn,
        'Figure6D':pc[['contrast','fold_ratio','fold_CI95_low','fold_CI95_high','p_adjusted']].to_dict(orient='list'),
        'Figure6E':{'genes':SECONDARY_GENES,'contrasts':SECONDARY_CONTRASTS,'log2_ratio':np.log2(mat.values)},'Figure6F':cck}


def draw_experiments(source):
    d=source['QPCR_BATCH_DELTACt'];c=source['QPCR_CONTRASTS_AND_CI'];bc=source['QPCR_BATCH_CONTRASTS'];cc=source['CCK8_BATCH_SUMMARY']
    fig=plt.figure(figsize=(180/25.4,225/25.4))
    # Reclaimed schematic height is distributed to six data axes; no empty slot remains.
    gs=fig.add_gridspec(3,2,height_ratios=[1,1,1.15],left=.12,right=.97,top=.94,bottom=.10,hspace=.72,wspace=.63)
    ax=fig.add_subplot(gs[0,0]);label(fig,ax,'A','ITGB1 transcript');x=np.arange(6)
    for i,g in enumerate(GROUPS):
        sub=d[(d.Gene=='ITGB1')&(d.Group_ID==g)].sort_values('Batch_ID');v=np.log2(sub.fold_vs_same_batch_G1.values)
        ax.scatter(i+np.linspace(-.1,.1,4),v,color=COLORS[i//2],marker='o' if i%2==0 else 's',s=15,zorder=3)
        ax.errorbar(i,v.mean(),yerr=stats.t.ppf(.975,3)*stats.sem(v),color='black',fmt='_',capsize=3,lw=.9)
    ax.axhline(0,color='.7',lw=.6);ax.set_xticks(x);ax.set_xticklabels(GROUP_LABELS,rotation=45,ha='right');ax.set_ylabel('log2 expression / same-batch G1')
    ax=fig.add_subplot(gs[0,1]);label(fig,ax,'B','PLOD2 in all six groups')
    for b,sub in d[d.Gene=='PLOD2'].groupby('Batch_ID'):
        v=sub.set_index('Group_ID').loc[GROUPS,'fold_vs_same_batch_G1'];ax.plot(x,np.log2(v),color='.77',lw=.6,zorder=1)
    for i,g in enumerate(GROUPS):
        v=np.log2(d[(d.Gene=='PLOD2')&(d.Group_ID==g)].fold_vs_same_batch_G1);ax.scatter(np.full(4,i),v,color=COLORS[i//2],marker='o' if i%2==0 else 's',s=15,zorder=3)
    ax.axhline(0,color='.7',lw=.6);ax.set_xticks(x);ax.set_xticklabels(GROUP_LABELS,rotation=45,ha='right');ax.set_ylabel('log2 expression / same-batch G1')
    ax=fig.add_subplot(gs[1,0]);label(fig,ax,'C','PLOD2 FN response by batch')
    for i,(b,sub) in enumerate(bc[(bc.gene=='PLOD2')&bc.contrast.isin(['FN_siNC','FN_si1','FN_si2'])].groupby('Batch_ID')):
        v=sub.set_index('contrast').loc[['FN_siNC','FN_si1','FN_si2'],'fold_ratio'];ax.plot(range(3),v,marker=BATCH_MARKERS[i],color='.35',lw=.8,ms=4,label=b)
    ax.axhline(1,color='.7',ls='--',lw=.7);ax.set_xticks(range(3));ax.set_xticklabels(['siNC','si1','si2']);ax.set_ylabel('FN / no added FN');ax.legend(ncol=2,handlelength=1,columnspacing=.5,loc='upper right');ax.set_ylim(.5,2.5)
    ax=fig.add_subplot(gs[1,1]);label(fig,ax,'D','Primary PLOD2 interactions')
    pc=c[(c.gene=='PLOD2')&c.contrast.str.startswith('interaction')].sort_values('contrast')
    for i,(_,r) in enumerate(pc.iterrows()):
        ax.errorbar(r.fold_ratio,1-i,xerr=[[r.fold_ratio-r.fold_CI95_low],[r.fold_CI95_high-r.fold_ratio]],fmt='o',color=COLORS[i+1],capsize=3,ms=5,lw=1.2)
        ax.text(.02,.82-i*.42,'Holm P = '+format(r.p_adjusted,'.3f'),transform=ax.transAxes,fontsize=8)
    ax.axvline(1,color='.5',ls='--',lw=.7);ax.set_yticks([1,0]);ax.set_yticklabels(['si1 / siNC','si2 / siNC']);ax.set_ylim(-.65,1.65);ax.set_xlim(.2,1.4);ax.set_xlabel('FN-response ratio of ratios (95% CI)')
    ax=fig.add_subplot(gs[2,0]);label(fig,ax,'E','Secondary transcript responses')
    sub=c[c.gene.isin(SECONDARY_GENES)];mat=sub.pivot(index='gene',columns='contrast',values='fold_ratio').loc[SECONDARY_GENES,SECONDARY_CONTRASTS];v=np.log2(mat.values)
    im=ax.imshow(v,aspect='auto',cmap='BrBG',norm=colors.TwoSlopeNorm(vmin=-1.5,vcenter=0,vmax=1.5))
    for i in range(4):
        for j in range(5):ax.text(j,i,format(v[i,j],'.2f'),ha='center',va='center',fontsize=8)
    ax.set_yticks(range(4));ax.set_yticklabels(SECONDARY_GENES);ax.set_xticks(range(5));ax.set_xticklabels(['FN NC','FN si1','FN si2','Int si1','Int si2'],rotation=55,ha='right')
    cb=fig.colorbar(im,ax=ax,fraction=.05,pad=.03);cb.set_label('log2 ratio')
    ax=fig.add_subplot(gs[2,1]);label(fig,ax,'F','CCK-8 metabolic signal')
    for i,g in enumerate(GROUPS):
        sub=cc[cc.Group_ID==g];summ=sub.groupby('actual_read_hours_after_reseeding').ratio_vs_same_batch_time_G1.agg(['mean','sem']);ls='-' if i%2 else '--'
        ax.errorbar(summ.index,summ['mean'],yerr=summ['sem'],color=COLORS[i//2],linestyle=ls,lw=1,label=GROUP_LABELS[i],capsize=2)
        for b,q in sub.groupby('Batch_ID'):ax.scatter(q.actual_read_hours_after_reseeding+(int(b[-1])-2.5)*.45,q.ratio_vs_same_batch_time_G1,color=COLORS[i//2],marker='s' if i%2 else 'o',s=8,alpha=.75)
    ax.set_xticks([6,26,50]);ax.set_xlabel('Read time after reseeding (h)');ax.set_ylabel('Matched-blank OD450 / G1');ax.legend(ncol=3,loc='upper center',handlelength=1,columnspacing=.5);ax.set_ylim(.65,1.43)
    fig.text(.5,.020,'INTERNAL REVIEW — experimental source unresolved',ha='center',fontsize=9,color='#8B4B30')
    return fig


def artist_fingerprint(fig):
    """Data-coordinate signature; geometry/relabelled panel letters excluded."""
    records=[]
    for ax in fig.axes:
        rec={'xlim':ax.get_xlim(),'ylim':ax.get_ylim(),'xlabel':ax.get_xlabel(),'ylabel':ax.get_ylabel(),
             'lines':[],'collections':[],'images':[],'bars':[]}
        for line in ax.lines:
            rec['lines'].append({'x':np.asarray(line.get_xdata()).tolist(),'y':np.asarray(line.get_ydata()).tolist(),
                'color':str(line.get_color()),'marker':str(line.get_marker()),'linestyle':str(line.get_linestyle())})
        for collection in ax.collections:
            c={'offsets':np.ma.filled(collection.get_offsets(),np.nan),'array':None if collection.get_array() is None else np.ma.filled(collection.get_array(),np.nan)}
            if hasattr(collection,'get_sizes'):c['sizes']=collection.get_sizes()
            if hasattr(collection,'get_segments'):c['segments']=collection.get_segments()
            rec['collections'].append(c)
        for image in ax.images:rec['images'].append({'values':np.ma.filled(image.get_array(),np.nan),'clim':image.get_clim(),'cmap':image.get_cmap().name})
        for patch in ax.patches:
            if isinstance(patch,Rectangle):rec['bars'].append([patch.get_x(),patch.get_y(),patch.get_width(),patch.get_height()])
        records.append(rec)
    return json_value(records)


def audit_and_export(fig,stem,out):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();W,H=fig.canvas.get_width_height()
    visible=[];outside=[];small=[];overlaps=[]
    for a in fig.findobj(Text):
        if not a.get_visible() or not a.get_text().strip():continue
        bb=a.get_window_extent(renderer);visible.append((a,bb))
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>W+.5 or bb.y1>H+.5:outside.append({'text':a.get_text(),'bbox_px':list(map(float,bb.extents))})
        if a.get_fontsize()<8-1e-9:small.append({'text':a.get_text(),'font_pt':a.get_fontsize()})
    for a,bb in visible:
        if getattr(a,'_qa_role','')!='panel_title':continue
        for b,bc in visible:
            if a is not b and bb.overlaps(bc):overlaps.append({'title':a.get_text(),'overlap':b.get_text()})
    for ext,kw in [('svg',{}),('pdf',{}),('png',{'dpi':300}),('tiff',{'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}})]:
        fig.savefig(out/(stem+'.'+ext),**kw)
    import fitz
    with fitz.open(out/(stem+'.pdf')) as doc:
        page=doc[0]
        for dpi in [100,200]:page.get_pixmap(matrix=fitz.Matrix(dpi/72,dpi/72),alpha=False).save(out/(stem+'_PDF_'+str(dpi)+'dpi.png'))
        pdf_text=page.get_text();spans=[s for b in page.get_text('dict')['blocks'] if 'lines'in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        pdf_min=min(s['size'] for s in spans)
        dimensions=[page.rect.width/72*25.4,page.rect.height/72*25.4]
    qa={'figure':stem,'status':'PASS_AUTOMATED' if not(outside or small or overlaps) and pdf_min>=7.99 else 'REVISE',
        'width_mm':dimensions[0],'height_mm':dimensions[1],'declared_font':'DejaVu Sans','min_pdf_text_pt':pdf_min,
        'minimum_artist_text_pt':min(a.get_fontsize() for a,_ in visible),'outside_canvas_text':outside,'small_text':small,'title_text_overlaps':overlaps,
        'editable_pdf_text_characters':len(pdf_text),'manual_visual_status':'PENDING_ACTUAL_REVIEW','source_status':SOURCE_STATUS if stem=='Figure6' else 'COMPUTATIONAL_VERIFIED',
        'no_raster_artwork_reuse':True,'no_analysis_refit':True}
    (out/(stem+'_AUTOMATED_QA.json')).write_text(json.dumps(qa,indent=2),encoding='utf8')
    return qa


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);ap.add_argument('--output-root',type=Path,required=True)
    ap.add_argument('--figures',choices=['both','2','6'],default='both');a=ap.parse_args()
    src,out=independent_output(a.source_root,a.output_root)
    reject_output_links(out)
    out.mkdir(parents=True,exist_ok=True);setup();manifest=[];qa=[];signatures=[]
    def read(path):
        frame=pd.read_csv(path,sep='\t');manifest.append({'source':path.relative_to(src).as_posix(),'sha256_before':sha(path),'bytes':path.stat().st_size,'rows':len(frame)});return frame
    if a.figures in ['both','2']:
        source={key:read(src/'03_tables/human'/(key+'_source.tsv')) for key in ['Figure1B','Figure1C','Figure1D','Figure1E']}
        payload=human_payloads(source);fig=draw_human(source)
        signatures.extend({'old_panel':'Figure1'+old,'new_panel':new,'data_signature':signature(value),'comparison':'unchanged original panel-source fields and original drawing transformations'} for old,(new,value) in zip('BCDE',payload.items()))
        qa.append(audit_and_export(fig,'Figure2',out));plt.close(fig)
    if a.figures in ['both','6']:
        s=src/'03_tables/experimental/internal/source_unresolved'
        source={key:read(s/(key+'.tsv')) for key in ['QPCR_BATCH_DELTACt','QPCR_CONTRASTS_AND_CI','QPCR_BATCH_CONTRASTS','CCK8_BATCH_SUMMARY']}
        if 'SOURCE_STATUS' in source['QPCR_CONTRASTS_AND_CI']:
            assert set(source['QPCR_CONTRASTS_AND_CI'].SOURCE_STATUS)=={SOURCE_STATUS}
        payload=experiment_payloads(source);fig=draw_experiments(source)
        signatures.extend({'old_panel':'Figure5'+old,'new_panel':new,'data_signature':signature(value),'comparison':'unchanged original panel-source fields, same log transforms, errorbar formulas, contrast selections and axis limits'} for old,(new,value) in zip('BCDEFG',payload.items()))
        qa.append(audit_and_export(fig,'Figure6',out));plt.close(fig)
    for record in manifest:
        record['sha256_after']=sha(src/record['source']);record['input_unchanged']=record['sha256_before']==record['sha256_after']
    assert all(m['input_unchanged'] for m in manifest)
    pd.DataFrame(manifest).to_csv(out/'REPLOT_INPUT_HASH_AUDIT.tsv',sep='\t',index=False)
    pd.DataFrame(signatures).to_csv(out/'REPLOT_PANEL_NUMERIC_SIGNATURES.tsv',sep='\t',index=False)
    (out/'REPLOT_COMPLETION.json').write_text(json.dumps({'figures':qa,'all_inputs_byte_unchanged':True,'new_analysis':False,'supplementary_exports':False,'raw_source_exported':False,'legacy_schematic_panels_removed':['Figure1A','Figure5A']},indent=2),encoding='utf8')
    print(json.dumps({'outputs':str(out),'figures':[{k:v for k,v in q.items() if k in ['figure','status','width_mm','height_mm','min_pdf_text_pt']} for q in qa]},indent=2))


if __name__=='__main__': main()
