"""Internal-only experimental figures; this script cannot release provenance."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy import stats

COLORS=['#526578','#1C7F85','#BE6A32'];MARKERS=['o','s','^','D']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7})
def label(ax,letter,title):
 pos=ax.get_position();y=1.07 if pos.y1>.9 and pos.width<.5 else 1.14
 ax.set_title(title,loc='left',pad=9);ax.text(-.07 if pos.width>.5 else -.16,y,letter,transform=ax.transAxes,fontweight='bold',fontsize=12,va='bottom')
def save(fig,dest,source):
 dest.parent.mkdir(parents=True,exist_ok=True)
 fig.text(.5,.016,'INTERNAL REVIEW — experimental source unresolved',ha='center',fontsize=9,color='#8B4B30')
 fig.canvas.draw();rend=fig.canvas.get_renderer();fb=fig.bbox;overflow=[]
 for t in fig.findobj(matplotlib.text.Text):
  if not t.get_visible() or not t.get_text():continue
  b=t.get_window_extent(rend)
  if b.x0<0 or b.y0<0 or b.x1>fb.x1 or b.y1>fb.y1:overflow.append(t.get_text())
 for ext,kw in [('pdf',{}),('svg',{}),('png',{'dpi':300}),('tiff',{'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}})]:fig.savefig(str(dest)+'.'+ext,**kw)
 source.write_text(json.dumps({'width_mm':180,'height_mm':fig.get_size_inches()[1]*25.4,'minimum_text_pt':8,'source_status':'INTERNAL_REVIEW_SOURCE_UNRESOLVED','text_outside_canvas':overflow,'bbox_status':'PASS' if not overflow else 'REVISE'},indent=2),encoding='utf-8')
 plt.close(fig)
def main():
 p=argparse.ArgumentParser();p.add_argument('--project',required=True);a=p.parse_args();root=Path(a.project);s=root/'03_tables/experimental/internal/source_unresolved'
 d=pd.read_csv(s/'QPCR_BATCH_DELTACt.tsv',sep='\t');c=pd.read_csv(s/'QPCR_CONTRASTS_AND_CI.tsv',sep='\t');bc=pd.read_csv(s/'QPCR_BATCH_CONTRASTS.tsv',sep='\t');cc=pd.read_csv(s/'CCK8_BATCH_SUMMARY.tsv',sep='\t');tech=pd.read_csv(s/'QPCR_TECHNICAL_MEANS.tsv',sep='\t');bl=pd.read_csv(s/'CCK8_MATCHED_BLANKS.tsv',sep='\t')
 fdir=root/'04_figures/Figure5_internal/source_unresolved';fdir.mkdir(parents=True,exist_ok=True);(fdir/'source_data').mkdir(exist_ok=True)
 for n in ['EXPERIMENT_DESIGN.tsv','QPCR_BATCH_DELTACt.tsv','QPCR_CONTRASTS_AND_CI.tsv','QPCR_BATCH_CONTRASTS.tsv','CCK8_BATCH_SUMMARY.tsv','CCK8_BLANK_CORRECTED.tsv']:shutil.copy2(s/n,fdir/'source_data'/n)
 fig=plt.figure(figsize=(180/25.4,238/25.4));gs=fig.add_gridspec(4,2,height_ratios=[.65,1.15,1.15,1.3],left=.11,right=.97,top=.955,bottom=.085,hspace=.65,wspace=.58)
 ax=fig.add_subplot(gs[0,:]);label(ax,'A','Six groups in four labelled culture/transfection batches');ax.axis('off')
 for x,title,text in [(.14,'Transfection','siNC / siITGB1_1 / siITGB1_2\n10 nM; two source wells pooled'),(.49,'Reseeding at 48 h','Equal cell input, ± plasma FN coating\n10 µg/mL; nominal 2 µg/cm²'),(.84,'Readouts','RNA: 72 h after transfection\nCCK: 6, 26, 50 h after reseeding')]:
  ax.text(x,.64,title,ha='center',fontweight='bold',fontsize=8);ax.text(x,.22,text,ha='center',va='center',fontsize=8)
 ax.annotate('',(.35,.60),(.25,.60),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.8});ax.annotate('',(.69,.60),(.60,.60),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.8})
 ax=fig.add_subplot(gs[1,0]);label(ax,'B','ITGB1 transcript');x=np.arange(6)
 for i,g in enumerate(['G1','G2','G3','G4','G5','G6']):
  sub=d[(d.Gene=='ITGB1')&(d.Group_ID==g)].sort_values('Batch_ID');v=np.log2(sub.fold_vs_same_batch_G1.values);ax.scatter(i+np.linspace(-.1,.1,4),v,color=COLORS[i//2],marker='o' if i%2==0 else 's',s=15,zorder=3);ax.errorbar(i,v.mean(),yerr=stats.t.ppf(.975,3)*stats.sem(v),color='black',fmt='_',capsize=3,lw=.9)
 ax.axhline(0,color='.7',lw=.6);ax.set_xticks(x);ax.set_xticklabels(['NC−','NC+','si1−','si1+','si2−','si2+'],rotation=45,ha='right');ax.set_ylabel('log2 expression / same-batch G1')
 ax=fig.add_subplot(gs[1,1]);label(ax,'C','PLOD2 in all six groups')
 for b,sub in d[d.Gene=='PLOD2'].groupby('Batch_ID'):
  v=sub.set_index('Group_ID').loc[['G1','G2','G3','G4','G5','G6'],'fold_vs_same_batch_G1'];ax.plot(x,np.log2(v),color='.77',lw=.6,zorder=1)
 for i,g in enumerate(['G1','G2','G3','G4','G5','G6']):
  v=np.log2(d[(d.Gene=='PLOD2')&(d.Group_ID==g)].fold_vs_same_batch_G1);ax.scatter(np.full(4,i),v,color=COLORS[i//2],marker='o' if i%2==0 else 's',s=15,zorder=3)
 ax.axhline(0,color='.7',lw=.6);ax.set_xticks(x);ax.set_xticklabels(['NC−','NC+','si1−','si1+','si2−','si2+'],rotation=45,ha='right');ax.set_ylabel('log2 expression / same-batch G1')
 ax=fig.add_subplot(gs[2,0]);label(ax,'D','PLOD2 FN response by batch')
 for i,(b,sub) in enumerate(bc[(bc.gene=='PLOD2')&bc.contrast.isin(['FN_siNC','FN_si1','FN_si2'])].groupby('Batch_ID')):
  v=sub.set_index('contrast').loc[['FN_siNC','FN_si1','FN_si2'],'fold_ratio'];ax.plot(range(3),v,marker=MARKERS[i],color='.35',lw=.8,ms=4,label=b)
 ax.axhline(1,color='.7',ls='--',lw=.7);ax.set_xticks(range(3));ax.set_xticklabels(['siNC','si1','si2']);ax.set_ylabel('FN / no added FN');ax.legend(ncol=2,handlelength=1,columnspacing=.5,loc='upper right');ax.set_ylim(.5,2.5)
 ax=fig.add_subplot(gs[2,1]);label(ax,'E','Primary PLOD2 interactions')
 pc=c[(c.gene=='PLOD2')&c.contrast.str.startswith('interaction')].sort_values('contrast')
 for i,(_,r) in enumerate(pc.iterrows()):
  ax.errorbar(r.fold_ratio,1-i,xerr=[[r.fold_ratio-r.fold_CI95_low],[r.fold_CI95_high-r.fold_ratio]],fmt='o',color=COLORS[i+1],capsize=3,ms=5,lw=1.2);ax.text(.02,.82-i*.42,'Holm P = '+format(r.p_adjusted,'.3f'),transform=ax.transAxes,fontsize=8)
 ax.axvline(1,color='.5',ls='--',lw=.7);ax.set_yticks([1,0]);ax.set_yticklabels(['si1 / siNC','si2 / siNC']);ax.set_ylim(-.65,1.65);ax.set_xlim(.2,1.4);ax.set_xlabel('FN-response ratio of ratios (95% CI)')
 ax=fig.add_subplot(gs[3,0]);label(ax,'F','Secondary transcript responses')
 genes=['PRSS23','CCN2','COL1A1','ACAN'];cts=['FN_siNC','FN_si1','FN_si2','interaction_si1','interaction_si2'];sub=c[c.gene.isin(genes)];mat=sub.pivot(index='gene',columns='contrast',values='fold_ratio').loc[genes,cts];v=np.log2(mat.values);im=ax.imshow(v,aspect='auto',cmap='BrBG',norm=colors.TwoSlopeNorm(vmin=-1.5,vcenter=0,vmax=1.5))
 for i in range(4):
  for j in range(5):ax.text(j,i,format(v[i,j],'.2f'),ha='center',va='center',fontsize=8)
 ax.set_yticks(range(4));ax.set_yticklabels(genes);ax.set_xticks(range(5));ax.set_xticklabels(['FN NC','FN si1','FN si2','Int si1','Int si2'],rotation=55,ha='right');cb=fig.colorbar(im,ax=ax,fraction=.05,pad=.03);cb.set_label('log2 ratio')
 ax=fig.add_subplot(gs[3,1]);label(ax,'G','CCK-8 metabolic signal')
 for i,g in enumerate(['G1','G2','G3','G4','G5','G6']):
  sub=cc[cc.Group_ID==g];summ=sub.groupby('actual_read_hours_after_reseeding').ratio_vs_same_batch_time_G1.agg(['mean','sem']);ls='-' if i%2 else '--';ax.errorbar(summ.index,summ['mean'],yerr=summ['sem'],color=COLORS[i//2],linestyle=ls,lw=1,label=['NC−','NC+','si1−','si1+','si2−','si2+'][i],capsize=2)
  for b,q in sub.groupby('Batch_ID'):ax.scatter(q.actual_read_hours_after_reseeding+(int(b[-1])-2.5)*.45,q.ratio_vs_same_batch_time_G1,color=COLORS[i//2],marker='s' if i%2 else 'o',s=8,alpha=.75)
 ax.set_xticks([6,26,50]);ax.set_xlabel('Read time after reseeding (h)');ax.set_ylabel('Matched-blank OD450 / G1');ax.legend(ncol=3,loc='upper center',handlelength=1,columnspacing=.5);ax.set_ylim(.65,1.43)
 save(fig,fdir/'Figure5_Cell_Perturbation_INTERNAL',fdir/'BBOX_QA.json')
 # Supplementary quality/data displays remain under the same source gate.
 sout=root/'06_supplementary/internal/source_unresolved';sout.mkdir(parents=True,exist_ok=True)
 fig,axs=plt.subplots(2,2,figsize=(180/25.4,185/25.4));fig.subplots_adjust(left=.12,right=.97,top=.94,bottom=.12,hspace=.45,wspace=.45)
 ax=axs[0,0];label(ax,'A','Technical Ct dispersion');gsym=sorted(tech.Gene.unique());ax.boxplot([tech[tech.Gene==g].technical_sd for g in gsym],positions=range(len(gsym)),showfliers=True,patch_artist=False);ax.set_xticks(range(len(gsym)));ax.set_xticklabels(gsym,rotation=60,ha='right');ax.set_ylabel('Within-sample technical Ct SD')
 ax=axs[0,1];label(ax,'B','Reference gene Ct values')
 for j,g in enumerate(['RPLP0','TBP']):
  for i,gr in enumerate(['G1','G2','G3','G4','G5','G6']):
   v=tech[(tech.Gene==g)&(tech.Group_ID==gr)].mean_Ct;ax.scatter(np.full(4,i)+j*.15-.075,v,s=12,color=COLORS[j],label=g if i==0 else None)
 ax.set_xticks(range(6));ax.set_xticklabels(['G1','G2','G3','G4','G5','G6']);ax.set_ylabel('Mean technical Ct');ax.legend()
 ax=axs[1,0];label(ax,'C','Secondary gene expression')
 for j,g in enumerate(genes):
  sub=d[d.Gene==g];xx=np.arange(6)+(j-1.5)*.15;means=sub.groupby('Group_ID').fold_vs_same_batch_G1.apply(lambda v:np.mean(np.log2(v)));ax.plot(xx,means.values,marker=MARKERS[j],lw=.8,label=g,ms=3)
 ax.set_xticks(range(6));ax.set_xticklabels(['G1','G2','G3','G4','G5','G6']);ax.set_ylabel('Mean batch log2 expression / G1');ax.set_ylim(-.4,.9);ax.legend(fontsize=8,ncol=2,loc='upper right')
 ax=axs[1,1];label(ax,'D','Recorded control status');ctrl=pd.read_csv(s/'QPCR_NTC_NO_RT_RECORDED.tsv',sep='\t');tab=ctrl.groupby(['Gene','Sample_type']).size().unstack(fill_value=0);tab.plot.bar(ax=ax,color=['#9AA3AE','#B4C8CD'],width=.8);ax.set_ylabel('Recorded wells');ax.set_xlabel('');ax.tick_params(axis='x',rotation=60);ax.legend(title=None,fontsize=8)
 save(fig,sout/'Supplementary_Figure_S5',sout/'S5_BBOX_QA.json')
 fig,axs=plt.subplots(2,2,figsize=(180/25.4,185/25.4));fig.subplots_adjust(left=.12,right=.97,top=.94,bottom=.12,hspace=.46,wspace=.45)
 ax=axs[0,0];label(ax,'A','Matched backgrounds')
 for fn,sub in bl.groupby('FN_condition'):
  for b,q in sub.groupby('Batch_ID'):ax.plot(q.FN_hours+2,q.matched_blank_mean_OD450,lw=.6,marker='s' if fn.startswith('FN_') else 'o',ms=3,color=COLORS[1] if fn.startswith('FN_') else COLORS[0],label=fn if b=='B01' else None)
 ax.set_ylabel('Mean blank OD450');ax.set_xlabel('Actual read time after reseeding (h)');ax.set_xticks([6,26,50]);ax.legend(labels=['FN','Carrier'],fontsize=8)
 ax=axs[0,1];label(ax,'B','Blank-corrected batch signals')
 for i,g in enumerate(['G1','G2','G3','G4','G5','G6']):
  sub=cc[cc.Group_ID==g];ss=sub.groupby('actual_read_hours_after_reseeding').mean_corrected_OD450.mean();ax.plot(ss.index,ss,color=COLORS[i//2],ls='-' if i%2 else '--',label=g)
 ax.set_ylabel('Mean corrected OD450');ax.set_xlabel('Actual read time after reseeding (h)');ax.set_xticks([6,26,50]);ax.legend(ncol=3)
 ax=axs[1,0];label(ax,'C','All batch ratios at each time')
 for i,g in enumerate(['G1','G2','G3','G4','G5','G6']):
  sub=cc[cc.Group_ID==g]
  for b,q in sub.groupby('Batch_ID'):ax.plot(q.actual_read_hours_after_reseeding,q.ratio_vs_same_batch_time_G1,color=COLORS[i//2],ls='-' if i%2 else '--',lw=.65,marker=MARKERS[int(b[-1])-1],ms=3)
 ax.set_ylabel('Corrected OD450 / same-time G1');ax.set_xlabel('Actual read time after reseeding (h)');ax.set_xticks([6,26,50])
 ax=axs[1,1];label(ax,'D','Auxiliary change from own 4 h')
 for i,g in enumerate(['G1','G2','G3','G4','G5','G6']):
  sub=cc[cc.Group_ID==g];ss=sub.groupby('actual_read_hours_after_reseeding').ratio_vs_own_4h.mean();ax.plot(ss.index,ss,color=COLORS[i//2],ls='-' if i%2 else '--',marker='o',ms=3,label=g)
 ax.set_ylabel('Corrected OD450 / own 4 h');ax.set_xlabel('Actual read time after reseeding (h)');ax.set_xticks([6,26,50]);ax.legend(ncol=3)
 save(fig,sout/'Supplementary_Figure_S6',sout/'S6_BBOX_QA.json')
 print('WROTE Figure5,S5,S6: internal only')
if __name__=='__main__':main()
