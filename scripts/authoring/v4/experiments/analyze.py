"""Batch-level numerical analysis of supplied Ct/OD, source status independent."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

STATUS='INTERNAL_REVIEW_SOURCE_UNRESOLVED'
def adjust(p,method):
 p=np.asarray(p,float);order=np.argsort(p);out=np.empty(len(p));q=p[order];n=len(p)
 if method=='Holm':v=np.maximum.accumulate(q*(n-np.arange(n)))
 else:v=np.minimum.accumulate((q*n/np.arange(1,n+1))[::-1])[::-1]
 out[order]=np.minimum(1,v);return out
def summarize(vals):
 v=np.asarray(vals,float);v=v[np.isfinite(v)];n=len(v);mean=float(v.mean());se=float(stats.sem(v));ci=stats.t.interval(.95,n-1,loc=mean,scale=se);p=float(stats.ttest_1samp(v,0).pvalue)
 return dict(n_batches=n,mean_deltaCt_contrast=mean,deltaCt_CI95_low=ci[0],deltaCt_CI95_high=ci[1],fold_ratio=2**(-mean),fold_CI95_low=2**(-ci[1]),fold_CI95_high=2**(-ci[0]),p_two_sided=p,df=n-1)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--project',required=True);a=ap.parse_args();root=Path(a.project);out=root/'03_tables/experimental/internal/source_unresolved';raw=json.loads((out/'RAW_EXTRACT.json').read_text(encoding='utf-8'))
 # Fixed before examining inferential outputs; all requested outcomes retained.
 plan={'primary':'PLOD2 two FN×siRNA difference-in-differences; two-sided one-sample t,df3;Holm2','secondary':'PRSS23,CCN2,COL1A1,ACAN each FN response(siNC,si1,si2) and2interactions;one fixed BH20 family','supporting':'PLOD2 FN responses3 and ITGB1 knockdown4;unadjusted explicitly supporting','CCK8':'descriptive matched-blank corrected batch means;no equivalence/toxicity or time-search testing','outlier_removal':'none','imputation':'none','source_status':STATUS}
 (out/'FIXED_ANALYSIS_PLAN.json').write_text(json.dumps(plan,indent=2))
 def save(d,name):
  d=d.copy();d['ANALYTIC_STATUS']='NUMERIC_RECALCULATION_DONE';d['SOURCE_STATUS']=STATUS;d.to_csv(out/name,sep='\t',index=False)
 q=pd.concat([pd.DataFrame(raw['siRNA_qPCR_Raw']),pd.DataFrame(raw['Program_qPCR_Raw'])],ignore_index=True)
 s=q[q.Sample_type=='Sample'].copy();s['Ct_numeric']=pd.to_numeric(s.Raw_Ct,errors='coerce');assert s.Ct_numeric.notna().all();assert len(s)==4*6*8*3
 means=s.groupby(['Batch_ID','Group_ID','siRNA','FN_condition','Gene'],as_index=False).agg(mean_Ct=('Ct_numeric','mean'),technical_n=('Ct_numeric','size'),technical_sd=('Ct_numeric','std'),technical_min=('Ct_numeric','min'),technical_max=('Ct_numeric','max'))
 assert means.technical_n.eq(3).all();save(means,'QPCR_TECHNICAL_MEANS.tsv')
 keys=['Batch_ID','Group_ID'];refs=means[means.Gene.isin(['RPLP0','TBP'])].pivot(index=keys,columns='Gene',values='mean_Ct').reset_index();refs['mean_reference_Ct']=(refs.RPLP0+refs.TBP)/2
 d=means[~means.Gene.isin(['RPLP0','TBP'])].merge(refs,on=keys,validate='many_to_one');d['deltaCt']=d.mean_Ct-d.mean_reference_Ct
 g1=d[d.Group_ID=='G1'][['Batch_ID','Gene','deltaCt']].rename(columns={'deltaCt':'same_batch_G1_deltaCt'});d=d.merge(g1,on=['Batch_ID','Gene'],validate='many_to_one');d['deltaDeltaCt']=d.deltaCt-d.same_batch_G1_deltaCt;d['fold_vs_same_batch_G1']=2**(-d.deltaDeltaCt);save(d,'QPCR_BATCH_DELTACt.tsv')
 contrasts=[];batchcontrasts=[]
 for gene,sub in d.groupby('Gene'):
  mat=sub.pivot(index='Batch_ID',columns='Group_ID',values='deltaCt');fn={'FN_siNC':mat.G2-mat.G1,'FN_si1':mat.G4-mat.G3,'FN_si2':mat.G6-mat.G5};fn.update({'interaction_si1':fn['FN_si1']-fn['FN_siNC'],'interaction_si2':fn['FN_si2']-fn['FN_siNC']})
  if gene=='ITGB1':fn.update({'KD_si1_noFN':mat.G3-mat.G1,'KD_si1_FN':mat.G4-mat.G2,'KD_si2_noFN':mat.G5-mat.G1,'KD_si2_FN':mat.G6-mat.G2})
  for name,vals in fn.items():
   kind='primary_Holm2' if gene=='PLOD2' and name.startswith('interaction') else ('secondary_BH20' if gene in ['PRSS23','CCN2','COL1A1','ACAN'] else 'supporting_unadjusted')
   r=dict(gene=gene,contrast=name,family=kind,**summarize(vals));contrasts.append(r)
   for batch,value in vals.items():batchcontrasts.append(dict(gene=gene,contrast=name,Batch_ID=batch,deltaCt_contrast=value,fold_ratio=2**(-value),family=kind))
 c=pd.DataFrame(contrasts);c['p_adjusted']=np.nan;c['adjustment']='none_supporting'
 for family,method in [('primary_Holm2','Holm'),('secondary_BH20','BH')]:
  mask=c.family==family;c.loc[mask,'p_adjusted']=adjust(c.loc[mask,'p_two_sided'],method);c.loc[mask,'adjustment']=method+'_'+str(mask.sum())
 assert (c.family=='primary_Holm2').sum()==2;assert (c.family=='secondary_BH20').sum()==20
 save(c,'QPCR_CONTRASTS_AND_CI.tsv');save(pd.DataFrame(batchcontrasts),'QPCR_BATCH_CONTRASTS.tsv')
 cells=pd.DataFrame(raw['CCK8_Raw']);cells['Raw_OD450']=pd.to_numeric(cells.Raw_OD450,errors='raise');blank=cells[cells.Well_type=='Blank'].groupby(['Batch_ID','FN_hours','FN_condition','Plate_ID','Blank_key'],as_index=False).agg(matched_blank_mean_OD450=('Raw_OD450','mean'),blank_n=('Raw_OD450','size'))
 assert len(blank)==24 and blank.blank_n.eq(3).all()
 cc=cells[cells.Well_type=='Cells'].merge(blank,on=['Batch_ID','FN_hours','FN_condition','Plate_ID','Blank_key'],validate='many_to_one');cc['corrected_OD450']=cc.Raw_OD450-cc.matched_blank_mean_OD450;cc['actual_read_hours_after_reseeding']=cc.FN_hours+2;save(cc,'CCK8_BLANK_CORRECTED.tsv');save(blank,'CCK8_MATCHED_BLANKS.tsv')
 cb=cc.groupby(['Batch_ID','Group_ID','siRNA','FN_condition','FN_hours','actual_read_hours_after_reseeding'],as_index=False).agg(mean_corrected_OD450=('corrected_OD450','mean'),tech_n=('corrected_OD450','size'),mean_raw_OD450=('Raw_OD450','mean'),matched_blank_mean_OD450=('matched_blank_mean_OD450','first'))
 ref=cb[cb.Group_ID=='G1'][['Batch_ID','FN_hours','mean_corrected_OD450']].rename(columns={'mean_corrected_OD450':'same_batch_time_G1_corrected_OD450'});cb=cb.merge(ref,on=['Batch_ID','FN_hours'],validate='many_to_one');cb['ratio_vs_same_batch_time_G1']=cb.mean_corrected_OD450/cb.same_batch_time_G1_corrected_OD450
 own4=cb[cb.FN_hours==4][['Batch_ID','Group_ID','mean_corrected_OD450']].rename(columns={'mean_corrected_OD450':'own_group_4h_corrected_OD450'});cb=cb.merge(own4,on=['Batch_ID','Group_ID'],validate='many_to_one');cb['ratio_vs_own_4h']=cb.mean_corrected_OD450/cb.own_group_4h_corrected_OD450
 assert cb.tech_n.eq(3).all();assert len(cb)==72;save(cb,'CCK8_BATCH_SUMMARY.tsv')
 des=d[['Batch_ID','Group_ID','siRNA','FN_condition']].drop_duplicates();des['RNA_hours_after_transfection']=72;des['RNA_hours_after_reseeding']=24;des['siRNA_nM']=10;save(des,'EXPERIMENT_DESIGN.tsv')
 qc=[('Original workbook supplied','YES','Private raw Ct/OD flat tables; observation counts are in the generated receipt'),('Independent instrument run exports','NO','Provenance/source conflict unresolved'),('NTC and No-RT main table entries','YES','Recorded status only;no curves'),('Standard dilution curves and efficiency raw records','NO','Protocol states thresholds;not supplied results'),('Amplification/melt traces','NO','No fabricated traces'),('ITGB1 second amplicon raw comparison','NO','Mainc2 overlaps si2target;independentc3QC missing'),('RPLP0/TBP observed technical Ct','YES','Descriptive values;stability validation not supplied'),('Cell authentication/mycoplasma/passage record','NO','Not invented'),('Donor or patient biological replication','NOT_APPLICABLE','One immortalized commercial cell source;4batchIDs'),('Nonpositive corrected CCK OD450',str(int((cc.corrected_OD450<=0).sum())),'Retained;no clipping')]
 save(pd.DataFrame(qc,columns=['item','availability','interpretation']),'EXPERIMENT_QC_AVAILABILITY.tsv')
 save(q[q.Sample_type!='Sample'],'QPCR_NTC_NO_RT_RECORDED.tsv')
 r=c[(c.gene=='PLOD2')&c.contrast.str.startswith('interaction')];report={'ANALYTIC_STATUS':'COMPLETE_BATCH_LEVEL_NUMERIC_RECALCULATION','SOURCE_STATUS':STATUS,'raw_qPCR_rows':len(q),'qPCR_sample_wells':len(s),'raw_CCK8_rows':len(cells),'independent_batch_ids':sorted(d.Batch_ID.unique()),'primary_interactions':r.to_dict('records'),'nonpositive_corrected_OD':int((cc.corrected_OD450<=0).sum()),'full_clean_experimental_source_verified':False}
 (out/'EXPERIMENT_ANALYSIS_STATUS.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
