import os,argparse,traceback
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='4')
import anndata as ad,scanpy as sc
from scipy import sparse
from common import *
p=argparse.ArgumentParser();p.add_argument('cohort',choices=COHORTS);args=p.parse_args();c=args.cohort
dest=ROOT/'work'/f'{c}_identity.h5ad'
if dest.exists():print('RESUME_IDENTITY',c);raise SystemExit
a=ad.read_h5ad(OLD/COHORTS[c]/'all_QC_counts_with_lineage.h5ad')
assert sparse.issparse(a.X) and a.obs_names.is_unique and a.var_names.is_unique
assert np.isfinite(a.X.data).all() and np.all(a.X.data>=0) and np.all(a.X.data==np.floor(a.X.data))
assert ((a.obs.n_features>=300)&(a.obs.n_features<=9000)&(a.obs.pct_mito<=25)&(a.obs.n_counts>0)).all()
logs=[];a.obs['doublet_score']=np.nan;a.obs['predicted_doublet']=False
for gsm in a.obs.gsm.unique():
 cache=ROOT/'work'/f'{gsm}_scrublet.tsv'
 if cache.exists(): d=pd.read_csv(cache,sep='\t',index_col=0);record=json.loads(cache.with_suffix('.json').read_text())
 else:
  b=a[a.obs.gsm==gsm].copy()
  record={'cohort':c,'gsm':gsm,'n_QC':len(b),'expected_doublet_rate':.06,'seed':SEED}
  try:
   sc.pp.scrublet(b,expected_doublet_rate=.06,random_state=SEED,verbose=True,n_prin_comps=30)
   d=b.obs[['doublet_score','predicted_doublet']].copy();record.update(status='completed',threshold=float(b.uns['scrublet']['threshold']),predicted=int(d.predicted_doublet.sum()))
  except Exception as e:
   traceback.print_exc();record.update(status='failed',error=repr(e),predicted=0);d=b.obs.assign(doublet_score=np.nan,predicted_doublet=False)[['doublet_score','predicted_doublet']]
  d.to_csv(cache,sep='\t');cache.with_suffix('.json').write_text(json.dumps(record))
 a.obs.loc[d.index,'doublet_score']=d.doublet_score;a.obs.loc[d.index,'predicted_doublet']=d.predicted_doublet.astype(bool)
 logs.append(record);print('SCRUBLET_SAMPLE',record,flush=True)
write(pd.DataFrame(logs),f'tables/{c}_doublet_summary.tsv')
write(a.obs.reset_index(names='cell_id'),f'tables/{c}_all_QC_doublet_audit.tsv.gz')
a=a[~a.obs.predicted_doublet].copy();a.layers['counts']=a.X.copy()
sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat',batch_key='gsm')
b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10);sc.pp.pca(b,n_comps=30,random_state=SEED)
a.obsm['X_pca']=b.obsm['X_pca'];sc.pp.neighbors(a,n_neighbors=15,n_pcs=30,random_state=SEED)
sc.tl.leiden(a,resolution=.5,random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
sc.tl.umap(a,random_state=SEED)
rows=[];top=[];markers=[]
globalmean=np.asarray(a.X.mean(axis=0)).ravel()
for cl in a.obs.leiden.cat.categories:
 mask=a.obs.leiden==cl;x=a.X[mask];mean=np.asarray(x.mean(axis=0)).ravel();frac=np.asarray((x>0).mean(axis=0)).ravel()
 r={'cohort':c,'cluster':cl,'n_cells':int(mask.sum()),'n_donors':a.obs.loc[mask,'donor_id'].nunique(),'largest_donor_fraction':float(a.obs.loc[mask,'donor_id'].value_counts(normalize=True).max()),'median_counts':float(a.obs.loc[mask,'n_counts'].median()),'median_mito':float(a.obs.loc[mask,'pct_mito'].median())}
 for group,genes in MARKERS.items():
  ix=a.var_names.get_indexer([g for g in genes if g in a.var_names]);r[group+'_mean']=float(mean[ix].mean());r[group+'_n_markers']=int((mean[ix]>.35).sum())
  for g in genes:
   j=a.var_names.get_indexer([g])[0];markers.append(dict(cohort=c,cluster=cl,marker_class=group,gene=g,mean_logCP10K=float(mean[j]) if j>=0 else np.nan,fraction=float(frac[j]) if j>=0 else np.nan,available=j>=0))
 # Coherent non-target dominance; otherwise NP-marker positive or other stromal.
 non=[g for g in ['immune','endothelial','mural','erythroid'] if r[g+'_n_markers']>=3 and r[g+'_mean']>max(r['NP_chondrocyte_like_mean'],r['stromal_mean'])*1.25]
 if non: label=max(non,key=lambda g:r[g+'_mean'])
 elif r['NP_chondrocyte_like_n_markers']>=2 and r['immune_mean']<max(.8,1.25*r['NP_chondrocyte_like_mean']): label='NP_chondrocyte_like'
 elif r['stromal_n_markers']>=3 and r['immune_mean']<max(.8,r['stromal_mean']):label='other_disc_stromal'
 else: label='mixed_ambiguous'
 r['identity']=label;rows.append(r)
 for rank,j in enumerate(np.argsort(mean-globalmean)[-50:][::-1],1):top.append(dict(cohort=c,cluster=cl,rank=rank,gene=a.var_names[j],mean_logCP10K=mean[j],fraction=frac[j],difference_from_all=mean[j]-globalmean[j]))
m={r['cluster']:r['identity'] for r in rows};a.obs['identity']=a.obs.leiden.astype(str).map(m)
a.obs['broad_stromal']=a.obs.identity.isin(['NP_chondrocyte_like','other_disc_stromal'])
a.obs['high_confidence_NP']=a.obs.identity=='NP_chondrocyte_like'
write(pd.DataFrame(rows),f'tables/{c}_cluster_identity.tsv');write(pd.DataFrame(markers),f'tables/{c}_cluster_markers.tsv');write(pd.DataFrame(top),f'tables/{c}_cluster_top50.tsv')
write(a.obs.reset_index(names='cell_id'),f'tables/{c}_cell_metadata.tsv.gz')
a.X=a.layers['counts'].copy();del a.layers['counts'];a.write_h5ad(dest,compression='gzip')
print('IDENTITY_COMPLETE',c,a.shape,pd.DataFrame(rows).to_string(index=False),flush=True)
