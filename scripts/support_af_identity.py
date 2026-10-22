#!/usr/bin/env python3
"""One fixed AF QC, per-library Scrublet and marker-based annotation pass."""
import os
os.environ.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='4')
from pathlib import Path
import datetime as dt, json, traceback
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import argparse, hashlib
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--raw-dir',type=Path,required=True)
p.add_argument('--out-dir',type=Path,required=True)
p.add_argument('--cache-dir',type=Path,required=True)
p.add_argument('--sample-manifest',type=Path,required=True)
p.add_argument('--seed',type=int,default=20261017)
args=p.parse_args()
RAW=args.raw_dir.resolve();OUT=args.out_dir.resolve();CACHE=args.cache_dir.resolve();SEED=args.seed
OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
def raw_file(gsm,suffix):
 candidates=list(RAW.glob(gsm+'*_'+suffix))+list((RAW/gsm).glob(gsm+'*_'+suffix))
 if len(candidates)!=1:raise ValueError(f'Exactly one raw file required: {gsm} {suffix}: {candidates}')
 return candidates[0]
def file_sha(path):
 h=hashlib.sha256()
 with path.open('rb') as handle:
  for block in iter(lambda:handle.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()
manifest=pd.read_csv(args.sample_manifest,sep='\t',dtype=str,keep_default_na=False)
afrows=manifest[(manifest.series=='GSE230808')&(manifest.tissue=='AF')]
if len(afrows)!=8 or afrows.donor_id.nunique()!=8:raise ValueError('Expected exactly eight independently verified paired AF libraries')
MARKERS={'NP_chondrocyte_like':['KRT8','KRT18','KRT19','CA12','SOX9','ACAN','COL2A1'],
 'stromal':['COL1A1','COL1A2','COL3A1','DCN','LUM','COL6A1','COL6A2'],
 'immune':['PTPRC','LST1','TYROBP','FCER1G','CD3D','CD79A','LYZ'],
 'endothelial':['PECAM1','VWF','CDH5','EMCN','KDR'],
 'mural':['RGS5','CSPG4','MCAM','MYH11','NOTCH3'],
 'erythroid':['HBB','HBA1','HBA2','ALAS2']}
def write(x,name):x.to_csv(OUT/name,sep='\t',index=False)
start=dt.datetime.now(dt.timezone.utc).isoformat()
dest=CACHE/'af_identity.h5ad'
if dest.exists():
 raise RuntimeError('AF identity output already exists; use cached mode with its real metadata or a fresh raw workspace')
af=afrows.to_dict('records')
objects=[];sample_records=[];qc_rows=[]
for d in af:
 gsm=d['gsm'];donor=d['donor_id'];barcode_file=raw_file(gsm,'barcodes.tsv.gz');prefix=barcode_file.name.replace('barcodes.tsv.gz','')
 raw_hashes={suffix:file_sha(raw_file(gsm,suffix)) for suffix in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']}
 cache=CACHE/f'af_{gsm}_QC_scrublet.h5ad'
 if cache.exists():
  a=ad.read_h5ad(cache); record=json.loads(cache.with_suffix('.json').read_text())
  if record.get('scrublet_status')!='COMPLETE' or record.get('seed')!=SEED or record.get('raw_sha256')!=raw_hashes:raise RuntimeError('Incomplete, different-input or different-seed AF cache; use a fresh workspace')
  print('REUSE_SAMPLE_QC',gsm,a.shape,flush=True)
 else:
  a=sc.read_10x_mtx(barcode_file.parent,var_names='gene_symbols',make_unique=True,prefix=prefix,cache=False)
  assert sparse.issparse(a.X)
  a.X=a.X.tocsr();assert np.isfinite(a.X.data).all() and np.all(a.X.data>=0) and np.all(a.X.data==np.floor(a.X.data))
  a.obs_names=[f'{gsm}:{barcode}' for barcode in a.obs_names]
  a.obs['gsm']=gsm;a.obs['donor_id']=donor;a.obs['condition']='diseased';a.obs['grade']=d['thompson_grade']
  a.obs['geo_series']='GSE230808';a.obs['tissue']='annulus_fibrosus';a.obs['disc']=d['disc']
  counts=np.asarray(a.X.sum(axis=1)).ravel();genes=np.asarray((a.X>0).sum(axis=1)).ravel()
  mt=a.var_names.str.upper().str.startswith('MT-');mito=np.asarray(a.X[:,mt].sum(axis=1)).ravel()
  a.obs['n_counts']=counts;a.obs['n_features']=genes;a.obs['pct_mito']=np.divide(mito,counts,out=np.zeros_like(mito),where=counts>0)*100
  a.obs['pass_qc']=(genes>=300)&(genes<=9000)&(a.obs.pct_mito<=25)&(counts>0)
  raw_n=a.n_obs
  a.obs[['gsm','donor_id','n_counts','n_features','pct_mito','pass_qc']].to_csv(OUT/f'af_{gsm}_raw_QC.tsv.gz',sep='\t',index_label='cell_id')
  a=a[a.obs.pass_qc].copy()
  record={'gsm':gsm,'donor_id':donor,'raw_cells':raw_n,'basic_QC_cells':a.n_obs,
          'expected_doublet_rate':.06,'n_prin_comps':30,'sim_doublet_ratio':2.,'seed':SEED,'raw_sha256':raw_hashes}
  try:
   sc.pp.scrublet(a,expected_doublet_rate=.06,random_state=SEED,n_prin_comps=30,verbose=True)
   record.update(scrublet_status='COMPLETE',scrublet_threshold=float(a.uns['scrublet']['threshold']),predicted_doublets=int(a.obs.predicted_doublet.sum()))
  except Exception as e:
   traceback.print_exc();a.obs['doublet_score']=np.nan;a.obs['predicted_doublet']=False
   record.update(scrublet_status='PARTIAL',scrublet_error=repr(e),predicted_doublets=0)
   cache.with_suffix('.json').write_text(json.dumps(record,indent=2))
   raise RuntimeError('AF Scrublet failed; no fallback singlet population is accepted') from e
  a.write_h5ad(cache,compression='gzip');cache.with_suffix('.json').write_text(json.dumps(record,indent=2))
 a.obs.to_csv(OUT/f'af_{gsm}_basic_QC_doublet_audit.tsv.gz',sep='\t',index_label='cell_id')
 record['retained_singlet_cells']=int((~a.obs.predicted_doublet).sum());sample_records.append(record)
 objects.append(a[~a.obs.predicted_doublet].copy())
 print('AF_SAMPLE_COMPLETE',record,flush=True)
write(pd.DataFrame(sample_records),'AF_QC_DOUBLET_SUMMARY.tsv')
a=ad.concat(objects,join='inner',merge='same');a.obs_names.name='cell_id'
assert a.obs_names.is_unique and a.var_names.is_unique
a.layers['counts']=a.X.copy();sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat',batch_key='gsm')
b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10)
sc.pp.pca(b,n_comps=30,random_state=SEED);a.obsm['X_pca']=b.obsm['X_pca']
sc.pp.neighbors(a,n_neighbors=15,n_pcs=30,random_state=SEED)
sc.tl.leiden(a,resolution=.5,random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
sc.tl.umap(a,random_state=SEED)
rows=[];markers=[];top=[];globalmean=np.asarray(a.X.mean(axis=0)).ravel()
for cluster in a.obs.leiden.cat.categories:
 mask=a.obs.leiden==cluster;x=a.X[mask];mean=np.asarray(x.mean(axis=0)).ravel();frac=np.asarray((x>0).mean(axis=0)).ravel()
 r={'cluster':str(cluster),'n_cells':int(mask.sum()),'n_donors':int(a.obs.loc[mask,'donor_id'].nunique()),
    'largest_donor_fraction':float(a.obs.loc[mask,'donor_id'].value_counts(normalize=True).max()),
    'median_counts':float(a.obs.loc[mask,'n_counts'].median()),'median_mito':float(a.obs.loc[mask,'pct_mito'].median())}
 for group,genes in MARKERS.items():
  ix=a.var_names.get_indexer([g for g in genes if g in a.var_names]);r[group+'_mean']=float(mean[ix].mean());r[group+'_n_markers']=int((mean[ix]>.35).sum())
  for g in genes:
   j=a.var_names.get_indexer([g])[0];markers.append({'cluster':str(cluster),'marker_class':group,'gene':g,'mean_logCP10K':float(mean[j]) if j>=0 else np.nan,'fraction':float(frac[j]) if j>=0 else np.nan,'available':j>=0})
 # One predeclared classification rule: a strong coherent non-stromal marker set
 # takes priority. Chondrocytic markers do not establish an NP identity in AF tissue.
 non=[g for g in ['immune','endothelial','mural','erythroid'] if r[g+'_n_markers']>=3 and r[g+'_mean']>.5 and r[g+'_mean']>2*r['NP_chondrocyte_like_mean']]
 if non:
  g=max(non,key=lambda g:r[g+'_mean']);label=g if r[g+'_mean']>=r['stromal_mean'] else 'mixed_ambiguous'
 elif r['stromal_n_markers']>=3:label='AF_tissue_stromal'
 elif r['NP_chondrocyte_like_n_markers']>=2:label='AF_chondrocytic_marker_positive'
 else:label='mixed_ambiguous'
 r['identity']=label;rows.append(r)
 for rank,j in enumerate(np.argsort(mean-globalmean)[-50:][::-1],1):
  top.append({'cluster':str(cluster),'rank':rank,'gene':a.var_names[j],'mean_logCP10K':float(mean[j]),'fraction':float(frac[j]),'difference_from_all':float(mean[j]-globalmean[j])})
mapping={r['cluster']:r['identity'] for r in rows};a.obs['identity']=a.obs.leiden.astype(str).map(mapping)
a.obs['broad_stromal']=a.obs.identity.isin(['AF_tissue_stromal','AF_chondrocytic_marker_positive'])
a.obs['high_confidence_NP']=False
# Candidate-free fixed lineage module summaries help subsequent sender definition.
for group,genes in MARKERS.items():
 ix=a.var_names.get_indexer([g for g in genes if g in a.var_names]);a.obs[group+'_marker_mean_logCP10K']=np.asarray(a.X[:,ix].mean(axis=1)).ravel()
write(pd.DataFrame(rows),'AF_CLUSTER_IDENTITY.tsv');write(pd.DataFrame(markers),'AF_CLUSTER_MARKERS.tsv');write(pd.DataFrame(top),'AF_CLUSTER_TOP50.tsv')
write(a.obs.reset_index(),'AF_CELL_METADATA.tsv.gz')
write(a.obs.groupby(['gsm','donor_id','identity'],observed=True).size().reset_index(name='n_cells'),'AF_DONOR_IDENTITY_COUNTS.tsv')
write(pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP1','UMAP2']).reset_index(),'AF_UMAP_COORDINATES.tsv.gz')
a.X=a.layers['counts'].copy();del a.layers['counts'];assert np.all(a.X.data==np.floor(a.X.data))
a.uns['af_audit_provenance']={'source_manifest_sha256':hashlib.sha256(args.sample_manifest.read_bytes()).hexdigest(),'seed':SEED,'analysis_type':'one fixed AF annotation pass','scrublet_expected_rate':.06,
 'source':'GSE230808 eight exact AF/NP donor pairs','not_NP_label':'AF chondrocytic markers do not imply NP identity',
 'qc':'300<=genes<=9000;mitochondrial_percent<=25;positive_counts','leiden_resolution':.5,'HVG':2000,
 'caveat':'Marker-based working labels; low doublet call rates or donor-dominant clusters remain limitations.'}
a.write_h5ad(dest,compression='gzip')
note=f'''# AF input and cell identity audit\n\nEight exact patient/disc-matched diseased AF samples were newly obtained from official GSE230808 supplementary URLs. No existing NP or external count file was re-downloaded.\n\nBasic QC was fixed at 300–9000 detected genes, mitochondrial percentage ≤25%, positive counts. Official Scanpy Scrublet was run once per sample (expected rate 0.06, simulation ratio 2, 30 PCs, seed {SEED}, automatic thresholds). Predicted doublets were excluded from the analysis object while all raw/basic-QC audit records were retained.\n\nA single AF-specific normalize/log1p, batch-aware 2,000-HVG, scaled 30-PC, 15-neighbour, Leiden resolution 0.5 / two-iteration graph and UMAP was run. No resolution search occurred. Fixed marker sets are written in the script and full marker means/detection fractions are exported. AF cells with cartilage markers remain labelled by AF tissue origin; they are not renamed NP. Coherent immune/endothelial/mural/erythroid dominance takes precedence; uncertain mixed cells remain in the object but are excluded from the broad-stromal primary set. The rules are in the script. These are working labels, not lineage or transdifferentiation proof.\n\nRetained cells: {a.n_obs}; genes: {a.n_vars}; broad-stromal cells: {int(a.obs.broad_stromal.sum())}; donors: {a.obs.donor_id.nunique()}.\n\nRaw integer sparse counts are stored in `work/af_identity.h5ad` X; cell/sample/donor IDs use the existing NP convention `GSM:barcode`, `gsm`, `donor_id`. Per-cell QC, doublets, donor composition, cluster marker summaries and UMAP coordinates are exported in `00_audit`. Full count objects remain on the server.\n'''
(OUT/'AF_IDENTITY_AUDIT.md').write_text(note)
(OUT/'af_IDENTITY_RUN.json').write_text(json.dumps({'started_utc':start,'completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'shape':list(a.shape),'broad_stromal':int(a.obs.broad_stromal.sum()),'donors':a.obs.donor_id.nunique(),'seed':SEED,'doublets':int(sum(d['predicted_doublets'] for d in sample_records)), 'scrublet_all_libraries_complete':all(d['scrublet_status']=='COMPLETE' for d in sample_records), 'manifest_sha256':hashlib.sha256(args.sample_manifest.read_bytes()).hexdigest(), 'versions':{'scanpy':sc.__version__,'anndata':ad.__version__,'numpy':np.__version__,'pandas':pd.__version__}},indent=2))
print('AF_IDENTITY_COMPLETE',a.shape,flush=True)
print(pd.DataFrame(rows).to_string(index=False),flush=True)
