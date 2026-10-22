"""Sparse integer counts, fixed QC, descriptive lineage audit, balanced fitting cells."""
import argparse
import gzip
import json
import os
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('MKL_NUM_THREADS','2')
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits

ROOT=Path(os.environ['IVDD_DISCOVERY_ROOT']).resolve()
SEED=20261012
MARKERS={
 'disc_matrix':['ACAN','COL2A1','SOX9','COMP','FMOD','CHAD','PRELP','CILP','DCN','LUM','COL1A1','COL1A2','COL3A1','POSTN'],
 'immune':['PTPRC','CD3D','CD3E','TRBC1','TRBC2','CD79A','MS4A1','LST1','TYROBP','FCER1G','LYZ'],
 'endothelial':['PECAM1','VWF','CDH5','EMCN','KDR','PLVAP'],
 'erythroid':['HBB','HBA1','HBA2','ALAS2'],
 'mural':['RGS5','MCAM','CSPG4','NOTCH3','MYH11'],
 'stress':['FOS','JUN','JUNB','HSPA1A','HSPA1B','DDIT3','ATF3'],
}
ALIASES={'CTGF':'CCN2','CYR61':'CCN1','WISP1':'CCN4','WISP2':'CCN5','WISP3':'CCN6'}

def read_table(path):
    return pd.read_csv(path,sep='\t',header=None,dtype=str,keep_default_na=False)

def load_sample(row):
    series,gsm=row.geo_series,row.gsm
    d=ROOT/'02_geo_downloads'/series
    if series=='GSE230808':
        matrix=next(d.glob(gsm+'*_matrix.mtx.gz'))
        features=next(d.glob(gsm+'*_features.tsv.gz'))
        barcodes=next(d.glob(gsm+'*_barcodes.tsv.gz'))
    else:
        archive=next(d.glob(gsm+'*.tar.gz'))
        base=archive.parent/archive.stem
        matrices=list(base.rglob('*matrix.mtx*'))
        if len(matrices)!=1: raise ValueError(f'Ambiguous matrix structure for {gsm}: {matrices}')
        matrix=matrices[0]
        features=next(iter(list(matrix.parent.glob('*features.tsv*'))+list(matrix.parent.glob('*genes.tsv*'))))
        barcodes=next(matrix.parent.glob('*barcodes.tsv*'))
    with gzip.open(matrix,'rb') if matrix.name.endswith('.gz') else matrix.open('rb') as f:
        x=io.mmread(f).T.tocsr()
    if x.data.size and (not np.isfinite(x.data).all() or x.data.min()<0 or not np.equal(x.data,np.floor(x.data)).all()):
        raise ValueError('Not nonnegative integer counts '+gsm)
    x=x.astype(np.float32)
    feat=read_table(features); bc=read_table(barcodes)[0].to_numpy()
    if x.shape!=(len(bc),len(feat)): raise ValueError('Dimension mismatch '+gsm)
    symbols=(feat[1] if feat.shape[1]>1 else feat[0]).map(lambda s:ALIASES.get(s,s)).to_numpy()
    symbols=np.where(symbols=='',feat[0].to_numpy(),symbols)
    unique=pd.Index(pd.unique(symbols))
    mapping=sparse.csr_matrix((np.ones(len(symbols),dtype=np.float32),(np.arange(len(symbols)),unique.get_indexer(symbols))),shape=(len(symbols),len(unique)))
    x=(x@mapping).tocsr();x.eliminate_zeros()
    var=pd.DataFrame(index=unique)
    obs=pd.DataFrame(index=[gsm+':'+b for b in bc])
    for k in ['gsm','donor_id','condition','grade','geo_series']: obs[k]=str(row[k])
    a=ad.AnnData(x,obs=obs,var=var)
    a.obs['n_counts']=np.asarray(x.sum(axis=1)).ravel()
    a.obs['n_features']=np.diff(x.indptr)
    mt=unique.str.upper().str.startswith('MT-')
    a.obs['pct_mito']=100*np.asarray(x[:,mt].sum(axis=1)).ravel()/np.maximum(a.obs.n_counts,1)
    a.obs['pass_qc']=(a.obs.n_features>=300)&(a.obs.n_features<=9000)&(a.obs.pct_mito<=25)&(a.obs.n_counts>0)
    report=dict(gsm=gsm,donor_id=row.donor_id,condition=row.condition,grade=row.grade,n_cells_original=a.n_obs,n_genes_original=len(feat),n_unique_symbols=a.n_vars,duplicate_symbols_collapsed=len(feat)-a.n_vars,n_cells_qc=int(a.obs.pass_qc.sum()),median_counts=float(a.obs.n_counts.median()),median_features=float(a.obs.n_features.median()),median_pct_mito=float(a.obs.pct_mito.median()),matrix_type='processed_nonnegative_integer_UMI_counts',doublet_assessment='not_evaluated; upper_feature_filter_not_a_doublet_detector')
    print('SAMPLE_QC',json.dumps(report),flush=True)
    return a,report


def audit_lineage(a,out):
    b=a.copy()
    sc.pp.normalize_total(b,target_sum=10000)
    sc.pp.log1p(b)
    gene_total=np.asarray(a.X.sum(axis=0)).ravel()
    b=b[:,gene_total>0].copy()
    sc.pp.highly_variable_genes(b,n_top_genes=2000,flavor='seurat')
    sc.pp.pca(b,n_comps=30,mask_var='highly_variable',svd_solver='arpack',random_state=SEED)
    labels=KMeans(n_clusters=20,n_init=10,random_state=SEED).fit_predict(b.obsm['X_pca'])
    a.obs['audit_cluster']=labels.astype(str)
    means=[];top=[]
    global_mean=np.asarray(b.X.mean(axis=0)).ravel()
    for cluster in sorted(np.unique(labels)):
        mask=labels==cluster
        mean=np.asarray(b.X[mask].mean(axis=0)).ravel()
        row=dict(cluster=str(cluster),n_cells=int(mask.sum()),n_donors=int(a.obs.loc[mask,'donor_id'].nunique()))
        for name,genes in MARKERS.items():
            present=[g for g in genes if g in b.var_names]
            row[name+'_mean_log1p']=float(mean[b.var_names.get_indexer(present)].mean()) if present else 0.
        for g in sum(MARKERS.values(),[]):
            row[g]=float(mean[b.var_names.get_loc(g)]) if g in b.var_names else 0.
        # Exclude only coherent, clearly dominant non-disc marker programs.
        immune=sum(row[g]>.5 for g in MARKERS['immune'])>=3 and row['immune_mean_log1p']>row['disc_matrix_mean_log1p']*1.5
        endo=sum(row[g]>.5 for g in MARKERS['endothelial'])>=3 and row['endothelial_mean_log1p']>row['disc_matrix_mean_log1p']*1.5
        eryth=sum(row[g]>1 for g in MARKERS['erythroid'])>=3 and row['erythroid_mean_log1p']>row['disc_matrix_mean_log1p']*2
        mural=sum(row[g]>.5 for g in ['RGS5','CSPG4','NOTCH3','MYH11'])>=3 and row['mural_mean_log1p']>row['disc_matrix_mean_log1p']*1.5
        cls='immune' if immune else 'endothelial' if endo else 'erythroid' if eryth else 'mural' if mural else 'disc_lineage_NP_like' if sum(row[g]>.3 for g in MARKERS['disc_matrix'])>=3 else 'ambiguous'
        row['annotation']=cls;row['retain_for_NP_programs']=cls=='disc_lineage_NP_like'
        means.append(row)
        score=mean-global_mean
        for rank,idx in enumerate(np.argsort(score)[-30:][::-1],1):
            top.append(dict(cluster=str(cluster),rank=rank,gene=b.var_names[idx],mean_log1p=mean[idx],difference_from_global=score[idx]))
    pd.DataFrame(means).to_csv(out/'cluster_marker_audit.tsv',sep='\t',index=False)
    pd.DataFrame(top).to_csv(out/'cluster_top_markers_descriptive.tsv',sep='\t',index=False)
    mapping={r['cluster']:r['annotation'] for r in means}
    a.obs['lineage_annotation']=a.obs.audit_cluster.map(mapping)
    a.obsm['X_pca_audit']=b.obsm['X_pca']
    print('LINEAGE_AUDIT',pd.DataFrame(means)[['cluster','n_cells','n_donors','annotation']].to_string(index=False),flush=True)
    return a


def run(series):
    out=ROOT/('03_discovery_NP' if series=='GSE230808' else '05_external_GSE251686')
    marker=out/'qc_complete.json'
    if marker.exists(): print('RESUME_QC_COMPLETE',series,flush=True);return
    manifest_path=ROOT/'01_metadata_audit/sample_manifest.tsv'
    manifest=pd.read_csv(manifest_path,sep='\t',keep_default_na=False)
    rows=manifest[(manifest.geo_series==series)&(manifest.selected_for_analysis.astype(str).str.lower()=='true')]
    intended_samples=len(rows)
    if series=='GSE251686':
        integrity=ROOT/'01_metadata_audit/external_matrix_integrity.tsv'
        if not integrity.exists(): raise RuntimeError('External source integrity audit required before QC')
        check=pd.read_csv(integrity,sep='\t')
        unusable=check[check.usable_counts.astype(str).str.lower()!='true']
        unusable.to_csv(out/'unreadable_source_samples.tsv',sep='\t',index=False)
        rows=rows[~rows.gsm.isin(unusable.gsm)]
        print('INCOMPLETE_EXTERNAL_COHORT',len(rows),'of',intended_samples,'solely source integrity exclusions',flush=True)
    arrays,reports=[],[]
    with threadpool_limits(limits=2):
        for _,row in rows.iterrows():
            a,r=load_sample(row); reports.append(r)
            a.write_h5ad(out/(row.gsm+'.original_counts.h5ad'),compression='gzip')
            arrays.append(a[a.obs.pass_qc].copy())
            manifest.loc[manifest.gsm==row.gsm,['n_cells_original','matrix_type']]=[str(a.n_obs),r['matrix_type']]
        manifest.to_csv(manifest_path,sep='\t',index=False)
        qc=pd.DataFrame(reports)
        qc.to_csv(out/'QC_by_donor.tsv',sep='\t',index=False)
        a=ad.concat(arrays,join='outer',merge='first',fill_value=0)
        a.X=a.X.tocsr();a.X.eliminate_zeros()
        a=a[:,np.asarray(a.X.sum(axis=0)).ravel()>0].copy()
        a=audit_lineage(a,out)
        a.write_h5ad(out/'all_QC_counts_with_lineage.h5ad',compression='gzip')
        a.obs.to_csv(out/'cell_metadata.tsv.gz',sep='\t')
        np_cells=a[a.obs.lineage_annotation=='disc_lineage_NP_like'].copy()
        np_cells=np_cells[:,np.asarray(np_cells.X.sum(axis=0)).ravel()>0].copy()
        np_cells.write_h5ad(out/'NP_lineage_counts.h5ad',compression='gzip')
        summary=np_cells.obs.groupby(['donor_id','gsm'],observed=True).size().rename('n_NP_lineage').reset_index()
        summary.to_csv(out/'NP_lineage_cells_by_donor.tsv',sep='\t',index=False)
        if series=='GSE230808':
            if np_cells.obs.donor_id.nunique()<2: raise RuntimeError('BLOCKED: fewer than 2 trustworthy NP donors')
            rng=np.random.default_rng(SEED)
            groups=[np.where(np_cells.obs.donor_id.to_numpy()==d)[0] for d in sorted(np_cells.obs.donor_id.unique())]
            # Water filling maximizes coverage while respecting a 12000 total cap.
            allocated=np.zeros(len(groups),dtype=int)
            remaining=min(12000,np_cells.n_obs)
            while remaining:
                eligible=np.where(allocated<np.array([len(g) for g in groups]))[0]
                for i in eligible:
                    if remaining==0: break
                    allocated[i]+=1;remaining-=1
            selected=np.sort(np.concatenate([rng.choice(g,n,replace=False) for g,n in zip(groups,allocated)]))
            fit=np_cells[selected].copy()
            fit=fit[:,np.asarray(fit.X.sum(axis=0)).ravel()>0].copy()
            fit.write_h5ad(out/'cnmf_input_counts.h5ad',compression='gzip')
            fit.obs.to_csv(out/'cnmf_fitting_cell_list.tsv',sep='\t')
        stats=dict(series=series,intended_samples=intended_samples,loaded_samples=len(rows),n_cells_original=int(qc.n_cells_original.sum()),n_cells_qc=a.n_obs,n_NP_lineage=np_cells.n_obs,n_donors=int(np_cells.obs.donor_id.nunique()),qc=dict(min_features=300,max_features=9000,max_pct_mito=25),doublets='not_evaluated',seed=SEED,ambiguous='retained_in_audit_but_excluded_from_program_fitting',fit_cap=12000)
        marker.write_text(json.dumps(stats,indent=2))
        print('QC_COMPLETE',stats,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--study',choices=['GSE230808','GSE251686'],required=True)
    run(parser.parse_args().study)
