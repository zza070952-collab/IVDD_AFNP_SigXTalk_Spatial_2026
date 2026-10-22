"""Rebuild the predeclared AF/NP input from genuine public 10x integer counts."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
import argparse,json,gzip,hashlib,datetime
from pathlib import Path
import numpy as np,pandas as pd
from scipy import io,sparse
from anndata.utils import make_index_unique
p=argparse.ArgumentParser();p.add_argument('--project-root',required=True);p.add_argument('--metadata-root',required=True);a=p.parse_args()
W=Path(a.project_root);old=Path(a.metadata_root);out=W/'02_sigxtalk_official_new/input';out.mkdir(parents=True,exist_ok=True)
seed=20261017
legacy=pd.read_csv(W/'configs/LEGACY_PREDECLARED_INPUT_CELLS.tsv.gz',sep='\t').set_index('cell_id',drop=False)
n=pd.read_csv(old/'data_cache/discovery/discovery_cell_metadata.tsv.gz',sep='\t').set_index('cell_id',drop=False)
scores=pd.read_csv(old/'data_cache/discovery/discovery_programs_UCell.tsv.gz',sep='\t').set_index('cell_id')
af=pd.read_csv(old/'data_cache/legacy_bridge/00_audit/AF_CELL_METADATA.tsv.gz',sep='\t').set_index('cell_id',drop=False)
n['GEP2']=scores.loc[n.index,'GEP2_top25_UCell'].to_numpy()
def bo(x):return x.astype(str).str.lower().isin(['true','1'])
rng=np.random.default_rng(seed);parts=[];elig=[]
for donor,z0 in n.groupby('donor_id',sort=True,observed=True):
    z=z0[bo(z0.high_confidence_NP)];q=z.GEP2.quantile([.25,.75]);lo=z[z.GEP2<=q.iloc[0]];hi=z[z.GEP2>=q.iloc[1]]
    okay=len(lo)>=20 and len(hi)>=20 and q.iloc[1]>q.iloc[0]
    elig.append(dict(donor_id=donor,n_NP=len(z),n_low=len(lo),n_high=len(hi),q25=q.iloc[0],q75=q.iloc[1],eligible=bool(okay)))
    if okay:
        for group,z1 in [('NP_GEP2_low',lo),('NP_GEP2_high',hi)]:
            cells=sorted(z1.index);sel=sorted(rng.choice(cells,min(150,len(cells)),replace=False));tmp=z1.loc[sel].copy();tmp['bridge_group']=group;tmp['tissue']='NP';parts.append(tmp)
for donor in [x['donor_id'] for x in elig if x['eligible']]:
    z=af[(af.donor_id==donor)&bo(af.broad_stromal)];assert len(z)>=20,donor
    cells=sorted(z.index);sel=sorted(rng.choice(cells,min(150,len(cells)),replace=False));tmp=z.loc[sel].copy();tmp['bridge_group']='AF_stromal';tmp['tissue']='AF';parts.append(tmp)
meta=pd.concat(parts,axis=0);assert meta.index.is_unique
pd.DataFrame(elig).to_csv(out/'RECOMPUTED_DONOR_ELIGIBILITY.tsv',sep='\t',index=False)
selection={'seed':seed,'n_selected':len(meta),'eligible_donors':sorted(meta.donor_id.unique()),'legacy_cells':len(legacy),'same_cell_set':set(meta.index)==set(legacy.index),'same_cell_order':list(meta.index)==list(legacy.index),'selection_criteria':'within-donor marker-positive NP UCell Q1/Q4; min20 and q75>q25; cap150; broad-stromal AF paired donors; seeded sorted-barcode sampling'}
if set(meta.index)==set(legacy.index):assert (meta.bridge_group==legacy.loc[meta.index,'bridge_group']).all()
print('RECOMPUTED_SELECTION',json.dumps(selection),flush=True)
(out/'BARCODE_SELECTION_PROOF.json').write_text(json.dumps(selection,indent=2))
source=pd.read_csv(W/'configs/GSE230808_SOURCE_MANIFEST.tsv',sep='\t').set_index('gsm')
matrices={};genes=None;audits=[];hashes=[]
for gsm in sorted(meta.gsm.unique()):
    d=W/'data_cache/GSE230808'/gsm;f=next(d.glob('*_features.tsv.gz'));b=next(d.glob('*_barcodes.tsv.gz'));m=next(d.glob('*_matrix.mtx.gz'))
    ft=pd.read_csv(f,sep='\t',header=None);bc=pd.read_csv(b,sep='\t',header=None)[0].astype(str)
    with gzip.open(m,'rb') as stream:x=io.mmread(stream).tocsc()
    assert x.shape==(len(ft),len(bc)) and np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
    keep=(ft[2].astype(str)=='Gene Expression').to_numpy() if ft.shape[1]>=3 else np.ones(len(ft),bool)
    gg=make_index_unique(pd.Index(ft.loc[keep,1].astype(str)));x=x[keep].tocsc()
    if genes is None:genes=gg
    assert list(gg)==list(genes),'different gene feature ordering requires explicit alignment'
    cells=meta[meta.gsm==gsm].copy();assert len(cells)>0
    assert (cells.donor_id==source.loc[gsm,'donor_id']).all()
    ids=pd.Index(gsm+':'+bc);idx=ids.get_indexer(cells.index);assert np.all(idx>=0),'unresolved barcode'
    y=x[:,idx].astype(np.int64);matrices[gsm]=(y,list(cells.index))
    sums=np.asarray(y.sum(0)).ravel();nfeat=np.asarray((y>0).sum(0)).ravel()
    for i,cell in enumerate(cells.index):
        r=cells.loc[cell];expected=float(r.n_counts);expectedgenes=int(r.n_features)
        audits.append(dict(cell_id=cell,gsm=gsm,donor_id=r.donor_id,tissue=r.tissue,bridge_group=r.bridge_group,barcode=cell.split(':',1)[1],raw_barcode_found=True,raw_total_UMI=int(sums[i]),cached_total_UMI=expected,total_UMI_equal=bool(sums[i]==expected),raw_detected_genes=int(nfeat[i]),cached_detected_genes=expectedgenes,detected_genes_equal=bool(nfeat[i]==expectedgenes),source_count_file=str(m),metadata_source='recovered barcode-resolved frozen QC/identity and UCell scores',GEP2_top25_UCell=float(r.GEP2) if pd.notna(r.get('GEP2')) else np.nan))
    for path in [m,f,b]:
        h=hashlib.sha256()
        with path.open('rb') as st:
            for chunk in iter(lambda:st.read(4*1024*1024),b''):h.update(chunk)
        hashes.append(dict(accession=gsm,dataset='GSE230808',donor_id=source.loc[gsm,'donor_id'],tissue=source.loc[gsm,'tissue'],path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest(),source=source.loc[gsm,'matrix_url'].rsplit('/',1)[0]+'/'+path.name,reused_or_downloaded='downloaded_this_run'))
    print('RAW_MATRIX_VERIFIED',gsm,x.shape,'selected',y.shape,flush=True)
    del x
# Restore the independently regenerated seeded cell order rather than concatenating by GSM.
by_gsm=[];column_ids=[]
for gsm in sorted(matrices):y,ids=matrices[gsm];by_gsm.append(y);column_ids.extend(ids)
x=sparse.hstack(by_gsm,format='csc');order=pd.Index(column_ids).get_indexer(meta.index);assert np.all(order>=0);x=x[:,order]
with gzip.open(out/'counts.mtx.gz','wb') as f:io.mmwrite(f,x,field='integer')
pd.Series(genes).to_csv(out/'genes.tsv',sep='\t',index=False,header=False)
meta['GEP2_top25_UCell']=meta.get('GEP2');meta.to_csv(out/'meta.tsv',sep='\t',index=False)
aud=pd.DataFrame(audits).set_index('cell_id').loc[meta.index].reset_index();aud.to_csv(W/'02_sigxtalk_official_new/NEW_INPUT_CELL_DONOR_BARCODE_AUDIT.tsv',sep='\t',index=False)
assert aud.raw_barcode_found.all()
# Changed raw gene universes are explicitly retained, never reverse engineered from old outputs.
summary=dict(selection,shape_genes_by_cells=list(x.shape),nnz=int(x.nnz),raw_UMI_equal_to_cached=int(aud.total_UMI_equal.sum()),raw_features_equal_to_cached=int(aud.detected_genes_equal.sum()),count_source='GEO 10x GSE230808 integer raw matrices',gene_universe='all common Gene Expression features made unique; no legacy processed gene filtering',completed_at=datetime.datetime.now().astimezone().isoformat())
(out/'NEW_RAW_INPUT_PROOF.json').write_text(json.dumps(summary,indent=2));pd.DataFrame(hashes).to_csv(W/'00_audit/NEW_RAW_INPUT_FILES.tsv',sep='\t',index=False)
meta.groupby(['donor_id','gsm','bridge_group'],observed=True).size().rename('n_cells').reset_index().to_csv(out/'NEW_DONOR_GROUP_COUNTS.tsv',sep='\t',index=False)
print('TRUE_RAW_INPUT_COMPLETE',json.dumps(summary),flush=True)
