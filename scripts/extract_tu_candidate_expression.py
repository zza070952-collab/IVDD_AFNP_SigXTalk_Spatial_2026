"""Summarize measured candidate expression by true donor and working cell identity."""
import argparse
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tu-sparse-root',type=Path,required=True);p.add_argument('--program-mapping',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
 p.add_argument('--genes',nargs='+',default='FN1 ITGB1 PRSS23 PLOD2 CCN2 THBS1 LOX TGFBI IGFBP3 GJA1'.split());a=p.parse_args()
 mapping=pd.read_csv(a.program_mapping,sep='\t').fillna({'mapped_gene':''});rows=[]
 folders=sorted(a.tu_sparse_root.glob('GSM*'));assert len(folders)==8,'Expected actual eight-donor Tu cache'
 for folder in folders:
  x=sparse.load_npz(folder/'expression_gene_by_cell.npz').tocsr();features=pd.read_csv(folder/'genes.tsv',sep='\t').gene.tolist();index={g:i for i,g in enumerate(features)}
  meta=pd.read_csv(folder/'metadata.tsv.gz',sep='\t');passed=meta.qc_pass.astype(str).str.lower().isin(['true','1']);assert x.shape==(len(features),len(meta))
  for gene in a.genes:
   feature=gene
   if gene not in index:
    allowed=mapping[(mapping.gsm==folder.name)&(mapping.original_gene==gene)].mapped_gene
    present=[v for v in set(allowed) if v in index]
    if len(present)==1:feature=present[0]
   found=feature in index;v=x.getrow(index[feature]).toarray().ravel() if found else None
   for identity,g in meta.loc[passed].groupby('identity'):
    take=g.index.to_numpy();rows.append({'gsm':folder.name,'donor_id':g.donor_id.iloc[0],'identity':identity,'gene':gene,'matched_released_feature':feature if found else '', 'n_cells':len(take),'feature_present':found,'mean_log1p_published_expression':float(np.log1p(v[take]).mean()) if found else np.nan,'fraction_detected':float((v[take]>0).mean()) if found else np.nan,'evidence_type':'published normalized expression; descriptive donor-level expression; clinical grade unresolved'})
 a.out.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(rows).to_csv(a.out,sep='\t',index=False,na_rep='NA');print('Candidate donor-identity rows:',len(rows))
if __name__=='__main__':main()
