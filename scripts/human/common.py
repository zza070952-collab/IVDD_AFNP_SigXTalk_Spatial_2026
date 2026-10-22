from pathlib import Path
import hashlib,json,os
import numpy as np
import pandas as pd
ROOT=Path(os.environ['IVDD_DISCOVERY_REBUILD_ROOT']).resolve()
OLD=Path(os.environ['IVDD_DISCOVERY_ROOT']).resolve()
SEED=20261016
COHORTS={'discovery':'03_discovery_NP','external':'05_external_GSE251686'}
MARKERS={'NP_chondrocyte_like':['KRT8','KRT18','KRT19','CA12','SOX9','ACAN','COL2A1'],
 'stromal':['COL1A1','COL1A2','COL3A1','DCN','LUM','COL6A1','COL6A2'],
 'immune':['PTPRC','LST1','TYROBP','FCER1G','CD3D','CD79A','LYZ'],
 'endothelial':['PECAM1','VWF','CDH5','EMCN','KDR'],
 'mural':['RGS5','CSPG4','MCAM','MYH11','NOTCH3'],
 'erythroid':['HBB','HBA1','HBA2','ALAS2']}
ALIASES={'CTGF':'CCN2','CYR61':'CCN1','WISP1':'CCN4','WISP2':'CCN5','WISP3':'CCN6'}
def read(p):return pd.read_csv(ROOT/p,sep='\t')
def write(df,p):
 p=ROOT/p;p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,sep='\t',index=False)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(x,p): (ROOT/p).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
