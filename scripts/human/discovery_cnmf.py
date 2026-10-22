"""One fixed cNMF run; retain all three K choices for a documented main-K review."""
import json
import os
import subprocess
import sys
from pathlib import Path

for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[key]='1'
os.environ['MPLBACKEND']='Agg'
import anndata as ad
import numpy as np
import pandas as pd
from cnmf import cNMF
from cnmf.cnmf import load_df_from_npz
from threadpoolctl import threadpool_limits

ROOT=Path(os.environ['IVDD_DISCOVERY_ROOT']).resolve()
OUT=ROOT/'04_cNMF'
NAME='IVDD_NP_seed20261012'
obj=cNMF(output_dir=str(OUT),name=NAME)
PARAMS=dict(K_candidates=[5,7,9],n_initializations_per_K=20,seed=20261012,num_highvar_genes=2000,max_fitting_cells=12000,max_NMF_iter=1000,beta_loss='frobenius',density_threshold=.5,local_neighborhood_size=.3,workers=4,threads_per_worker=1,normalization='official_cNMF; no batch correction; counts input',main_K_rule='Inspect official stability/error and descriptive marker composition; prefer parsimonious stable solution without reference to THBS1 or external outcomes. No new K or parameter search.')
(OUT/'parameters_pre_run.json').write_text(json.dumps(PARAMS,indent=2))

def factor(worker):
    obj.factorize(worker_i=worker,total_workers=4,skip_completed_runs=True)

def export_candidates():
    candidate_genes=[]
    for k in [5,7,9]:
        usage,scores,tpm,top=obj.load_results(K=k,density_threshold=.5,n_top_genes=50)
        assert scores.shape[1]==k and usage.shape[1]==k, 'Unexpected official load_results orientation'
        usage.to_csv(OUT/f'candidate_K{k}_cell_usage.tsv.gz',sep='\t')
        scores.to_csv(OUT/f'candidate_K{k}_gene_spectra_score.tsv.gz',sep='\t')
        for program in scores.columns:
            for rank,(gene,score) in enumerate(scores[program].sort_values(ascending=False).head(50).items(),1):
                candidate_genes.append(dict(K=k,program=str(program),rank=rank,gene=gene,gene_spectra_score=score))
    target=OUT/'all_candidate_K_top50_genes.tsv'
    if target.exists() and len(pd.read_csv(target,sep='\t'))!=1050:
        target.rename(OUT/'superseded_orientation_error_gene_export.tsv')
    pd.DataFrame(candidate_genes).to_csv(target,sep='\t',index=False)

def main():
    if (OUT/'candidate_Ks_complete.json').exists():
        export_candidates();print('RESUME_CANDIDATES_COMPLETE_EXPORT_VALIDATED');return
    with threadpool_limits(limits=1):
        if not Path(obj.paths['nmf_replicate_parameters']).exists():
            obj.prepare(counts_fn=str(ROOT/'03_discovery_NP/cnmf_input_counts.h5ad'),components=[5,7,9],n_iter=20,seed=20261012,num_highvar_genes=2000,max_NMF_iter=1000)
        print('PREPARED',flush=True)
        procs=[]
        for i in range(4):
            log=(ROOT/'logs'/f'04_factorize_worker{i}.log').open('a')
            procs.append((subprocess.Popen([sys.executable,'-u',str(Path(__file__).resolve()),'--worker',str(i)],stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy()),log))
        exits=[p.wait() for p,_ in procs]
        for _,f in procs: f.close()
        if any(exits): raise RuntimeError('cNMF worker failure '+str(exits))
        obj.combine(skip_missing_files=False)
        obj.k_selection_plot(close_fig=True)
        stats=load_df_from_npz(obj.paths['k_selection_stats'])
        stats['n_initializations']=20;stats['seed']=20261012;stats['n_HVG']=2000
        stats.to_csv(OUT/'K_selection_and_parameters.tsv',sep='\t',index=False)
        for k in [5,7,9]:
            print('CONSENSUS',k,flush=True)
            obj.consensus(k=k,density_threshold=.5,local_neighborhood_size=.3,show_clustering=True,close_clustergram_fig=True,build_ref=False)
        export_candidates()
        (OUT/'candidate_Ks_complete.json').write_text(json.dumps({'status':'main_K_review_required','parameters':PARAMS},indent=2))
        print('ALL_CANDIDATE_KS_COMPLETE',stats.to_string(index=False),flush=True)

if __name__=='__main__':
    if '--worker' in sys.argv: factor(int(sys.argv[sys.argv.index('--worker')+1]))
    else: main()
