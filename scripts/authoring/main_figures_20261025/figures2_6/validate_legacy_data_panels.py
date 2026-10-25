#!/usr/bin/env python3
"""Compare legacy data artists in memory with the new source-level composition.

Only the selected original data-panel drawing statements are executed. Original
main(), save(), manuscript, supplement and source-data export statements never run.
The extraction is not a sandbox: use trusted local code only and explicitly opt in
with --trust-legacy-code. Nothing is downloaded or run on import.
"""
import argparse
import ast
import hashlib
import json
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import replot_figures2_6 as new
from _io_guard import independent_output, reject_output_links


def between(text,start,end):
    return textwrap.dedent(text[text.index(start):text.index(end,text.index(start))])


def legacy_helpers(text,names,env):
    main=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='main')
    for name in names:
        node=next(n for n in main.body if isinstance(n,ast.FunctionDef) and n.name==name)
        exec(compile(textwrap.dedent(ast.get_source_segment(text,node)),'<legacy-data-helper>','exec'),env)


def equal(a,b,path=''):
    if isinstance(a,dict):
        if a.keys()!=b.keys(): return [path+': dict fields differ']
        return sum((equal(a[k],b[k],path+'/'+str(k)) for k in a),[])
    if isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):return [path+': lengths differ']
        return sum((equal(x,y,path+'/'+str(i)) for i,(x,y) in enumerate(zip(a,b))),[])
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        return [] if np.isclose(a,b,rtol=0,atol=1e-12,equal_nan=True) else [path+': numeric differs']
    return [] if a==b else [path+': value differs']


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True,help='Local audit JSON outside source and legacy-code trees.')
    ap.add_argument('--legacy-code-root',type=Path,help='Directory containing the two trusted legacy plotting scripts; defaults to SOURCE_ROOT/07_public_code_delta.')
    ap.add_argument('--trust-legacy-code',action='store_true',help='Required acknowledgement: selected source statements execute with current Python privileges; inspect and trust the local scripts first.')
    a=ap.parse_args()
    if not a.trust_legacy_code:ap.error('--trust-legacy-code is required; this oracle executes selected trusted local statements and is not a sandbox.')
    src,out=independent_output(a.source_root,a.output.parent)
    legacy=(a.legacy_code_root or src/'07_public_code_delta').resolve()
    independent_output(legacy,out);reject_output_links(out)
    output=a.output.resolve()
    if output.parent!=out:raise ValueError('Output file resolves outside the independent output directory.')
    inputs=[src/'03_tables/human'/(k+'_source.tsv') for k in ['Figure1B','Figure1C','Figure1D','Figure1E']]
    inputs.extend(src/'03_tables/experimental/internal/source_unresolved'/(k+'.tsv') for k in ['QPCR_BATCH_DELTACt','QPCR_CONTRASTS_AND_CI','QPCR_BATCH_CONTRASTS','CCK8_BATCH_SUMMARY'])
    inputs.extend(legacy/k for k in ['build_human_figures.py','plot_experiments.py'])
    input_hashes={path:new.sha(path) for path in inputs}
    new.setup();receipts=[]
    source={key:pd.read_csv(src/'03_tables/human'/(key+'_source.tsv'),sep='\t') for key in ['Figure1B','Figure1C','Figure1D','Figure1E']}
    oldpath=legacy/'build_human_figures.py';text=oldpath.read_text(encoding='utf8')
    old=plt.figure(figsize=(180/25.4,230/25.4));env=dict(vars(new));env['fig']=old
    legacy_helpers(text,['panel','bubble','cbar','dot_legend'],env)
    env['bdata']=source['Figure1B']
    exec(compile(between(text,"    ax=panel(fig,[.16,.555",'    tu=cnt['),'<legacy-Figure1B>','exec'),env)
    env['tu']=source['Figure1C'].pivot(index='donor_id',columns='identity',values='n_cells').reindex(index=new.TUD,columns=new.IDS).fillna(0)
    exec(compile(between(text,"    ax=panel(fig,[.62,.555",'    md=m.groupby'),'<legacy-Figure1C>','exec'),env)
    env['ddata']=source['Figure1D'];env['markers']=new.MARKERS
    exec(compile(between(text,"    ax=panel(fig,[.20,.245",'    edata=ds.reset_index'),'<legacy-Figure1D>','exec'),env)
    env['ds']=source['Figure1E'].pivot(index='donor_id',columns='identity',values='n_cells').reindex(new.DONS)
    exec(compile(between(text,"    ax=panel(fig,[.17,.035","    save(fig,'Figure1_Cohort"),'<legacy-Figure1E>','exec'),env)
    replacement=new.draw_human(source)
    lhs=new.artist_fingerprint(old);rhs=new.artist_fingerprint(replacement);diff=equal(lhs,rhs)
    receipts.append({'old_figure':'Figure1 data panels B/C/D/E','new_figure':'Figure2 A/B/C/D','comparison':'data-coordinate artists: arrays, points, bars, error segments, colors, markers, line styles, axis limits and labels','axes_compared':len(lhs),'numeric_tolerance':1e-12,'differences':diff,'status':'PASS' if not diff else 'FAIL','legacy_code_sha256':new.sha(oldpath)})
    plt.close(old);plt.close(replacement)
    s=src/'03_tables/experimental/internal/source_unresolved'
    source={key:pd.read_csv(s/(key+'.tsv'),sep='\t') for key in ['QPCR_BATCH_DELTACt','QPCR_CONTRASTS_AND_CI','QPCR_BATCH_CONTRASTS','CCK8_BATCH_SUMMARY']}
    oldpath=legacy/'plot_experiments.py';text=oldpath.read_text(encoding='utf8')
    old=plt.figure(figsize=(180/25.4,238/25.4));gs=old.add_gridspec(4,2,height_ratios=[.65,1.15,1.15,1.3],left=.11,right=.97,top=.955,bottom=.085,hspace=.65,wspace=.58)
    env=dict(vars(new));env.update({'fig':old,'gs':gs,'d':source['QPCR_BATCH_DELTACt'],'c':source['QPCR_CONTRASTS_AND_CI'],'bc':source['QPCR_BATCH_CONTRASTS'],'cc':source['CCK8_BATCH_SUMMARY'],'MARKERS':new.BATCH_MARKERS})
    node=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='label')
    exec(compile(ast.get_source_segment(text,node),'<legacy-label>','exec'),env)
    exec(compile(between(text," ax=fig.add_subplot(gs[1,0]);label(ax,'B'"," save(fig,fdir/'Figure5_Cell"),'<legacy-Figure5B-to-G>','exec'),env)
    replacement=new.draw_experiments(source)
    lhs=new.artist_fingerprint(old);rhs=new.artist_fingerprint(replacement);diff=equal(lhs,rhs)
    receipts.append({'old_figure':'Figure5 data panels B/C/D/E/F/G','new_figure':'Figure6 A/B/C/D/E/F','comparison':'data-coordinate artists: arrays, points, error segments, colors, markers, line styles, axis limits and labels','axes_compared':len(lhs),'numeric_tolerance':1e-12,'differences':diff,'status':'PASS' if not diff else 'FAIL','legacy_code_sha256':new.sha(oldpath)})
    plt.close(old);plt.close(replacement)
    if any(new.sha(path)!=expected for path,expected in input_hashes.items()):raise RuntimeError('A source table or trusted legacy script changed during validation.')
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps({'status':'PASS' if all(r['status']=='PASS' for r in receipts) else 'FAIL','checks':receipts,'legacy_whole_main_executed':False,'legacy_export_functions_executed':False,'analysis_or_inference_reexecuted':False,'all_inputs_byte_unchanged':True},indent=2),encoding='utf8')
    print(output.read_text(encoding='utf8'))
    if any(r['status']!='PASS' for r in receipts):raise SystemExit(1)


if __name__=='__main__':main()
