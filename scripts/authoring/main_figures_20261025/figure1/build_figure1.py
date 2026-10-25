"""Original study-design schematic from explicit, audited factual inputs.

No assay measurements, inferred embeddings or result heatmaps are generated.
Mini matrices, cells and spatial points are geometry-only illustrations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import (Circle, Ellipse, Rectangle, FancyBboxPatch,
                               FancyArrowPatch, PathPatch, Arc)
from matplotlib.path import Path as MplPath
from matplotlib import font_manager
import numpy as np

P = {'ink':'#26353E','muted':'#536671','edge':'#BDCBD0','np':'#267B7C',
     'np_light':'#E7F2EF','af':'#BF723B','af_light':'#FAEEE4',
     'spatial':'#7C739C','spatial_light':'#EEEAF4','cell':'#5C776C',
     'cell_light':'#ECF1EE','grey':'#F6F8F9','white':'#FFFFFF'}
WIDTH, HEIGHT = 180, 220
TEXTS=[]

def text(ax,x,y,s,fs=8,color=None,weight='normal',ha='left',va='center',**kw):
    t=ax.text(x,y,s,fontsize=fs,color=color or P['ink'],fontweight=weight,
              ha=ha,va=va,linespacing=1.2,**kw)
    TEXTS.append(t);return t

def panel(ax,x,y,w,h,fill='white',edge=None,r=2,lw=.75):
    p=FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',
                    facecolor=fill,edgecolor=edge or P['edge'],linewidth=lw)
    ax.add_patch(p);return p

def arrow(ax,x,y,xx,yy,dashed=False,color=None,lw=.85,head=True,rad=0):
    ax.add_patch(FancyArrowPatch((x,y),(xx,yy),arrowstyle='-|>' if head else '-',
                                connectionstyle=f'arc3,rad={rad}',
                                mutation_scale=7,color=color or P['muted'],
                                linewidth=lw,linestyle=(0,(3,2)) if dashed else '-'))

def cell(ax,x,y,r=3.5,color=None,angle=0):
    color=color or P['np']
    verts=[]
    for i in range(37):
        a=i*2*math.pi/36;rr=r*(1+.08*math.cos(5*a+angle))
        verts.append((x+rr*math.cos(a),y+.76*rr*math.sin(a)))
    ax.add_patch(PathPatch(MplPath(verts,[MplPath.MOVETO]+[MplPath.LINETO]*35+[MplPath.CLOSEPOLY]),
                           facecolor=matplotlib.colors.to_rgba(color,.12),edgecolor=color,lw=.85))
    ax.add_patch(Ellipse((x+.3,y-.15),r*.8,r*.65,facecolor=matplotlib.colors.to_rgba(color,.45),edgecolor=color,lw=.5))

def donor(ax,x,y,color=None,scale=1):
    color=color or P['np']
    ax.add_patch(Circle((x,y-2*scale),.85*scale,facecolor=color,edgecolor='none'))
    ax.add_patch(FancyBboxPatch((x-1*scale,y-.75*scale),2*scale,2.7*scale,
                               boxstyle=f'round,pad=0,rounding_size={.5*scale}',facecolor=color,edgecolor='none'))
    ax.plot([x-.55*scale,x-.7*scale],[y+1.7*scale,y+3*scale],color=color,lw=1.3*scale)
    ax.plot([x+.55*scale,x+.7*scale],[y+1.7*scale,y+3*scale],color=color,lw=1.3*scale)

def disc(ax,x,y,rx=7,ry=4,with_labels=False):
    for k,alpha in [(1,.12),(.8,.08),(.58,.04)]:
        ax.add_patch(Ellipse((x,y),rx*2*k,ry*2*k,facecolor=matplotlib.colors.to_rgba(P['af'],alpha),
                             edgecolor=P['af'],lw=.7))
    ax.add_patch(Ellipse((x,y),rx*.85,ry*1.05,facecolor=P['np_light'],edgecolor=P['np'],lw=.8))
    if with_labels:
        text(ax,x,y,'NP',fs=8,ha='center',color=P['np'],weight='bold')

def matrix(ax,x,y,cols=5,rows=4,step=2.5):
    # A fixed, non-data icon. Uniform ranks cannot be mistaken for inferred values.
    shades=['#E5EEEE','#D0E1DF','#A6C8C4','#6DA09B']
    for r in range(rows):
        for c in range(cols):
            ax.add_patch(Rectangle((x+c*step,y+r*step),step-.35,step-.35,
                                   facecolor=shades[(r+2*c)%4],edgecolor='white',lw=.25))

def spatial_map(ax,x,y):
    for i in range(9):
        for j in range(6):
            xx=(i-4)*2.1;yy=(j-2.5)*2.1
            if (xx/10)**2+(yy/6.5)**2 < 1:
                inner=(xx/4.5)**2+(yy/3.1)**2 < 1
                ax.add_patch(Circle((x+xx,y+yy),.62,facecolor=P['np'] if inner else P['af'],
                                    alpha=.75,edgecolor='none'))
    ax.add_patch(Ellipse((x,y),21,14,fill=False,edgecolor=P['spatial'],lw=.6))
    ax.add_patch(Ellipse((x,y),9.2,6.6,fill=False,edgecolor=P['np'],lw=.8,linestyle=(0,(2,2))))

def plate(ax,x,y):
    ax.add_patch(Ellipse((x,y),14,5.4,facecolor=P['cell_light'],edgecolor=P['cell'],lw=.8))
    ax.add_patch(Arc((x,y+1),14,5.4,theta1=0,theta2=180,color=P['cell'],lw=.8))
    for dx,dy in [(-4,0),(-1,1),(2,-.5),(4,.8)]:
        cell(ax,x+dx,y+dy,r=1.2,color=P['cell'])
    ax.plot([x-5.6,x+5.6],[y+2.5,y+2.5],color=P['af'],lw=1.5)

def card(ax,x,y,w,h,title,accession,lines,fill,icon):
    panel(ax,x,y,w,h,fill=fill,edge='#C7D3D7',r=1.4)
    icon(ax,x+4.5,y+8)
    text(ax,x+10,y+6.5,title,fs=8.2,weight='bold')
    text(ax,x+10,y+10.7,accession,fs=8)
    text(ax,x+2.5,y+17,'\n'.join(lines),fs=8,va='top')

def check_qa(fig,out):
    fig.canvas.draw();ren=fig.canvas.get_renderer()
    allrec=[]
    for i,t in enumerate(TEXTS):
        bb=t.get_window_extent(ren);data=fig.transFigure.inverted().transform([[bb.x0,bb.y0],[bb.x1,bb.y1]])
        allrec.append({'text':t.get_text(),'fontsize_pt':t.get_fontsize(),
                       'bbox_mm':(data*np.array([WIDTH,HEIGHT])).round(3).tolist(),
                       'within_canvas':bool((data>=0).all() and (data<=1).all())})
    # Report intersections; the reviewer distinguishes intended label/diagram proximities.
    intersections=[]
    for i,a in enumerate(TEXTS):
        ba=a.get_window_extent(ren)
        for j in range(i+1,len(TEXTS)):
            b=TEXTS[j];bb=b.get_window_extent(ren)
            if ba.overlaps(bb):
                intersections.append([a.get_text(),b.get_text()])
    receipt={'dimensions_mm':[WIDTH,HEIGHT],'minimum_font_pt':min(t.get_fontsize() for t in TEXTS),
             'font_requested':'DejaVu Sans','font_resolved':font_manager.findfont('DejaVu Sans'),
             'editable_svg_text':True,'schematic_only':True,'quantitative_source_data':'not applicable',
             'all_text_inside_canvas':all(x['within_canvas'] for x in allrec),
             'text_text_intersections':intersections,'texts':allrec,
             'manual_visual_review':'PENDING until separately documented from actual rendered images'}
    (out/'Figure1_AUTOMATED_LAYOUT_QA.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    return receipt

def draw(facts,out,dpi):
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
                         'font.size':8,'svg.fonttype':'none','pdf.fonttype':42,
                         'savefig.facecolor':'white','figure.facecolor':'white'})
    fig=plt.figure(figsize=(WIDTH/25.4,HEIGHT/25.4))
    ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,WIDTH);ax.set_ylim(HEIGHT,0);ax.axis('off')

    # A: compact, role-specific cohort cards and the explicitly human cell model.
    text(ax,5,5.5,'A',fs=10.3,weight='bold');text(ax,11,5.5,'Data sources',fs=9.4,weight='bold')
    card(ax,5,11,41,38,'Human NP support','GSE165722',
         [f"{facts['independent_human_np_donors']} adult NP donors",'Independent expression','background'],P['np_light'],lambda a,x,y:cell(a,x,y,r=3))
    card(ax,48,11,41,38,'Human discovery','GSE230808',
         [f"{facts['discovery_np_donors']} adult NP donors + AF",f"{facts['communication_eligible_pairs']} eligible AF–NP pairs",'for communication'],P['af_light'],lambda a,x,y:donor(a,x,y,P['af'],.8))
    card(ax,91,11,41,38,'Adult mouse disc','GSE287499',
         ['Spatial transcriptomics','Anatomical expression','context'],P['spatial_light'],lambda a,x,y:disc(a,x,y,3.4,2.3))
    card(ax,134,11,41,38,'Human NP cells','CP-H097Y',
         ['Immortalized model','FN coating / ITGB1 siRNA','qPCR / CCK-8'],P['cell_light'],lambda a,x,y:cell(a,x,y,r=3,color=P['cell']))
    arrow(ax,26,49.5,26,56.5)
    arrow(ax,51,49.5,51,56.5)
    arrow(ax,86,49.5,92,56.5)
    arrow(ax,112,49.5,150,56.5)
    # The cell-model card is an input to E; it is not fed into the spatial branch.
    ax.plot([175.5,178,178],[38,38,192],color=P['cell'],lw=.7)
    arrow(ax,178,192,175,192,color=P['cell'])

    # B: frozen expression programs; separate score names and scales are explicit.
    panel(ax,5,57,54,87,fill='white',edge='#AFCFC9')
    panel(ax,5,57,54,10,fill=P['np_light'],edge=P['np_light'])
    text(ax,8,62,'B',fs=10.3,weight='bold');text(ax,14,62,'Frozen NP programs',fs=9.1,weight='bold',color=P['np'])
    text(ax,9,72,'cNMF K = 5; fixed Top25',fs=8)
    for i in range(5):
        panel(ax,10+i*8.5,77,7,7,fill=P['np'] if i==1 else '#EDF3F2',
              edge=P['np'] if i==1 else '#CBD9D6',r=.8)
        text(ax,13.5+i*8.5,80.5,str(i+1),fs=8,ha='center',weight='bold' if i==1 else 'normal',color='white' if i==1 else P['muted'])
    text(ax,9,90,'GEP2: ECM-remodelling',fs=8.2,weight='bold')
    matrix(ax,10,96,cols=5,rows=4,step=2.5)
    text(ax,10,109,'schematic',fs=8,color=P['muted'])
    text(ax,27,98,'UCell: discovery',fs=8)
    text(ax,27,105,'AUCell:\nindependent NP',fs=8)
    text(ax,9,117,'Distinct scoring scales',fs=8,color=P['muted'])
    text(ax,9,123,'Gene coverage + expression',fs=8)
    for i in range(4):donor(ax,11+i*3.8,134,P['np'],.65)
    text(ax,29,133,'Donor-level\nsummaries',fs=8)

    # C: LR, intracellular SSC and target are distinct icon and semantic layers.
    panel(ax,63,57,58,87,fill='white',edge='#DAC2B1')
    panel(ax,63,57,58,10,fill=P['af_light'],edge=P['af_light'])
    text(ax,66,62,'C',fs=10.3,weight='bold');text(ax,72,62,'Predicted AF → NP',fs=9.1,weight='bold',color=P['af'])
    cell(ax,77,74,r=3.7,color=P['af']);cell(ax,107,74,r=3.7,color=P['np'])
    arrow(ax,82,74,102,74)
    text(ax,77,81,'AF sender',fs=8,ha='center')
    text(ax,106,81,'GEP2-high NP',fs=8,ha='center')
    text(ax,92,88,'CellChat LR + complex keys',fs=8,ha='center')
    # Receptor membrane and ligand are generic shapes, not borrowed illustrations.
    ax.plot([95,116],[95,95],color='#7E989B',lw=1.1)
    ax.plot([95,116],[97,97],color='#7E989B',lw=1.1)
    ax.add_patch(Circle((78,96),1.45,facecolor=P['af'],edgecolor='none'))
    for dx in [-1.2,1.2]:
        ax.plot([104+dx,104+dx],[92.3,97.8],color=P['np'],lw=1.5)
    ax.plot([102.8,105.2],[94,94],color=P['np'],lw=1.5)
    arrow(ax,82,96,99,96,dashed=True,head=False)
    text(ax,92,103,'FN1–integrin context',fs=8,ha='center')
    text(ax,92,111,'Official SigXTalk HGNN',fs=8.3,ha='center',weight='bold')
    text(ax,92,116,'PA / PRS → candidate paths',fs=8,ha='center')
    for xx,role,gene in [(76,'receptor','ITGB1'),(92,'SSC','HIF1A'),(108,'target','PLOD2')]:
        text(ax,xx,123,role,fs=8,ha='center',color=P['muted'])
        if role=='receptor':ax.add_patch(Rectangle((xx-1.4,127-1.4),2.8,2.8,facecolor=P['np'],edgecolor='none'))
        elif role=='SSC':ax.add_patch(Circle((xx,127),1.45,facecolor=P['spatial'],edgecolor='none'))
        else:ax.add_patch(Rectangle((xx-1.4,127-1.4),2.8,2.8,facecolor='#6A8296',edgecolor='none'))
        text(ax,xx,132,gene,fs=8,ha='center')
    arrow(ax,79,127,89,127,dashed=True);arrow(ax,95,127,105,127,dashed=True)
    text(ax,92,139,'Predicted path; SSC activity\nnot measured',fs=8,ha='center',color=P['muted'])

    # D: non-data anatomy, core ROI concept, and fixed-depth sensitivity context.
    panel(ax,125,57,50,87,fill='white',edge='#C9C0D8')
    panel(ax,125,57,50,10,fill=P['spatial_light'],edge=P['spatial_light'])
    text(ax,128,62,'D',fs=10.3,weight='bold');text(ax,134,62,'Adult spatial context',fs=9.1,weight='bold',color=P['spatial'])
    text(ax,150,72,'H&E-readable disc tissue',fs=8,ha='center')
    disc(ax,150,82,rx=10,ry=5.8,with_labels=True)
    text(ax,165,82,'AF',fs=8,ha='center',color=P['af'],weight='bold')
    text(ax,150,93,'Frozen NP / AF core ROIs',fs=8,ha='center')
    spatial_map(ax,150,105)
    text(ax,150,117,'48 μm bins; schematic',fs=8,ha='center')
    text(ax,150,124,'GEP2 / candidate orthologs',fs=8,ha='center')
    text(ax,150,131,'Fixed depth-sensitivity check',fs=8,ha='center')
    text(ax,150,138.5,'Anatomical context;\ndepth-sensitive',fs=8,ha='center',color=P['muted'])

    # E: no composite score and no asserted molecular causal transfer.
    arrow(ax,32,145,50,151);arrow(ax,92,145,92,151);arrow(ax,150,145,135,151)
    panel(ax,5,152,170,28,fill=P['grey'],edge='#C5D1D3')
    text(ax,8,157.8,'E',fs=10.3,weight='bold')
    text(ax,14,157.8,'Integrated interpretation and cell assessment',fs=9.3,weight='bold')
    text(ax,90,164,'Complementary evidence • no composite score',fs=8.2,ha='center')
    panel(ax,10,169,76,7,fill='white',edge='#BFD0CE',r=1)
    panel(ax,90,169,80,7,fill='white',edge='#BFD0CE',r=1)
    text(ax,48,172.5,'FN1–integrin candidate context',fs=8.2,ha='center')
    text(ax,130,172.5,'PLOD2-associated NP transcriptional response',fs=8.2,ha='center')
    arrow(ax,90,181,90,185.5)
    panel(ax,5,186,170,26,fill=P['cell_light'],edge='#C4D3CC')
    text(ax,9,190.5,'Cell-based assessment',fs=8.7,weight='bold',color=P['cell'])
    plate(ax,17,201)
    text(ax,28,200,'FN surface\ncoating',fs=8)
    arrow(ax,51,201,61,201)
    for i,xx in enumerate([68,81,94]):cell(ax,xx,198,r=2.2,color=P['cell'],angle=i)
    text(ax,83,206,'siNC / siITGB1-1 / siITGB1-2',fs=8,ha='center')
    arrow(ax,107,201,119,201)
    text(ax,124,198,'RT-qPCR: PLOD2 primary',fs=8)
    text(ax,124,205,'CCK-8: supportive readout',fs=8)
    text(ax,90,216.3,'Solid arrows: analysis flow; dashed links: candidate associations',fs=8,ha='center',color=P['muted'])

    qa=check_qa(fig,out)
    stem=out/'Figure1_Study_Design_Workflow'
    for ext in ['svg','pdf']:
        fig.savefig(str(stem)+'.'+ext)
    fig.savefig(str(stem)+'.png',dpi=dpi)
    fig.savefig(str(stem)+'.tiff',dpi=dpi,pil_kwargs={'compression':'tiff_lzw'})
    # Exact 180-mm export, without a changing tight bbox.
    fig.savefig(out/'Figure1_180mm_QA_preview.png',dpi=180)
    plt.close(fig)
    manifest={}
    for f in sorted(out.glob('Figure1*')):
        if f.is_file():manifest[f.name]={'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
    (out/'Figure1_EXPORT_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return qa

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--facts',type=Path,required=True)
    ap.add_argument('--output-root','--output',dest='output',type=Path,required=True);ap.add_argument('--dpi',type=int,default=600)
    args=ap.parse_args();facts=json.loads(args.facts.read_text(encoding='utf-8'))
    if args.dpi<300:raise ValueError('Final raster exports must be at least 300 dpi')
    for key,value in [('discovery_np_donors',8),('independent_human_np_donors',8),('communication_eligible_pairs',6)]:
        if facts[key]!=value:raise ValueError('Facts differ from the requested frozen design; review input provenance')
    args.output.mkdir(parents=True,exist_ok=True)
    receipt=draw(facts,args.output,args.dpi)
    print(json.dumps({k:receipt[k] for k in ['dimensions_mm','minimum_font_pt','all_text_inside_canvas','text_text_intersections']},indent=2))
