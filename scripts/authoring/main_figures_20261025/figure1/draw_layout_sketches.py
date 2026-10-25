"""Two original vector layout studies; geometry only, no scientific measurements."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-root', type=Path, required=True)
OUT = parser.parse_args().output_root.resolve()
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
                     'font.size':8, 'svg.fonttype':'none', 'pdf.fonttype':42})
COLORS = ['#DDEEEA','#FAE8D9','#E8E7F1','#E8EFEB']

def box(ax, x,y,w,h,text,c):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=2',
                              facecolor=c,edgecolor='#65737B',linewidth=.7))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',linespacing=1.45)
def arrow(ax,x,y,xx,yy):
    ax.add_patch(FancyArrowPatch((x,y),(xx,yy),arrowstyle='-|>',mutation_scale=8,
                                color='#566570',linewidth=.8))

for variant in ['01_four_tier','02_side_spine']:
    fig=plt.figure(figsize=(180/25.4,200/25.4))
    ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,180);ax.set_ylim(200,0);ax.axis('off')
    if variant.startswith('01'):
        box(ax,5,5,170,34,'', '#F5F7F8')
        ax.text(90,12,'A  Data sources: three cohorts + cell model',ha='center',va='center')
        for i,t in enumerate(['Human discovery','Human NP support','Adult mouse spatial','Human NP cells']):
            box(ax,9+i*41,19,39,16,t,COLORS[i])
        for x,w,t,c in [(5,52,'B  Frozen programs\n\nFive programs\nFixed gene membership\nSeparate score methods\nDonor summaries',COLORS[0]),
                        (61,58,'C  Predicted AF → NP\n\nLR inference\nOfficial HGNN\nReceptor / SSC / target\nCandidate paths',COLORS[1]),
                        (123,52,'D  Spatial context\n\nDisc schematic\nNP / AF core ROIs\nBins and orthologs\nDepth sensitivity',COLORS[2])]:
            box(ax,x,48,w,88,t,c);arrow(ax,x+w/2,40,x+w/2,47)
            arrow(ax,x+w/2,137,90,145)
        box(ax,5,146,170,17,'E  Complementary evidence → candidate interpretation','#F2F5F4')
        arrow(ax,90,164,90,170)
        box(ax,16,171,148,24,'Cell assessment\nFN coating × two independent ITGB1 siRNAs → qPCR + CCK-8',COLORS[3])
    else:
        box(ax,5,5,170,29,'A  Cohort inputs + cell model','#F5F7F8')
        for i,(t,c) in enumerate(zip(['B  Human programs','C  AF → NP prediction','D  Spatial context'],COLORS)):
            box(ax,5,45+i*42,43,32,t,c)
            box(ax,56,45+i*42,82,32,'Step 1 → schematic → output',c)
            arrow(ax,49,61+i*42,55,61+i*42)
            arrow(ax,139,61+i*42,151,126)
        box(ax,145,45,30,112,'E\nEvidence\nintegration\n\nCandidate\ncontext', '#F2F5F4')
        arrow(ax,160,158,160,169)
        box(ax,5,170,170,25,'Cell assessment: qPCR primary + CCK-8 supportive',COLORS[3])
    fig.savefig(OUT/(variant+'.svg'));fig.savefig(OUT/(variant+'.pdf'))
    fig.savefig(OUT/(variant+'.png'),dpi=180);plt.close(fig)
