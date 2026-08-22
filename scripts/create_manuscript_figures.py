from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch

OUT = Path('/workspace/scratch/f4b527e783e3/manuscript_assets')
OUT.mkdir(exist_ok=True)
BLUE='#2186B5'; ORANGE='#C67820'; NAVY='#173B57'; TEAL='#318D86'; GRAY='#65717C'; LIGHT='#EEF3F6'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10})

def save(fig, name):
    fig.savefig(OUT/name, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

# Figure 1: study design
fig, ax = plt.subplots(figsize=(10,5.2)); ax.set_xlim(0,10); ax.set_ylim(0,5.2); ax.axis('off')
boxes=[(0.3,3.25,2.2,1.15,'Mammalian evidence','Human • mouse • rat\nPublished studies'),
       (2.9,3.25,2.2,1.15,'Expression atlases','GTEx v8 (human)\nBodyMap (mouse)'),
       (5.5,3.25,2.2,1.15,'Methylation atlases','EPIC: GSE213478\nMM285: GSE290585'),
       (2.0,1.15,2.5,1.15,'Gene-level integration','Promoter and gene body\nAcross-tissue correlations'),
       (5.3,1.15,2.5,1.15,'Tissue specificity','Tau; deciles and\nfixed thresholds'),
       (8.0,1.15,1.7,1.15,'Evolution','1:1 orthologues\nTissue overlap')]
for x,y,w,h,title,sub in boxes:
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.03,rounding_size=0.08',fc=LIGHT,ec=NAVY,lw=1.4))
    ax.text(x+w/2,y+h*0.67,title,ha='center',va='center',weight='bold',color=NAVY,fontsize=11)
    ax.text(x+w/2,y+h*0.30,sub,ha='center',va='center',color=GRAY,fontsize=9)
for a,b in [((2.5,3.82),(2.9,3.82)),((5.1,3.82),(5.5,3.82)),((4.0,3.25),(3.3,2.30)),((6.5,3.25),(6.55,2.30)),((4.5,1.72),(5.3,1.72)),((7.8,1.72),(8.0,1.72))]:
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,lw=1.2,color=GRAY))
ax.text(5,4.95,'Study design and analytical progression',ha='center',fontsize=15,weight='bold',color=NAVY)
ax.text(5,0.35,'Rat informed the comparative evidence base; quantitative cross-atlas integration was performed for human and mouse.',ha='center',fontsize=9.5,color=GRAY)
save(fig,'Figure1_study_design.png')

# Figure 2: TSI distributions and fixed group sizes
h=pd.read_csv('/tmp/fixed_tsi/human_all_genes.csv'); m=pd.read_csv('/tmp/fixed_tsi/mouse_all_genes.csv')
for d in (h,m): d['TSI_tau']=pd.to_numeric(d['TSI_tau'],errors='coerce')
fig,axs=plt.subplots(1,2,figsize=(10,4.2),gridspec_kw={'wspace':0.27})
bins=np.linspace(0,1,41)
axs[0].hist(h.TSI_tau.dropna(),bins=bins,density=True,alpha=.55,color=BLUE,label='Human')
axs[0].hist(m.TSI_tau.dropna(),bins=bins,density=True,alpha=.45,color=ORANGE,label='Mouse')
axs[0].axvline(.4,color=GRAY,ls='--'); axs[0].axvline(.8,color=GRAY,ls='--')
axs[0].set(xlabel='Tau tissue-specificity index',ylabel='Density',title='A  TSI distributions'); axs[0].legend(frameon=False)
groups=['TSI < 0.4','TSI > 0.8']; hv=[162,27889]; mv=[523,18718]; x=np.arange(2); w=.34
axs[1].bar(x-w/2,hv,w,color=BLUE,label='Human'); axs[1].bar(x+w/2,mv,w,color=ORANGE,label='Mouse')
axs[1].set_xticks(x,groups); axs[1].set_yscale('log'); axs[1].set_ylabel('Genes (log scale)'); axs[1].set_title('B  Fixed-threshold group sizes'); axs[1].legend(frameon=False)
for i,v in enumerate(hv): axs[1].text(i-w/2,v*1.13,f'{v:,}',ha='center',fontsize=9,color=BLUE)
for i,v in enumerate(mv): axs[1].text(i+w/2,v*1.13,f'{v:,}',ha='center',fontsize=9,color=ORANGE)
for a in axs: a.spines[['top','right']].set_visible(False)
save(fig,'Figure2_tsi_distribution.png')

# Figure 3: fixed threshold signed and absolute correlations
rows=[('Human','Gene body',-.338,-.030,.378,.374),('Human','Promoter',-.286,-.249,.370,.404),
      ('Mouse','Gene body',.095,-.119,.339,.420),('Mouse','Promoter',.010,-.208,.329,.418)]
d=pd.DataFrame(rows,columns=['Species','Region','Low_r','High_r','Low_abs','High_abs'])
fig,axs=plt.subplots(1,2,figsize=(10,4.5),gridspec_kw={'wspace':.28}); labels=[f'{s}\n{r}' for s,r in zip(d.Species,d.Region)]; x=np.arange(4); w=.34
axs[0].bar(x-w/2,d.Low_r,w,color='#A7CFE3',label='TSI < 0.4'); axs[0].bar(x+w/2,d.High_r,w,color=BLUE,label='TSI > 0.8'); axs[0].axhline(0,color='#333',lw=.8)
axs[0].set_xticks(x,labels); axs[0].set_ylabel('Median Pearson r'); axs[0].set_title('A  Correlation direction'); axs[0].legend(frameon=False)
axs[1].bar(x-w/2,d.Low_abs,w,color='#E9C9A4',label='TSI < 0.4'); axs[1].bar(x+w/2,d.High_abs,w,color=ORANGE,label='TSI > 0.8')
axs[1].set_xticks(x,labels); axs[1].set_ylabel('Median |Pearson r|'); axs[1].set_ylim(0,.48); axs[1].set_title('B  Correlation magnitude'); axs[1].legend(frameon=False)
for a in axs: a.spines[['top','right']].set_visible(False)
save(fig,'Figure3_correlations.png')

# Figure 4: threshold-dependent orthologue conservation
fig,axs=plt.subplots(1,3,figsize=(10,3.8),gridspec_kw={'width_ratios':[1,1,1.05],'wspace':.28})
def venn(ax,left,shared,right,title,j):
    ax.set_aspect('equal'); ax.axis('off'); ax.add_patch(Circle((.42,.5),.34,fc='#6DCBF488',ec=BLUE,lw=1.5)); ax.add_patch(Circle((.68,.5),.34,fc='#F4B66D88',ec=ORANGE,lw=1.5))
    ax.text(.22,.52,f'{left:,}',ha='center',weight='bold',fontsize=12); ax.text(.55,.52,f'{shared:,}',ha='center',weight='bold',fontsize=11); ax.text(.86,.52,f'{right:,}',ha='center',weight='bold',fontsize=12)
    ax.text(.22,.82,'Human',ha='center',color=BLUE,weight='bold'); ax.text(.86,.82,'Mouse',ha='center',color=ORANGE,weight='bold'); ax.set_title(f'{title}\nJaccard = {j:.3f}',fontsize=11,weight='bold')
venn(axs[0],1733,6219,1753,'TSI > 0.8',.641); venn(axs[1],116,32,448,'TSI < 0.4',.054)
cats=['Top/bottom\ndecile','Fixed\nthreshold']; top=[.209,.641]; low=[.329,.054]; xx=np.arange(2); ww=.34
axs[2].bar(xx-ww/2,top,ww,color=BLUE,label='High TSI'); axs[2].bar(xx+ww/2,low,ww,color=ORANGE,label='Low TSI'); axs[2].set_xticks(xx,cats); axs[2].set_ylim(0,.72); axs[2].set_ylabel('Jaccard overlap'); axs[2].set_title('Threshold sensitivity'); axs[2].legend(frameon=False,fontsize=8); axs[2].spines[['top','right']].set_visible(False)
save(fig,'Figure4_orthologue_overlap.png')

# Figure 5: tissue overlap
tissues=['Brain','Testis','Liver','Skeletal muscle','Heart','Stomach','Small intestine','Kidney','Lung','Spleen','Large intestine','Adrenal gland','Uterus','Ovary','Vagina']
shared=[923,1097,186,164,44,39,71,77,55,81,32,17,15,14,0]
jac=[.4303,.4147,.4035,.3952,.2366,.2229,.1935,.1930,.1322,.1294,.1000,.0625,.0479,.0342,0]
order=np.arange(len(tissues))[::-1]
fig,ax=plt.subplots(figsize=(8.2,6.2)); sizes=25+np.sqrt(np.array(shared))*5
ax.scatter(np.array(jac)[order],order,s=sizes[order],c=np.array(shared)[order],cmap='viridis',edgecolor='white',lw=.8)
ax.set_yticks(order,[tissues[i] for i in order]); ax.set_xlabel('Jaccard overlap of TSI > 0.8 one-to-one orthologues'); ax.set_title('Tissue-level conservation of high-TSI genes',weight='bold')
ax.grid(axis='x',alpha=.22); ax.spines[['top','right','left']].set_visible(False)
for idx in order: ax.text(jac[idx]+.009,idx,f'n={shared[idx]:,}',va='center',fontsize=8,color=GRAY)
save(fig,'Figure5_tissue_overlap.png')

print(OUT)
