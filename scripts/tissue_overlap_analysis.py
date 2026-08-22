import gzip, json, re, zipfile
from collections import Counter, defaultdict
from xml.etree.ElementTree import iterparse
import pandas as pd
from scipy.stats import fisher_exact

XLSX='outputs/tsi_public_data/Human_Mouse_TSI_Methylation_Correlation.xlsx'
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def read_sheet(sheet_no):
    rows=[]
    with zipfile.ZipFile(XLSX) as z, z.open(f'xl/worksheets/sheet{sheet_no}.xml') as f:
        for ev,el in iterparse(f,events=('end',)):
            if el.tag!=NS+'row': continue
            vals={}
            for c in el.findall(NS+'c'):
                ref=c.attrib.get('r','A1'); col=re.match(r'[A-Z]+',ref).group()
                if c.attrib.get('t')=='inlineStr':
                    node=c.find('.//'+NS+'t'); v=node.text if node is not None else ''
                else:
                    node=c.find(NS+'v'); v=node.text if node is not None else ''
                vals[col]=v
            rows.append(vals); el.clear()
    cols=sorted({c for r in rows for c in r},key=lambda x:(len(x),x))
    return pd.DataFrame([[r.get(c,'') for c in cols] for r in rows],columns=cols)

def table(sheet_no,header_row=None):
    x=read_sheet(sheet_no)
    if header_row is None:
        header_row=next(i for i,row in x.iterrows() if 'Species' in set(row.astype(str)))
    headers=x.iloc[header_row].tolist()
    keep=[i for i,h in enumerate(headers) if h!='']
    d=x.iloc[header_row+1:,keep].copy(); d.columns=[headers[i] for i in keep]
    return d

extreme=table(10)
for c in ['TSI_tau']:
    extreme[c]=pd.to_numeric(extreme[c],errors='coerce')
human=extreme[extreme.Species=='Human'].copy(); mouse=extreme[extreme.Species=='Mouse'].copy()
hcorr=set(table(8).Gene_symbol); mcorr=set(table(9).Gene_symbol)

def info(path):
    d=pd.read_csv(path,sep='\t',compression='gzip',dtype=str,usecols=['GeneID','Symbol'])
    return dict(zip(d.GeneID,d.Symbol))
# Ensembl Compara one-to-one orthologs, indexed by the stable IDs used in both expression atlases.
ens=pd.read_csv('/tmp/ensembl_human_mouse.tsv',sep='\t',dtype=str).fillna('')
ens=ens[ens['Mouse homology type']=='ortholog_one2one']
pairs=set(zip(ens['Gene stable ID'],ens['Mouse gene stable ID']))

# Homologous tissue groupings. Human GTEx subregions are aggregated to the corresponding mouse organ.
def hgroup(x):
    if x.startswith('Brain - '): return 'Brain'
    mp={'Adrenal Gland':'Adrenal gland','Colon - Sigmoid':'Large intestine','Colon - Transverse':'Large intestine',
        'Heart - Atrial Appendage':'Heart','Heart - Left Ventricle':'Heart','Kidney - Cortex':'Kidney',
        'Kidney - Medulla':'Kidney','Liver':'Liver','Lung':'Lung','Muscle - Skeletal':'Skeletal muscle',
        'Ovary':'Ovary','Small Intestine - Terminal Ileum':'Small intestine','Spleen':'Spleen','Stomach':'Stomach',
        'Testis':'Testis','Uterus':'Uterus','Vagina':'Vagina'}
    return mp.get(x)
def mgroup(x):
    mp={'mAg':'Adrenal gland','mBr':'Brain','mHe':'Heart','mKi':'Kidney','mLi':'Liver','mLin':'Large intestine',
        'mLu':'Lung','mMu':'Skeletal muscle','mOv':'Ovary','mSin':'Small intestine','mSp':'Spleen','mSt':'Stomach',
        'mTe':'Testis','mUt':'Uterus','mVg':'Vagina'}
    return mp.get(x)
human['Tissue_group']=human.Dominant_tissue.map(hgroup); mouse['Tissue_group']=mouse.Dominant_tissue.map(mgroup)

def analyze(correlation_subset=False):
    hd=human[human.TSI_group=='Top 10%']; md=mouse[mouse.TSI_group=='Top 10%']
    if correlation_subset:
        hd=hd[hd.Gene_symbol.isin(hcorr)]; md=md[md.Gene_symbol.isin(mcorr)]
    ht=dict(zip(hd.Gene_ID,hd.Tissue_group)); mt=dict(zip(md.Gene_ID,md.Tissue_group))
    hsym=dict(zip(hd.Gene_ID,hd.Gene_symbol)); msym=dict(zip(md.Gene_ID,md.Gene_symbol))
    tissues=sorted(set(x for x in ht.values() if x)&set(x for x in mt.values() if x))
    rows=[]; detail=[]
    for t in tissues:
        hp={(h,m) for h,m in pairs if ht.get(h)==t}
        mp={(h,m) for h,m in pairs if mt.get(m)==t}
        shared=hp&mp; union=hp|mp
        rows.append({'Tissue':t,'Human_top_orthologs':len(hp),'Mouse_top_orthologs':len(mp),
                     'Shared_same_tissue':len(shared),'Human_only':len(hp-shared),'Mouse_only':len(mp-shared),
                     'Jaccard_overlap':len(shared)/len(union) if union else None,
                     'Human_overlap_fraction':len(shared)/len(hp) if hp else None,
                     'Mouse_overlap_fraction':len(shared)/len(mp) if mp else None})
        for h,m in sorted(shared):detail.append({'Tissue':t,'Human_gene_ID':h,'Human_symbol':hsym.get(h,''),'Mouse_gene_ID':m,'Mouse_symbol':msym.get(m,''),'Subset':'Correlation data' if correlation_subset else 'All decile genes'})
    return pd.DataFrame(rows),pd.DataFrame(detail)

full,details1=analyze(False); corr,details2=analyze(True)
full.insert(0,'Subset','All decile genes');corr.insert(0,'Subset','Genes with correlation data')
summary=pd.concat([full,corr],ignore_index=True)
details=pd.concat([details1,details2],ignore_index=True)

def transition_table(correlation_subset=False):
    hd=human[human.TSI_group=='Top 10%'];md=mouse[mouse.TSI_group=='Top 10%']
    if correlation_subset:
        hd=hd[hd.Gene_symbol.isin(hcorr)];md=md[md.Gene_symbol.isin(mcorr)]
    hg=hd.set_index('Gene_ID');mg=md.set_index('Gene_ID')
    shared=[(h,m) for h,m in pairs if h in hg.index and m in mg.index]
    rows=[]
    for h,m in shared:
        rows.append({'Subset':'Genes with correlation data' if correlation_subset else 'All decile genes',
          'Human_gene_ID':h,'Human_symbol':hg.loc[h,'Gene_symbol'],'Human_dominant_tissue':hg.loc[h,'Dominant_tissue'],
          'Human_tissue_group':hg.loc[h,'Tissue_group'] or 'Unmatched human tissue',
          'Mouse_gene_ID':m,'Mouse_symbol':mg.loc[m,'Gene_symbol'],'Mouse_dominant_tissue':mg.loc[m,'Dominant_tissue'],
          'Mouse_tissue_group':mg.loc[m,'Tissue_group'] or 'Unmatched mouse tissue',
          'Same_homologous_tissue':bool(hg.loc[h,'Tissue_group'] and hg.loc[h,'Tissue_group']==mg.loc[m,'Tissue_group'])})
    return pd.DataFrame(rows)
transitions=pd.concat([transition_table(False),transition_table(True)],ignore_index=True)
transition_summary=(transitions.groupby(['Subset','Human_tissue_group','Mouse_tissue_group'],dropna=False).size().reset_index(name='Ortholog_pairs'))

# Overall benchmark from prior one-to-one analysis, recalculated here.
def overall(group,subset=False):
    hd=human[human.TSI_group==group];md=mouse[mouse.TSI_group==group]
    if subset:hd=hd[hd.Gene_symbol.isin(hcorr)];md=md[md.Gene_symbol.isin(mcorr)]
    hs=set(hd.Gene_ID);ms=set(md.Gene_ID);hp={(h,m) for h,m in pairs if h in hs};mp={(h,m) for h,m in pairs if m in ms};c=hp&mp
    return len(c),len((hp|mp)-c),len(c)/len(hp|mp)
bench=[]
for subset,label in [(False,'All decile genes'),(True,'Genes with correlation data')]:
    top=overall('Top 10%',subset);bot=overall('Bottom 10%',subset);od,p=fisher_exact([[top[0],top[1]],[bot[0],bot[1]]])
    bench.append({'Subset':label,'Top_shared':top[0],'Top_nonshared':top[1],'Top_Jaccard':top[2],
                  'Bottom_shared':bot[0],'Bottom_nonshared':bot[1],'Bottom_Jaccard':bot[2],'Odds_ratio':od,'Fisher_p':p})
benchmark=pd.DataFrame(bench)

summary.to_csv('/tmp/tissue_overlap_summary.csv',index=False)
details.to_csv('/tmp/tissue_overlap_genes.csv',index=False)
benchmark.to_csv('/tmp/tissue_overlap_benchmark.csv',index=False)
transitions.to_csv('/tmp/tissue_overlap_transitions.csv',index=False)
transition_summary.to_csv('/tmp/tissue_overlap_transition_summary.csv',index=False)
print('Unmapped human dominant:',Counter(human.loc[human.TSI_group=='Top 10%'].loc[human.Tissue_group.isna(),'Dominant_tissue']).most_common())
print(summary.to_string(index=False));print(benchmark.to_string(index=False));print('shared detail',len(details));print(transition_summary.to_string(index=False))
