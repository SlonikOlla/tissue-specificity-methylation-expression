import csv, gzip, json, math, re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

OUT = Path('/tmp/tsi_corr')
OUT.mkdir(exist_ok=True)

def deciles(path):
    d = pd.read_csv(path).sort_values(['TSI_tau','Gene_ID'], kind='mergesort')
    n = int(math.ceil(len(d) * 0.10))
    low = d.iloc[:n].copy(); high = d.iloc[-n:].copy()
    low['TSI_group']='Bottom 10%'; high['TSI_group']='Top 10%'
    return pd.concat([low, high], ignore_index=True), n

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p)
    x=p[ok]; order=np.argsort(x); ranked=x[order]*len(x)/np.arange(1,len(x)+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]; tmp=np.empty(len(x)); tmp[order]=np.minimum(ranked,1)
    out[np.where(ok)[0]]=tmp; return out

def stats_rows(species, genes, expr, meth, tissues):
    meta=genes.set_index('Gene_ID' if species=='Mouse' else 'Gene_symbol')
    rows=[]
    for (gene,region), mv in meth.items():
        if gene not in meta.index or gene not in expr: continue
        ev=expr[gene]
        xs=[]; ys=[]; used=[]
        for t in tissues:
            if t in mv and t in ev and np.isfinite(mv[t]) and np.isfinite(ev[t]):
                xs.append(mv[t]); ys.append(np.log2(ev[t]+1)); used.append(t)
        if len(xs)<5 or np.std(xs)==0 or np.std(ys)==0: continue
        pr,pp=pearsonr(xs,ys); sr,sp=spearmanr(xs,ys)
        m=meta.loc[gene]
        if isinstance(m,pd.DataFrame): m=m.iloc[0]
        rows.append(dict(Species=species,Gene_ID=m.get('Gene_ID',gene),Gene_symbol=m.get('Gene_symbol',gene),
                         TSI_group=m.TSI_group,TSI_tau=m.TSI_tau,Region=region,N_tissues=len(xs),
                         Tissues='; '.join(used),Pearson_r=pr,Pearson_p=pp,Spearman_rho=sr,Spearman_p=sp,
                         Mean_methylation=float(np.mean(xs)),Mean_log2_expression=float(np.mean(ys))))
    d=pd.DataFrame(rows)
    if len(d):
        d['Pearson_FDR']=bh(d.Pearson_p); d['Spearman_FDR']=bh(d.Spearman_p)
    return d

def human_meta():
    titles=tissues=None
    with gzip.open('/tmp/human_series.txt.gz','rt',errors='replace') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                titles=next(csv.reader([line],delimiter='\t'))[1:]
                titles=[x.strip('"') for x in titles]
            elif line.startswith('!Sample_characteristics_ch1'):
                vals=next(csv.reader([line],delimiter='\t'))[1:]
                if vals and vals[0].strip('"').startswith('tissue:'):
                    tissues=[x.strip('"').split(': ',1)[1] for x in vals]
    return dict(zip(titles,tissues))

def human_expr():
    with gzip.open('/tmp/gtex_median.gct.gz','rt') as f:
        f.readline(); f.readline(); hdr=f.readline().rstrip().split('\t')
    d=pd.read_csv('/tmp/gtex_median.gct.gz',sep='\t',skiprows=2)
    map_t={'Breast - Mammary Tissue':'Breast - Mammary Tissue','Colon - Transverse':'Colon - Transverse',
           'Kidney - Cortex':'Kidney - Cortex','Lung':'Lung','Muscle - Skeletal':'Muscle - Skeletal',
           'Ovary':'Ovary','Prostate':'Prostate','Testis':'Testis','Whole Blood':'Whole Blood'}
    result={}
    for _,r in d.iterrows():
        result[str(r['Description'])]={t:float(r[col]) for t,col in map_t.items() if col in d.columns}
    return result, list(map_t)

def human_mapping(genes):
    symbols=set(genes.Gene_symbol.dropna().astype(str)); mp=defaultdict(list)
    man=pd.read_csv('/tmp/aclust/EPIC.hg38.manifest.csv',usecols=['Probe_ID','UCSC_RefGene_Name','UCSC_RefGene_Group'])
    for row in man.itertuples(index=False):
        ns=str(row.UCSC_RefGene_Name).split(';'); gs=str(row.UCSC_RefGene_Group).split(';')
        for i,n in enumerate(ns):
            if n not in symbols: continue
            g=gs[i] if i<len(gs) else ''
            region='Promoter' if g in {'TSS1500','TSS200','5UTR','1stExon'} else ('Gene body' if g=='Body' else None)
            if region: mp[row.Probe_ID].append((n,region))
    return mp

def stream_human(mp):
    sm=human_meta()
    with gzip.open('/tmp/human_meth.txt.gz','rt') as f:
        hdr=next(csv.reader(f)); samples=hdr[1:]
        idx=defaultdict(list)
        for i,s in enumerate(samples):
            if s in sm: idx[sm[s]].append(i)
        sums=defaultdict(lambda: defaultdict(float)); counts=defaultdict(lambda: defaultdict(int))
        for line in f:
            probe=line.split(',',1)[0].strip('"')
            if probe not in mp: continue
            vals=np.fromstring(line[line.find(',')+1:],sep=',')
            pm={t:float(np.nanmean(vals[ii])) for t,ii in idx.items()}
            for key in mp[probe]:
                for t,v in pm.items(): sums[key][t]+=v; counts[key][t]+=1
    return {k:{t:s/counts[k][t] for t,s in v.items()} for k,v in sums.items()}

def mouse_expr():
    d=pd.read_excel('/tmp/mouse_bodymap.xls',sheet_name='Supplementary Table S6',header=1)
    code={'mAg':'Adrenal gland','mBr':'Brain','mHe':'Heart','mKi':'Kidney','mLi':'Liver','mLu':'Lung',
          'mSp':'Spleen','mSt':'Stomach','mTe':'Testis','mTh':'Thymus','mUt':'Uterus'}
    result={}
    for _,r in d.iterrows():
        ev={}
        for c,t in code.items():
            cols=[x for x in d.columns if str(x).startswith(c+'_')]
            ev[t]=float(pd.to_numeric(r[cols],errors='coerce').mean())
        result[str(r.gene_id)]=ev
    return result,list(code.values())

def mouse_mapping(genes):
    wanted=set(genes.Gene_ID); rec=[]
    pat=re.compile(r'gene_id "([^"]+)";.*gene_name "([^"]+)"')
    with gzip.open('/tmp/mm10.gtf.gz','rt') as f:
        for line in f:
            if line.startswith('#'): continue
            p=line.rstrip().split('\t')
            if len(p)<9 or p[2]!='gene': continue
            m=pat.search(p[8]);
            if not m or m.group(1) not in wanted: continue
            gid,name=m.groups(); start,end=int(p[3]),int(p[4]); strand=p[6]; tss=start if strand=='+' else end
            prom=(max(1,tss-1500),tss+500) if strand=='+' else (tss-500,tss+1500)
            body=(tss+501,end) if strand=='+' else (start,tss-501)
            rec.append((p[0] if p[0].startswith('chr') else 'chr'+p[0],gid,name,'Promoter',*prom))
            if body[0]<=body[1]: rec.append((p[0] if p[0].startswith('chr') else 'chr'+p[0],gid,name,'Gene body',*body))
    bins=defaultdict(list); B=10000
    for chrom,gid,name,region,a,b in rec:
        for z in range(a//B,b//B+1): bins[(chrom,z)].append((a,b,gid,region))
    man=pd.read_csv('/tmp/aclust/MM285.mm10.manifest.csv',usecols=['seqnames','start','Probe_ID'])
    mp=defaultdict(list)
    for r in man.itertuples(index=False):
        for a,b,gid,region in bins.get((r.seqnames,int(r.start)//B),[]):
            if a<=r.start<=b: mp[r.Probe_ID].append((gid,region))
    return mp

def mouse_meta():
    titles=tissues=types=None
    with gzip.open('/tmp/mouse_series.txt.gz','rt',errors='replace') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                vals=next(csv.reader([line],delimiter='\t'))[1:]
                titles=[re.search(r'\[([^]]+)\]',x).group(1) for x in vals]
            elif line.startswith('!Sample_characteristics_ch1'):
                vals=next(csv.reader([line],delimiter='\t'))[1:]
                if vals and vals[0].strip('"').startswith('tissue:'):
                    tissues=[x.strip('"').split(': ',1)[1] for x in vals]
                elif vals and vals[0].strip('"').startswith('cytosine_conversion_type:'):
                    types=[x.strip('"').split(': ',1)[1] for x in vals]
    return {s:t for s,t,k in zip(titles,tissues,types) if k=='BS'}

def stream_mouse(mp):
    sm=mouse_meta(); sums=defaultdict(lambda: defaultdict(float)); counts=defaultdict(lambda: defaultdict(int))
    tissue_harmonize={'Adrenal':'Adrenal gland','Cerebellum':'Brain','Subcortical_Brain':'Brain'}
    with gzip.open('/tmp/mouse_meth.csv.gz','rt') as f:
        hdr=next(csv.reader(f)); sample_cols=[]; tissues=[]
        for j in range(1,len(hdr),2):
            if hdr[j] in sm:
                sample_cols.append(j-1); tissues.append(tissue_harmonize.get(sm[hdr[j]],sm[hdr[j]]))
        groups=defaultdict(list)
        for pos,t in zip(sample_cols,tissues): groups[t].append(pos)
        for line in f:
            probe=line.split(',',1)[0].strip('"')
            if probe not in mp: continue
            vals=np.fromstring(line[line.find(',')+1:],sep=',')
            pm={t:float(np.nanmean(vals[ii])) for t,ii in groups.items()}
            for key in mp[probe]:
                for t,v in pm.items(): sums[key][t]+=v; counts[key][t]+=1
    return {k:{t:s/counts[k][t] for t,s in v.items()} for k,v in sums.items()}

def summarize(d):
    rows=[]
    for (sp,g,r),x in d.groupby(['Species','TSI_group','Region']):
        rows.append(dict(Species=sp,TSI_group=g,Region=r,Genes_with_data=x.Gene_ID.nunique(),Tests=len(x),
                         Median_Pearson_r=x.Pearson_r.median(),Mean_Pearson_r=x.Pearson_r.mean(),
                         Fraction_negative=(x.Pearson_r<0).mean(),Fraction_FDR_0_05=(x.Pearson_FDR<.05).mean(),
                         Median_abs_r=x.Pearson_r.abs().median()))
    return pd.DataFrame(rows)

if __name__ == '__main__':
    human,n_h=deciles('/tmp/tsi_work/human_genes.csv'); mouse,n_m=deciles('/tmp/tsi_work/mouse_genes.csv')
    human.to_csv(OUT/'human_deciles.csv',index=False); mouse.to_csv(OUT/'mouse_deciles.csv',index=False)
    he,ht=human_expr(); hm=stream_human(human_mapping(human)); hr=stats_rows('Human',human,he,hm,ht)
    hr.to_csv(OUT/'human_correlations.csv',index=False)
    me,mt=mouse_expr(); mm=stream_mouse(mouse_mapping(mouse)); mr=stats_rows('Mouse',mouse,me,mm,mt)
    mr.to_csv(OUT/'mouse_correlations.csv',index=False)
    allr=pd.concat([hr,mr],ignore_index=True); summarize(allr).to_csv(OUT/'summary.csv',index=False)
    with open(OUT/'run_meta.json','w') as f: json.dump({'human_decile_n':n_h,'mouse_decile_n':n_m,'human_tests':len(hr),'mouse_tests':len(mr)},f,indent=2)
    print(pd.read_csv(OUT/'summary.csv').to_string(index=False)); print(json.load(open(OUT/'run_meta.json')))
