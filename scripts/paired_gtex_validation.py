import csv, gzip, re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from run_extreme_tsi_correlation import bh, human_meta

OUT=Path('/tmp/tsi_validation'); OUT.mkdir(exist_ok=True)
SEED=20260817
MIN_N=20

def root_id(x):
    p=str(x).split('-')
    return '-'.join(p[:3]) if len(p)>=3 else str(x)

def choose_genes():
    d=pd.read_csv('/tmp/fixed_tsi/human_fixed_genes.csv')
    d=d[~d.Dominant_tissue.isin({'Testis','Ovary','Uterus','Vagina','Prostate','Breast - Mammary Tissue'})].copy()
    low=d[d.TSI_group.eq('TSI < 0.4')]
    high=d[d.TSI_group.eq('TSI > 0.8')]
    # Candidate high pool spans the expression range while keeping the matrix tractable.
    rng=np.random.default_rng(SEED)
    hi=high.iloc[rng.choice(len(high),min(5000,len(high)),replace=False)]
    return pd.concat([low,hi],ignore_index=True)

def manifest_map(genes):
    symbols=set(genes.Gene_symbol.dropna().astype(str)); mp=defaultdict(list)
    with gzip.open('/tmp/epic_manifest.csv.gz','rt',errors='replace') as f:
        for line in f:
            if line.startswith('IlmnID,Name,'):
                hdr=next(csv.reader([line])); break
        idx={x:i for i,x in enumerate(hdr)}
        for row in csv.reader(f):
            if len(row)<=idx['UCSC_RefGene_Group']: continue
            names=row[idx['UCSC_RefGene_Name']].split(';')
            groups=row[idx['UCSC_RefGene_Group']].split(';')
            for i,n in enumerate(names):
                if n not in symbols: continue
                g=groups[i] if i<len(groups) else ''
                region='Promoter' if g in {'TSS1500','TSS200','5UTR','1stExon'} else ('Gene body' if g=='Body' else None)
                if region: mp[row[idx['Name']]].append((n,region))
    return mp

def expression(genes, meth_samples):
    wanted=set(genes.Gene_ID.astype(str)); expr={}
    with gzip.open('/tmp/gtex_gene_tpm.gct.gz','rt') as f:
        f.readline(); f.readline(); hdr=f.readline().rstrip().split('\t')
        root_to_col={root_id(s):i-2 for i,s in enumerate(hdr[2:],2)}
        matched=[(j,root_to_col.get(root_id(s))) for j,s in enumerate(meth_samples) if root_id(s) in root_to_col]
        use_m=[a for a,b in matched]; use_e=[b for a,b in matched]
        for line in f:
            a,b,vals=line.rstrip().split('\t',2)
            gid=a.split('.')[0]
            if gid not in wanted: continue
            arr=np.fromstring(vals,sep='\t',dtype=np.float32)
            expr[gid]=arr[use_e]
    return expr,use_m

def main():
    genes=choose_genes(); meta=genes.set_index('Gene_symbol')
    sm=human_meta()
    with gzip.open('/tmp/human_meth.txt.gz','rt') as f:
        meth_samples=next(csv.reader(f))[1:]
    expr,use_m=expression(genes,meth_samples)
    samples=[meth_samples[i] for i in use_m]
    tissues=np.array([sm.get(s,'') for s in samples])
    keep_tissues={'Colon - Transverse','Kidney - Cortex','Lung','Muscle - Skeletal','Whole Blood'}
    mp=manifest_map(genes)
    keys=sorted({k for v in mp.values() for k in v})
    key_i={k:i for i,k in enumerate(keys)}
    sums=np.zeros((len(keys),len(samples)),np.float32)
    cnt=np.zeros((len(keys),len(samples)),np.uint16)
    with gzip.open('/tmp/human_meth.txt.gz','rt') as f:
        next(f)
        for line in f:
            probe=line.split(',',1)[0].strip('"')
            targets=mp.get(probe)
            if not targets: continue
            vals=np.fromstring(line[line.find(',')+1:],sep=',',dtype=np.float32)[use_m]
            ok=np.isfinite(vals)
            for k in targets:
                z=key_i[k]; sums[z,ok]+=vals[ok]; cnt[z,ok]+=1
    rows=[]
    for (symbol,region),z in key_i.items():
        mm=np.divide(sums[z],cnt[z],out=np.full(len(samples),np.nan,np.float32),where=cnt[z]>0)
        rec=meta.loc[symbol]
        if isinstance(rec,pd.DataFrame): rec=rec.iloc[0]
        ev=expr.get(str(rec.Gene_ID))
        if ev is None: continue
        for tissue in sorted(keep_tissues):
            ok=(tissues==tissue)&np.isfinite(mm)&np.isfinite(ev)
            if ok.sum()<MIN_N or np.std(mm[ok])==0 or np.std(ev[ok])==0: continue
            r,p=pearsonr(mm[ok],np.log2(ev[ok]+1))
            rows.append(dict(Tissue=tissue,Region=region,Gene_ID=rec.Gene_ID,Gene_symbol=symbol,
                             TSI_group=rec.TSI_group,TSI_tau=rec.TSI_tau,N_pairs=int(ok.sum()),
                             Pearson_r=r,Pearson_p=p))
    d=pd.DataFrame(rows)
    d['FDR']=np.nan
    for _,ix in d.groupby(['Tissue','Region']).groups.items(): d.loc[ix,'FDR']=bh(d.loc[ix,'Pearson_p'])
    d.to_csv(OUT/'paired_gtex_gene_correlations.csv',index=False)
    rng=np.random.default_rng(SEED); ss=[]
    for (t,r),x in d.groupby(['Tissue','Region']):
        lo=x[x.TSI_group.eq('TSI < 0.4')].Pearson_r.to_numpy(); hi=x[x.TSI_group.eq('TSI > 0.8')].Pearson_r.to_numpy()
        if not len(lo) or len(hi)<len(lo): continue
        dif=[]
        for _ in range(5000):
            hs=rng.choice(hi,len(lo),replace=False)
            dif.append(np.median(np.abs(hs))-np.median(np.abs(lo)))
        dif=np.array(dif)
        ss.append(dict(Tissue=t,Region=r,Low_n=len(lo),High_candidate_n=len(hi),
                       Low_median_r=np.median(lo),High_median_r=np.median(hi),
                       Low_median_abs_r=np.median(np.abs(lo)),High_median_abs_r=np.median(np.abs(hi)),
                       Size_matched_abs_difference=np.median(dif),CI_low=np.quantile(dif,.025),CI_high=np.quantile(dif,.975),
                       Probability_high_gt_low=np.mean(dif>0)))
    pd.DataFrame(ss).to_csv(OUT/'paired_gtex_summary.csv',index=False)
    pd.DataFrame({'Sample':samples,'Biospecimen_root':[root_id(x) for x in samples],'Tissue':tissues}).to_csv(OUT/'paired_gtex_samples.csv',index=False)
    print('matched samples',len(samples),'genes',len(d),'tests',len(rows))
    print(pd.DataFrame(ss).to_string(index=False))

if __name__=='__main__': main()
