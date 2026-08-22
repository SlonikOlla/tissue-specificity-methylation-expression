import gzip, json
import pandas as pd

def info(path):
    d=pd.read_csv(path,sep='\t',compression='gzip',dtype=str,usecols=['#tax_id','GeneID','Symbol'])
    return dict(zip(d.GeneID,d.Symbol))

hinfo=info('/tmp/human_gene_info.gz'); minfo=info('/tmp/mouse_gene_info.gz')
o=pd.read_csv('/tmp/gene_orthologs.gz',sep='\t',compression='gzip',dtype=str)
o=o[((o['#tax_id']=='9606')&(o['Other_tax_id']=='10090'))|((o['#tax_id']=='10090')&(o['Other_tax_id']=='9606'))]
pairs=set()
for r in o.itertuples(index=False):
    if r._0=='9606': h,m=r.GeneID,r.Other_GeneID
    else: h,m=r.Other_GeneID,r.GeneID
    if h in hinfo and m in minfo: pairs.add((hinfo[h],minfo[m]))
# Strict one-to-one NCBI orthologs.
hc={}; mc={}
for h,m in pairs: hc.setdefault(h,set()).add(m); mc.setdefault(m,set()).add(h)
pairs={(h,m) for h,m in pairs if len(hc[h])==1 and len(mc[m])==1}

hd=pd.read_csv('/tmp/tsi_corr/human_deciles.csv'); md=pd.read_csv('/tmp/tsi_corr/mouse_deciles.csv')
hr=pd.read_csv('/tmp/tsi_corr/human_correlations.csv'); mr=pd.read_csv('/tmp/tsi_corr/mouse_correlations.csv')

def sets(group, correlated=False):
    if correlated:
        hs=set(hr.loc[hr.TSI_group==group,'Gene_symbol'].dropna())
        ms=set(mr.loc[mr.TSI_group==group,'Gene_symbol'].dropna())
    else:
        hs=set(hd.loc[hd.TSI_group==group,'Gene_symbol'].dropna())
        ms=set(md.loc[md.TSI_group==group,'Gene_symbol'].dropna())
    hp={(h,m) for h,m in pairs if h in hs}; mp={(h,m) for h,m in pairs if m in ms}
    common=hp&mp
    return {'human_only':len(hp-common),'overlap':len(common),'mouse_only':len(mp-common),
            'human_total':len(hp),'mouse_total':len(mp),
            'overlap_genes':[{'human':h,'mouse':m} for h,m in sorted(common)]}

out={
 'top_full':sets('Top 10%',False),'bottom_full':sets('Bottom 10%',False),
 'top_correlated':sets('Top 10%',True),'bottom_correlated':sets('Bottom 10%',True),
 'ortholog_pairs':len(pairs)
}
json.dump(out,open('/tmp/tsi_corr/venn.json','w'),indent=2)
for k,v in out.items():
    if isinstance(v,dict): print(k,{x:y for x,y in v.items() if x!='overlap_genes'})
