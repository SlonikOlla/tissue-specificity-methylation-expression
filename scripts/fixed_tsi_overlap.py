from collections import Counter
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

OUT = Path('/tmp/fixed_tsi')


def hgroup(x):
    x = str(x)
    if x.startswith('Brain - '):
        return 'Brain'
    return {'Adrenal Gland':'Adrenal gland','Colon - Sigmoid':'Large intestine',
            'Colon - Transverse':'Large intestine','Heart - Atrial Appendage':'Heart',
            'Heart - Left Ventricle':'Heart','Kidney - Cortex':'Kidney',
            'Kidney - Medulla':'Kidney','Liver':'Liver','Lung':'Lung',
            'Muscle - Skeletal':'Skeletal muscle','Ovary':'Ovary',
            'Small Intestine - Terminal Ileum':'Small intestine','Spleen':'Spleen',
            'Stomach':'Stomach','Testis':'Testis','Uterus':'Uterus','Vagina':'Vagina'}.get(x)


def mgroup(x):
    return {'mAg':'Adrenal gland','mBr':'Brain','mHe':'Heart','mKi':'Kidney',
            'mLi':'Liver','mLin':'Large intestine','mLu':'Lung','mMu':'Skeletal muscle',
            'mOv':'Ovary','mSin':'Small intestine','mSp':'Spleen','mSt':'Stomach',
            'mTe':'Testis','mUt':'Uterus','mVg':'Vagina'}.get(str(x))


def overall(human, mouse, pairs, group, hcorr=None, mcorr=None):
    hd = human[human.TSI_group == group]
    md = mouse[mouse.TSI_group == group]
    if hcorr is not None:
        hd = hd[hd.Gene_ID.isin(hcorr)]
        md = md[md.Gene_ID.isin(mcorr)]
    hs, ms = set(hd.Gene_ID), set(md.Gene_ID)
    hp = {(h,m) for h,m in pairs if h in hs}
    mp = {(h,m) for h,m in pairs if m in ms}
    shared = hp & mp
    union = hp | mp
    return len(shared), len(union-shared), len(shared)/len(union) if union else None, len(hp), len(mp)


def tissue_summary(human, mouse, pairs, hcorr=None, mcorr=None):
    hd = human[human.TSI_group == 'TSI > 0.8'].copy()
    md = mouse[mouse.TSI_group == 'TSI > 0.8'].copy()
    subset = 'All threshold genes'
    if hcorr is not None:
        hd = hd[hd.Gene_ID.isin(hcorr)]
        md = md[md.Gene_ID.isin(mcorr)]
        subset = 'Genes with correlation data'
    ht, mt = dict(zip(hd.Gene_ID, hd.Tissue_group)), dict(zip(md.Gene_ID, md.Tissue_group))
    tissues = sorted((set(ht.values()) & set(mt.values())) - {None})
    rows = []
    for t in tissues:
        hp = {(h,m) for h,m in pairs if ht.get(h) == t}
        mp = {(h,m) for h,m in pairs if mt.get(m) == t}
        shared, union = hp & mp, hp | mp
        rows.append({'Subset': subset, 'Tissue': t, 'Human_high_orthologs': len(hp),
                     'Mouse_high_orthologs': len(mp), 'Shared_same_tissue': len(shared),
                     'Human_only': len(hp-shared), 'Mouse_only': len(mp-shared),
                     'Jaccard_overlap': len(shared)/len(union) if union else None,
                     'Human_overlap_fraction': len(shared)/len(hp) if hp else None,
                     'Mouse_overlap_fraction': len(shared)/len(mp) if mp else None})
    return pd.DataFrame(rows)


def main():
    human = pd.read_csv(OUT/'human_fixed_genes.csv')
    mouse = pd.read_csv(OUT/'mouse_fixed_genes.csv')
    human['Tissue_group'] = human.Dominant_tissue.map(hgroup)
    mouse['Tissue_group'] = mouse.Dominant_tissue.map(mgroup)
    ens = pd.read_csv('/tmp/ensembl_human_mouse.tsv', sep='\t', dtype=str).fillna('')
    ens = ens[ens['Mouse homology type'] == 'ortholog_one2one']
    pairs = set(zip(ens['Gene stable ID'], ens['Mouse gene stable ID']))
    hr = pd.read_csv(OUT/'human_fixed_correlations.csv')
    mr = pd.read_csv(OUT/'mouse_fixed_correlations.csv')
    hcorr, mcorr = set(hr.Gene_ID), set(mr.Gene_ID)

    tissue = pd.concat([tissue_summary(human, mouse, pairs),
                        tissue_summary(human, mouse, pairs, hcorr, mcorr)], ignore_index=True)
    tissue.to_csv(OUT/'fixed_tissue_overlap.csv', index=False)

    bench = []
    for label, hc, mc in [('All threshold genes', None, None),
                          ('Genes with correlation data', hcorr, mcorr)]:
        hi = overall(human, mouse, pairs, 'TSI > 0.8', hc, mc)
        lo = overall(human, mouse, pairs, 'TSI < 0.4', hc, mc)
        odds, p = fisher_exact([[hi[0], hi[1]], [lo[0], lo[1]]])
        bench.append({'Subset': label, 'High_shared': hi[0], 'High_nonshared': hi[1],
                      'High_Jaccard': hi[2], 'High_human_orthologs': hi[3],
                      'High_mouse_orthologs': hi[4], 'Low_shared': lo[0],
                      'Low_nonshared': lo[1], 'Low_Jaccard': lo[2],
                      'Low_human_orthologs': lo[3], 'Low_mouse_orthologs': lo[4],
                      'Odds_ratio_high_vs_low': odds, 'Fisher_p': p})
    pd.DataFrame(bench).to_csv(OUT/'fixed_ortholog_benchmark.csv', index=False)
    print(pd.DataFrame(bench).to_string(index=False))
    print(tissue.to_string(index=False))
    print('Unmapped high human:', Counter(human.loc[human.TSI_group.eq('TSI > 0.8') & human.Tissue_group.isna(), 'Dominant_tissue']).most_common(10))


if __name__ == '__main__':
    main()
