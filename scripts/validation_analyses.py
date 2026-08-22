import json
from pathlib import Path

import numpy as np
import pandas as pd


OUT = Path("/tmp/tsi_validation")
OUT.mkdir(exist_ok=True)
SEED = 20260817
N_ITER = 10_000


def mapped_sets(h, m, pairs, group, exclude_sex=False):
    hs = h[h.TSI_group.eq(group)].copy()
    ms = m[m.TSI_group.eq(group)].copy()
    if exclude_sex:
        hs = hs[~hs.Dominant_tissue.isin({"Testis", "Ovary", "Uterus", "Vagina", "Prostate", "Breast - Mammary Tissue"})]
        ms = ms[~ms.Dominant_tissue.isin({"mTe", "mOv", "mUt", "mVg", "mPr"})]
    hids, mids = set(hs.Gene_ID), set(ms.Gene_ID)
    hp = {(a, b) for a, b in pairs if a in hids}
    mp = {(a, b) for a, b in pairs if b in mids}
    return hp, mp


def metrics(hp, mp):
    shared = len(hp & mp)
    union = len(hp | mp)
    return dict(Human_pairs=len(hp), Mouse_pairs=len(mp), Shared=shared,
                Union=union, Jaccard=shared / union if union else np.nan)


def downsample_compare(high_h, high_m, low_h, low_m, rng):
    """Downsample high sets to low species-specific sizes; return null distribution."""
    nh, nm = len(low_h), len(low_m)
    c = len(high_h & high_m)
    # Number of shared-universe pairs captured by each species-specific draw.
    kh = rng.hypergeometric(c, len(high_h)-c, nh, size=N_ITER)
    km = rng.hypergeometric(c, len(high_m)-c, nm, size=N_ITER)
    shared = np.array([rng.hypergeometric(a, c-a, b) for a,b in zip(kh,km)])
    rows = np.column_stack([shared, shared/(nh+nm-shared)])
    return rows


def random_null(universe, nh, nm, rng):
    # Conditional on the human draw, intersection with the independent mouse draw.
    return rng.hypergeometric(nh, len(universe)-nh, nm, size=N_ITER).astype(float)


def q(x, p): return float(np.quantile(x, p))


def main():
    h = pd.read_csv('/tmp/fixed_tsi/human_fixed_genes.csv')
    m = pd.read_csv('/tmp/fixed_tsi/mouse_fixed_genes.csv')
    ens = pd.read_csv('/tmp/ensembl_human_mouse.tsv', sep='\t', dtype=str).fillna('')
    ens = ens[ens['Mouse homology type'].eq('ortholog_one2one')]
    pairs = set(zip(ens['Gene stable ID'], ens['Mouse gene stable ID']))
    rng = np.random.default_rng(SEED)
    summary, null_rows = [], []

    for exclusion in [False, True]:
        label = 'All organs' if not exclusion else 'Sex-organ-associated dominant genes excluded'
        sets = {}
        for group in ['TSI < 0.4', 'TSI > 0.8']:
            hp, mp = mapped_sets(h, m, pairs, group, exclusion)
            sets[group] = (hp, mp)
            obs = metrics(hp, mp)
            null = random_null(pairs, len(hp), len(mp), rng)
            expected = null.mean()
            summary.append(dict(Analysis=label, Group=group, **obs,
                                Expected_shared_random=float(expected),
                                Fold_enrichment=obs['Shared']/expected if expected else np.nan,
                                Null_shared_CI_low=q(null,.025), Null_shared_CI_high=q(null,.975),
                                Empirical_p_enrichment=(1+np.sum(null >= obs['Shared']))/(N_ITER+1)))
            for i, v in enumerate(null):
                null_rows.append(dict(Analysis=label, Null_type='Random orthologue universe',
                                      Group=group, Iteration=i+1, Shared=int(v), Jaccard=np.nan))

        low_h, low_m = sets['TSI < 0.4']
        high_h, high_m = sets['TSI > 0.8']
        ds = downsample_compare(high_h, high_m, low_h, low_m, rng)
        low_obs = metrics(low_h, low_m)
        summary.append(dict(Analysis=label, Group='High TSI downsampled to low sizes',
                            Human_pairs=len(low_h), Mouse_pairs=len(low_m),
                            Shared=float(ds[:,0].mean()), Union=np.nan,
                            Jaccard=float(ds[:,1].mean()), Expected_shared_random=np.nan,
                            Fold_enrichment=np.nan,
                            Null_shared_CI_low=q(ds[:,0],.025), Null_shared_CI_high=q(ds[:,0],.975),
                            Empirical_p_enrichment=(1+np.sum(ds[:,1] >= low_obs['Jaccard']))/(N_ITER+1)))
        for i, (s,j) in enumerate(ds):
            null_rows.append(dict(Analysis=label, Null_type='High TSI downsampled to low sizes',
                                  Group='TSI > 0.8', Iteration=i+1, Shared=int(s), Jaccard=float(j)))

    pd.DataFrame(summary).to_csv(OUT/'orthology_validation_summary.csv', index=False)
    pd.DataFrame(null_rows).to_csv(OUT/'orthology_resampling_distributions.csv', index=False)

    counts=[]
    for sp,d,sex in [('Human',h,{"Testis","Ovary","Uterus","Vagina","Prostate","Breast - Mammary Tissue"}),
                     ('Mouse',m,{"mTe","mOv","mUt","mVg","mPr"})]:
        for g,x in d.groupby('TSI_group'):
            nsex=x.Dominant_tissue.isin(sex).sum()
            counts.append(dict(Species=sp,Group=g,All_genes=len(x),Sex_organ_associated=nsex,
                               Retained=len(x)-nsex,Retained_fraction=(len(x)-nsex)/len(x)))
    pd.DataFrame(counts).to_csv(OUT/'sex_organ_exclusion_counts.csv',index=False)
    (OUT/'run_metadata.json').write_text(json.dumps({'seed':SEED,'iterations':N_ITER,
        'sex_organs_human':['Testis','Ovary','Uterus','Vagina','Prostate','Breast - Mammary Tissue'],
        'sex_organs_mouse':['mTe','mOv','mUt','mVg','mPr'],
        'note':'Exclusion is based on dominant-expression organ. Existing per-gene correlation coefficients were not recomputed after removing organs.'},indent=2))
    print(pd.DataFrame(summary).to_string(index=False))
    print(pd.DataFrame(counts).to_string(index=False))


if __name__ == '__main__': main()
