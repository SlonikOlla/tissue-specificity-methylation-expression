import json
from pathlib import Path

import pandas as pd

OUT = Path('/tmp/fixed_tsi')


def hgroup(x):
    x = str(x)
    if x.startswith('Brain - '): return 'Brain'
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


def counts(human, mouse, pairs, group, tissue=None):
    h = human[human.TSI_group == group]
    m = mouse[mouse.TSI_group == group]
    if tissue is not None:
        h = h[h.Tissue_group == tissue]
        m = m[m.Tissue_group == tissue]
    hs, ms = set(h.Gene_ID), set(m.Gene_ID)
    hp = {(a,b) for a,b in pairs if a in hs}
    mp = {(a,b) for a,b in pairs if b in ms}
    shared = hp & mp
    return {'human_total':len(hp),'mouse_total':len(mp),'shared':len(shared),
            'human_only':len(hp-shared),'mouse_only':len(mp-shared),
            'union':len(hp|mp),'jaccard':len(shared)/len(hp|mp) if hp|mp else 0}


human = pd.read_csv(OUT/'human_fixed_genes.csv')
mouse = pd.read_csv(OUT/'mouse_fixed_genes.csv')
human['Tissue_group'] = human.Dominant_tissue.map(hgroup)
mouse['Tissue_group'] = mouse.Dominant_tissue.map(mgroup)
ens = pd.read_csv('/tmp/ensembl_human_mouse.tsv',sep='\t',dtype=str).fillna('')
ens = ens[ens['Mouse homology type']=='ortholog_one2one']
pairs = set(zip(ens['Gene stable ID'],ens['Mouse gene stable ID']))
tissues = ['Adrenal gland','Brain','Heart','Kidney','Large intestine','Liver','Lung','Ovary',
           'Skeletal muscle','Small intestine','Spleen','Stomach','Testis','Uterus','Vagina']
data = {'overall':{},'tissues':{}}
for group in ['TSI > 0.8','TSI < 0.4']:
    data['overall'][group] = counts(human,mouse,pairs,group)
for tissue in tissues:
    data['tissues'][tissue] = {g:counts(human,mouse,pairs,g,tissue) for g in ['TSI > 0.8','TSI < 0.4']}
(OUT/'venn_counts.json').write_text(json.dumps(data,indent=2))
rows=[]
for label,groups in [('All tissues combined',data['overall']),*data['tissues'].items()]:
    for group,v in groups.items(): rows.append({'Scope':label,'Group':group,**v})
pd.DataFrame(rows).to_csv(OUT/'venn_counts.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
