import json, math, re, zipfile
from pathlib import Path
from xml.etree.ElementTree import iterparse

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

XLSX = Path('outputs/tsi_public_data/Human_Mouse_TSI_Tissue_Overlap.xlsx')
OUT = Path('/tmp/fixed_tsi')
OUT.mkdir(exist_ok=True)
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


def read_sheet(sheet_no):
    rows = []
    with zipfile.ZipFile(XLSX) as z, z.open(f'xl/worksheets/sheet{sheet_no}.xml') as f:
        for _, el in iterparse(f, events=('end',)):
            if el.tag != NS + 'row':
                continue
            vals = {}
            for c in el.findall(NS + 'c'):
                col = re.match(r'[A-Z]+', c.attrib.get('r', 'A1')).group()
                if c.attrib.get('t') == 'inlineStr':
                    node = c.find('.//' + NS + 't')
                else:
                    node = c.find(NS + 'v')
                vals[col] = node.text if node is not None else ''
            rows.append(vals)
            el.clear()
    cols = sorted({c for row in rows for c in row}, key=lambda x: (len(x), x))
    return pd.DataFrame([[row.get(c, '') for c in cols] for row in rows], columns=cols)


def table(sheet_no, header_token):
    x = read_sheet(sheet_no)
    header_row = next(i for i, row in x.iterrows() if header_token in set(row.astype(str)))
    headers = x.iloc[header_row].tolist()
    keep = [i for i, h in enumerate(headers) if h != '']
    d = x.iloc[header_row + 1:, keep].copy()
    d.columns = [headers[i] for i in keep]
    return d


def fixed_groups(d, species):
    d = d.copy()
    d['TSI_tau'] = pd.to_numeric(d['TSI_tau'], errors='coerce')
    d = d[d.TSI_tau.notna()]
    low = d[d.TSI_tau < 0.4].copy()
    high = d[d.TSI_tau > 0.8].copy()
    low['TSI_group'] = 'TSI < 0.4'
    high['TSI_group'] = 'TSI > 0.8'
    out = pd.concat([low, high], ignore_index=True)
    out.insert(0, 'Species', species)
    return out, {'Species': species, 'All_eligible_genes': len(d),
                 'TSI_below_0_4': len(low), 'Percent_below_0_4': len(low)/len(d),
                 'TSI_above_0_8': len(high), 'Percent_above_0_8': len(high)/len(d)}


def summarize_correlations(d):
    rows = []
    for (sp, group, region), x in d.groupby(['Species', 'TSI_group', 'Region']):
        rows.append({'Species': sp, 'TSI_group': group, 'Region': region,
                     'Genes_with_data': x.Gene_ID.nunique(), 'Tests': len(x),
                     'Median_Pearson_r': x.Pearson_r.median(),
                     'Mean_Pearson_r': x.Pearson_r.mean(),
                     'Fraction_negative': (x.Pearson_r < 0).mean(),
                     'Fraction_FDR_0_05': (x.Pearson_FDR < 0.05).mean(),
                     'Median_abs_r': x.Pearson_r.abs().median()})
    return pd.DataFrame(rows)


def group_tests(d):
    rows = []
    for (sp, region), x in d.groupby(['Species', 'Region']):
        low = x[x.TSI_group == 'TSI < 0.4'].Pearson_r.dropna()
        high = x[x.TSI_group == 'TSI > 0.8'].Pearson_r.dropna()
        if len(low) and len(high):
            u, p = mannwhitneyu(low.abs(), high.abs(), alternative='two-sided')
            us, ps = mannwhitneyu(low, high, alternative='two-sided')
        else:
            u = p = us = ps = np.nan
        rows.append({'Species': sp, 'Region': region, 'Low_n': len(low), 'High_n': len(high),
                     'Low_median_abs_r': low.abs().median(), 'High_median_abs_r': high.abs().median(),
                     'Abs_r_difference_high_minus_low': high.abs().median()-low.abs().median(),
                     'Mann_Whitney_abs_U': u, 'Mann_Whitney_abs_p': p,
                     'Low_median_signed_r': low.median(), 'High_median_signed_r': high.median(),
                     'Mann_Whitney_signed_U': us, 'Mann_Whitney_signed_p': ps})
    return pd.DataFrame(rows)


if __name__ == '__main__':
    human_all = table(3, 'Gene_ID')
    mouse_all = table(4, 'Gene_ID')
    human, hs = fixed_groups(human_all, 'Human')
    mouse, ms = fixed_groups(mouse_all, 'Mouse')
    human.to_csv(OUT/'human_fixed_genes.csv', index=False)
    mouse.to_csv(OUT/'mouse_fixed_genes.csv', index=False)
    human_all.to_csv(OUT/'human_all_genes.csv', index=False)
    mouse_all.to_csv(OUT/'mouse_all_genes.csv', index=False)
    pd.DataFrame([hs, ms]).to_csv(OUT/'group_sizes.csv', index=False)
    print(pd.DataFrame([hs, ms]).to_string(index=False))
