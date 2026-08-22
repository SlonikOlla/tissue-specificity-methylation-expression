import json
from pathlib import Path

import pandas as pd

from fixed_tsi_analysis import OUT, group_tests, summarize_correlations
from run_extreme_tsi_correlation import (
    human_expr, human_mapping, mouse_expr, mouse_mapping,
    stats_rows, stream_human, stream_mouse,
)


def main():
    human = pd.read_csv(OUT/'human_fixed_genes.csv')
    mouse = pd.read_csv(OUT/'mouse_fixed_genes.csv')

    he, ht = human_expr()
    hm = stream_human(human_mapping(human))
    hr = stats_rows('Human', human, he, hm, ht)
    hr.to_csv(OUT/'human_fixed_correlations.csv', index=False)

    me, mt = mouse_expr()
    mm = stream_mouse(mouse_mapping(mouse))
    mr = stats_rows('Mouse', mouse, me, mm, mt)
    mr.to_csv(OUT/'mouse_fixed_correlations.csv', index=False)

    allr = pd.concat([hr, mr], ignore_index=True)
    summarize_correlations(allr).to_csv(OUT/'fixed_correlation_summary.csv', index=False)
    group_tests(allr).to_csv(OUT/'fixed_group_comparison.csv', index=False)
    meta = {'human_tests': len(hr), 'mouse_tests': len(mr),
            'human_genes_with_data': int(hr.Gene_ID.nunique()),
            'mouse_genes_with_data': int(mr.Gene_ID.nunique())}
    (OUT/'fixed_run_meta.json').write_text(json.dumps(meta, indent=2))
    print(pd.read_csv(OUT/'fixed_correlation_summary.csv').to_string(index=False))
    print(pd.read_csv(OUT/'fixed_group_comparison.csv').to_string(index=False))
    print(meta)


if __name__ == '__main__':
    main()
