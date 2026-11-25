# Table WA5 (columns Papers to Share Z <= -1.96): newer trials, journal articles versus working papers and reports.
# Publication status as of October 2026 (data/paper_status.csv; working papers published since then count as published).
# Also writes data/pubstatus_*.csv for 2_random_effects.R. Run from this folder: python3 1_counts.py
import pandas as pd
st = pd.read_csv('data/paper_status.csv').set_index('study').published
out = []
for dom, f in [('behavior', 'C_full'), ('knowledge', 'C_know')]:
    D = pd.read_csv(f'data/{f}.csv'); D = D[~D.source.str.startswith('Kaiser')].copy(); D['published'] = D.study.map(st).astype(int)
    D[['study', 'cluster', 'yi', 'sei', 'published']].to_csv(f'data/pubstatus_{dom}.csv', index=False)
    for p, x in D.groupby('published'):
        z = x.yi / x.sei
        out.append(dict(domain=dom, source='Journal article or chapter' if p else 'Working paper or report', papers=x.study.nunique(),
                        trials=x.cluster.nunique(), estimates=len(x), share_Z_ge_1_96=(z >= 1.96).mean(), share_Z_lt_0=(z < 0).mean(),
                        share_Z_le_minus_1_96=(z <= -1.96).mean()))
r = pd.DataFrame(out); r.to_csv('output/table_WA5_counts.csv', index=False); print(r.round(3).to_string(index=False))
