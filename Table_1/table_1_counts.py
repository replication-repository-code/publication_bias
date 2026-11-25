# Table 1, row "Estimates (trials)" of each panel: the number of estimates and, in parentheses, of distinct trials in each region
# of Z = g / SE (combined behavior data), and in the reference region (note to table 1).
# Run from this folder: python3 table_1_counts.py   (writes output/table_1_counts.csv)
import pandas as pd, numpy as np
b = pd.read_csv('data/C_full.csv'); inf = np.inf
def cnt(d, cuts):
    z = d.yi / d.sei
    return [f"{int(((z >= a) & (z < c)).sum())} ({d[(z >= a) & (z < c)].cluster.nunique()})" for a, c in cuts]
names = ['Z < -1.96', '-1.96 <= Z < -1.645', '-1.645 <= Z < 0', '0 <= Z < 1.645', '1.645 <= Z < 1.96']
panels = {'Panel A (cutoff 1.96)': ([(-inf, -1.96), (-1.96, 0), (0, 1.96), (1.96, inf)], ['Z < -1.96', '-1.96 <= Z < 0', '0 <= Z < 1.96', 'Z >= 1.96 (reference)']),
          'Panel B (cutoff 1.645)': ([(-inf, -1.645), (-1.645, 0), (0, 1.645), (1.645, inf)], ['Z < -1.645', '-1.645 <= Z < 0', '0 <= Z < 1.645', 'Z >= 1.645 (reference)']),
          'Panel C (both cutoffs)': ([(-inf, -1.96), (-1.96, -1.645), (-1.645, 0), (0, 1.645), (1.645, 1.96), (1.96, inf)], names + ['Z >= 1.96 (reference)'])}
t1 = pd.concat([pd.DataFrame({'panel': p, 'region': n, 'estimates (trials)': cnt(b, c)}) for p, (c, n) in panels.items()])
t1.to_csv('output/table_1_counts.csv', index=False); print(t1.to_string(index=False))
