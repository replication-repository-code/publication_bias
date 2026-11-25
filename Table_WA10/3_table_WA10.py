"""Table WA10: additional publication-bias corrections, financial behavior and downstream outcomes, combined data (web appendix J).
Collects output/corrections.csv (1_corrections.R) and, for the 3PSM and 4PSM, which the table reports with t-distributed true
effects (Andrews and Kasy's t specification with the steps of these models; standard errors clustered by trial),
output/selection_models.csv (2_selection_models.R). Run from this folder after the two R scripts: python3 3_table_WA10.py
(writes output/table_WA10.csv)"""
import pandas as pd

c = pd.read_csv('output/corrections.csv'); c = c[c['sample'] == 'C_full'].set_index('method')
s = pd.read_csv('output/selection_models.csv')
s = s[(s['sample'] == 'C_full') & (s.true_effects == 't (Andrews and Kasy)')].set_index('selection')


def f3(v):
    t = f'{v:.3f}'
    return '−' + t[1:] if t.startswith('-') else t


rows = []
for label, x in [(m, c.loc[m]) for m in ['Uncorrected random-effects model', 'UWLS', 'WAAP', 'PET', 'PEESE', 'PET-PEESE']] + \
        [('Three-parameter selection model (3PSM)', s.loc['3PSM']), ('Four-parameter selection model (4PSM)', s.loc['4PSM']),
         ('Trim-and-fill', c.loc['Trim-and-fill'])]:
    rows.append(dict(method=label, g=x.g, ci_lo=x.lo, ci_hi=x.hi, cell=f'{f3(x.g)} [{f3(x.lo)}, {f3(x.hi)}]'))
x = c.loc['Uncorrected, leave one trial out (range)']
rows.append(dict(method='Uncorrected, leaving out one trial at a time (range)', g=float('nan'), ci_lo=x.lo, ci_hi=x.hi,
                 cell=f'{f3(x.lo)} to {f3(x.hi)}'))
t = pd.DataFrame(rows)
t.to_csv('output/table_WA10.csv', index=False)
print(t[['method', 'cell']].to_string(index=False))
