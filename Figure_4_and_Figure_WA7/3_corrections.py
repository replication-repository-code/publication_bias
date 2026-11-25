"""Precision-based corrections (UWLS, WAAP, PET, PEESE, PET-PEESE) and trim-and-fill, combined data (all outcomes) and Kaiser et
al.'s data alone, behavior and knowledge.
Same estimators as Table_WA10/1_corrections.R, without R packages (corrections_functions.py): weights 1/SE^2,
CR2 cluster-robust standard errors by trial with Satterthwaite degrees of freedom (as clubSandwich 0.7.0 for weighted lm fits);
trim-and-fill (L0) on the two-level REML random-effects model (as metafor 4.4.0). For behavior (C_full) and knowledge (C_know) the
script reproduces Table_WA10/output/corrections.csv to six decimals. Run from this folder: python3 3_corrections.py   (a few seconds)"""
import pandas as pd
from corrections_functions import precision_methods, trimfill

rows = []
for s in ['C_full', 'K_only', 'C_know', 'K_know']:
    d = pd.read_csv(f'data/{s}.csv')
    for k, v in precision_methods(d.yi, d.sei, d.cluster).items():
        rows.append(dict(sample=s, method=k, **v))
    tf = trimfill(d.yi, d.sei ** 2)
    rows.append(dict(sample=s, method='Trim-and-fill', g=tf['g'], se=tf['se'], lo=tf['lo'], hi=tf['hi'], k0=tf['k0'], side=tf['side']))
    rows.append(dict(sample=s, method='Two-level random effects (REML)', g=tf['g_before'], se=tf['se_before'], lo=tf['lo_before'],
                     hi=tf['hi_before']))
r = pd.DataFrame(rows)
# check against Table_WA10/output/corrections.csv
REF = {('C_full', 'UWLS'): 0.027573, ('C_full', 'WAAP'): 0.004759, ('C_full', 'PET'): 0.003817, ('C_full', 'PEESE'): 0.022369,
       ('C_full', 'Trim-and-fill'): 0.088624, ('C_know', 'UWLS'): 0.156613, ('C_know', 'WAAP'): 0.146936, ('C_know', 'PET'): 0.104911,
       ('C_know', 'PEESE'): 0.138282, ('C_know', 'Trim-and-fill'): 0.200607}
for (s, m), g in REF.items():
    assert abs(r[(r['sample'] == s) & (r.method == m)].g.iloc[0] - g) < 2e-6, (s, m)
r.to_csv('output/precision_and_trimfill.csv', index=False)
pd.set_option('display.width', 200)
print(r[['sample', 'method', 'g', 'se', 'lo', 'hi', 'df', 'n_estimates', 'n_trials', 'rule', 'k0', 'side']].round(4).to_string(index=False))
