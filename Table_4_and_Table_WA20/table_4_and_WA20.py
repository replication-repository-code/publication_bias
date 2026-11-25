"""Tables 4 and WA20: types of corrections for publication bias, what they assume, and the corrected effects on financial behavior
and downstream outcomes (table 4) and on financial knowledge (table WA20, web appendix R), combined data.
The files in data/ are copies of the outputs of other folders (see the README of the package); this script only collects them.
The Andrews and Kasy models, the 3PSM, and the 4PSM use t-distributed true effects unless the row says otherwise. For the Andrews
and Kasy models, the table gives one line per choice of cutoffs (Z = 1.645, 1.96, and both). For knowledge, the GMM row uses the
version of equation 2 that pools all negative estimates in one region (only 4 knowledge estimates are significantly negative).
Intervals: 95% confidence intervals (credible intervals for RoBMA-PSMA).
Run from this folder: python3 table_4_and_WA20.py   (writes output/table_4.csv and output/table_WA20.csv, one line per row of the
table and, for the Andrews and Kasy models, per cutoff; the column cell is the text of the table)"""
import os
import pandas as pd

os.makedirs('output', exist_ok=True)
rd = lambda f: pd.read_csv(f'data/{f}')
Z = 1.959964
DOM = {'behavior': 'C_full', 'knowledge': 'C_know'}
CUTS = [('1.645', '1645'), ('1.96', '196'), ('Both', 'both')]


def ak(file, model, eq):
    """One line per cutoff: Andrews and Kasy maximum likelihood, equation 1 (sym) or 2 (asym)."""
    d = rd(file); d = d[d.model == model].set_index('spec')
    return [(lab, x.mu, x.mu - Z * x.se_mu, x.mu + Z * x.se_mu) for lab, c in CUTS for x in [d.loc[eq + c]]]


def gmm(dom):
    d = rd('gmm.csv'); d = d[d['sample'] == DOM[dom]].set_index('spec')
    end = '' if dom == 'behavior' else '_3'                       # knowledge: three-region equation 2
    return [(lab, x.gmm_g, x.gmm_ci_lo, x.gmm_ci_hi) for lab, c in CUTS for x in [d.loc['asym' + c + end]]]


def selection(dom, model):
    x = rd('selection_models.csv')
    x = x[(x['sample'] == DOM[dom]) & (x.selection == model) & (x.true_effects == 't (Andrews and Kasy)')].iloc[0]
    return [('', x.g, x.lo, x.hi)]


def other(dom, method):
    x = rd('corrections.csv')
    x = x[(x['sample'] == DOM[dom]) & (x.method == method)].iloc[0]
    return [('', x.g, x.lo, x.hi)]


def random_effects(dom):
    x = rd(f'random_effects_{dom}.csv').iloc[0]
    return [('', x.mu, x.lo, x.hi)]


def robma(dom, version):
    x = rd('robma.csv').set_index('model').loc[f'{DOM[dom]}_{version}_default']
    return [('', x.mu, x.mu_lo, x.mu_hi)]


ROWS = [  # correction, publication may depend on, distribution assumed, estimates
    ('Random-effects model', 'Not modeled', 'Normal', random_effects),
    ('Andrews and Kasy, equation 1', 'Significance, not sign', 't', lambda d: ak(f'andrews_kasy_t_{d}.csv', 't', 'sym')),
    ('Andrews and Kasy, equation 2', 'Sign and significance', 't', lambda d: ak(f'andrews_kasy_t_{d}.csv', 't', 'asym')),
    ('Andrews and Kasy, equation 2', 'Sign and significance', 'Normal',
     lambda d: ak(f'andrews_kasy_t_and_normal_{d}.csv', 'normal', 'asym')),
    ('Andrews and Kasy, equation 2 nonparametric (GMM)', 'Sign and significance', 'Unrestricted', gmm),
    ('Three-parameter selection model (3PSM)', 'Positive and significant, or not', 't', lambda d: selection(d, '3PSM')),
    ('Four-parameter selection model (4PSM)', 'Sign, and significance of positive estimates', 't', lambda d: selection(d, '4PSM')),
    ('PET', 'Not modeled; bias rises with SE', 'Not specified', lambda d: other(d, 'PET')),
    ('PEESE', 'Not modeled; bias rises with SE²', 'Not specified', lambda d: other(d, 'PEESE')),
    ('UWLS', 'Not modeled', 'Not specified', lambda d: other(d, 'UWLS')),
    ('WAAP', 'Not modeled', 'Not specified', lambda d: other(d, 'WAAP')),
    ('RoBMA-PSMA, multilevel', 'Six selection models, PET, PEESE, or none', 'Normal', lambda d: robma(d, 'ml')),
    ('RoBMA-PSMA, trial averages', 'Six selection models, PET, PEESE, or none', 'Normal', lambda d: robma(d, 'avg')),
]


def f3(v):
    s = f'{v:.3f}'
    return '−' + s[1:] if s.startswith('-') else s


pd.set_option('display.width', 200); pd.set_option('display.max_colwidth', 60)
for dom, name in [('behavior', 'table_4'), ('knowledge', 'table_WA20')]:
    out = []
    for k, (corr, dep, dist, fun) in enumerate(ROWS, 1):
        for cut, g, lo, hi in fun(dom):
            cell = (f'{cut}: ' if cut else '') + f'{f3(g)} ({f3(lo)} to {f3(hi)})'
            out.append(dict(row=k, correction=corr, publication_may_depend_on=dep, distribution_assumed=dist, cutoffs=cut,
                            g=round(g, 4), ci_lo=round(lo, 4), ci_hi=round(hi, 4), cell=cell))
    t = pd.DataFrame(out)
    t.to_csv(f'output/{name}.csv', index=False)
    print(f'{name} ({dom})'); print(t[['correction', 'distribution_assumed', 'cell']].to_string(index=False)); print()
