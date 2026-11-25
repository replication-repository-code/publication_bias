"""Summaries of the specification curves of figures 4 and WA7 (corrected estimates only; the uncorrected estimates are benchmarks),
quoted in the body of the paper and in web appendix N. Run from this folder: python3 6_summaries.py   (writes
output/spec_curve_summary.csv, output/spec_curve_paper_estimates_rank.csv, and output/spec_curve_variance_decomposition.csv)"""
import os, itertools, math
import numpy as np, pandas as pd

AK, GMM = 'Andrews and Kasy (maximum likelihood)', 'Andrews and Kasy (GMM)'
d = pd.read_csv('output/spec_curve_estimates.csv'); OUT = 'output'
MAIN = {'behavior': 'C_full', 'knowledge': 'C_know'}
unc = d[d.method == 'Random effects, estimates nested in trials'].set_index('sample').estimate
C = d[d.corrected].copy()
C['family'] = np.select(
    [(C.method == AK) & (C.true_effects == 't'), (C.method == AK) & (C.true_effects == 'Normal'), C.method == GMM,
     C.method.isin(['3PSM', '4PSM']), C.method.isin(['PET', 'PEESE', 'UWLS', 'WAAP'])],
    ['Andrews and Kasy, t true effects (ML)', 'Andrews and Kasy, normal true effects (ML)', 'Andrews and Kasy, GMM',
     '3PSM and 4PSM, t true effects', 'Precision-based (PET, PEESE, UWLS, WAAP)'], default=C.method)
C['eq'] = C.selection.map({'sym1645': 'Eq. 1', 'sym196': 'Eq. 1', 'symboth': 'Eq. 1', 'asym1645': 'Eq. 2', 'asym196': 'Eq. 2',
                           'asymboth': 'Eq. 2'})
u = C['sample'].map(unc)
C['below_uncorrected'] = C.estimate < u
C['below_half_uncorrected'] = C.estimate < 0.5 * u
C['interval_excludes_zero'] = C.ci_lo > 0
C['interval_below_uncorrected'] = C.ci_hi < u          # the whole interval lies below the uncorrected estimate of the same data
C['interval_below_0.10'] = C.ci_hi < 0.10              # excludes Kaiser et al.'s (2022a) 0.10


def summ(x):
    return pd.Series(dict(specifications=len(x), median=x.estimate.median(), q25=x.estimate.quantile(.25), q75=x.estimate.quantile(.75),
                          min=x.estimate.min(), max=x.estimate.max(), share_interval_excludes_zero=x.interval_excludes_zero.mean(),
                          share_below_uncorrected=x.below_uncorrected.mean(), share_below_half_uncorrected=x.below_half_uncorrected.mean(),
                          share_interval_below_uncorrected=x.interval_below_uncorrected.mean(),
                          share_interval_below_0_10=x['interval_below_0.10'].mean()))


# share of the variance of the corrected estimates explained by each choice: Shapley decomposition of the R^2 of a regression of
# the estimates on indicators (main effects only; 1 - R^2 is the part due to interactions)
def r2(y, cols, x):
    if not cols: return 0.0
    X = pd.get_dummies(x[cols].astype(str), drop_first=True).astype(float)
    X.insert(0, 'const', 1.0)
    b, *_ = np.linalg.lstsq(X.values, y, rcond=None)
    return 1 - np.sum((y - X.values @ b) ** 2) / np.sum((y - y.mean()) ** 2)


def shapley(y, fac, x):
    return {f: sum(math.factorial(k) * math.factorial(len(fac) - k - 1) / math.factorial(len(fac)) *
                   (r2(y, list(S) + [f], x) - r2(y, list(S), x))
                   for k in range(len(fac)) for S in itertools.combinations([g for g in fac if g != f], k)) for f in fac}


pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
rows, pr, vd = [], [], []
for o, x in C.groupby('outcome'):
    rows.append(dict(outcome=o, subset='All corrected specifications', **summ(x)))
    rows.append(dict(outcome=o, subset='Without Andrews and Kasy with normal true effects',
                     **summ(x[x.family != 'Andrews and Kasy, normal true effects (ML)'])))
    if x['sample'].nunique() > 1:
        for s, y in x.groupby('sample', sort=False):
            rows.append(dict(outcome=o, subset=f'Data: {s} (uncorrected g = {unc[s]:.3f})', **summ(y)))
    for f, y in x.groupby('family', sort=False):
        rows.append(dict(outcome=o, subset=f'Correction: {f}', **summ(y)))
    for e, y in x.groupby('eq'):
        rows.append(dict(outcome=o, subset=f'Andrews and Kasy (ML and GMM), {e}', **summ(y)))
    # where the estimates reported in the body of the paper fall among all corrected specifications
    m = x[(x['sample'] == MAIN[o]) & (x.method == AK) & (x.true_effects == 't') & x.selection.isin(['asym1645', 'asym196', 'asymboth']) &
          (x.unit == 'estimates')]
    for r in m.itertuples():
        pr.append(dict(outcome=o, specification=r.variant, estimate=r.estimate, rank=int((x.estimate < r.estimate).sum() + 1),
                       of=len(x), percentile=100 * (x.estimate < r.estimate).mean()))
    # variance decomposition
    if x['sample'].nunique() > 1:
        xx = x.assign(correction=x.method + ' | ' + x.variant.fillna(''))
        y = xx.estimate.values; fac = ['sample', 'correction']; sh = shapley(y, fac, xx); full = r2(y, fac, xx)
        for f, lab in zip(fac, ['Data', 'Correction']):
            vd.append(dict(outcome=o, set='All corrected specifications', specifications=len(xx), factor=lab, r2_alone=r2(y, [f], xx),
                           shapley=sh[f], r2_all_factors=full))
    mm = x[(x.method == AK) & (x.unit == 'estimates')]; y = mm.estimate.values
    fac = (['sample'] if mm['sample'].nunique() > 1 else []) + ['selection', 'true_effects']
    sh = shapley(y, fac, mm); full = r2(y, fac, mm)
    for f in fac:
        vd.append(dict(outcome=o, set='Andrews and Kasy (maximum likelihood)', specifications=len(mm),
                       factor={'sample': 'Data', 'selection': 'Selection function', 'true_effects': 'True effects (t or normal)'}[f],
                       r2_alone=r2(y, [f], mm), shapley=sh[f], r2_all_factors=full))
tab, pr, vd = pd.DataFrame(rows), pd.DataFrame(pr), pd.DataFrame(vd)
tab.to_csv(os.path.join(OUT, 'spec_curve_summary.csv'), index=False)
pr.to_csv(os.path.join(OUT, 'spec_curve_paper_estimates_rank.csv'), index=False)
vd.to_csv(os.path.join(OUT, 'spec_curve_variance_decomposition.csv'), index=False)
print(tab.round(3).to_string(index=False)); print(pr.round(3).to_string(index=False))
print(vd.round(3).to_string(index=False))
