"""Table WA3: excess mass of Z statistics just above the significance thresholds 1.645 and 1.96 in the combined data (Kaiser et
al.'s estimates plus the newer trials), for behavior and for knowledge.

Excess mass (Allen, Dechow, Pope, and Wu 2017, NBER WP 20343, following Chetty et al. 2011). Z statistics are grouped in bins of
width h whose edges fall on the threshold c, with h = 0.05, 0.10 or 0.40 and 8 bins on each side of c (windows of 0.4, 0.8 and 3.2).
The bunching region is the first bin above c, [c, c + h). A cubic polynomial is fit to the bin counts in the window, leaving out the
bunching region and the bin that covers the bunching region of the other conventional threshold (1.645 or 1.96) when it falls in the
window (the bin whose midpoint lies in [o, o + h)). As in Allen et al., the counterfactual is shifted so that it has the same number of
estimates over the window as the data. Excess mass =
(observed - counterfactual count in the bunching region) / counterfactual count. Standard error of the excess count: 500 bootstrap
samples that resample trials (clusters) with replacement; p one-sided (excess > 0) from the bootstrap t statistic. Each test draws
its own bootstrap samples (seed 1).
With bins of 0.40, the first bin above 1.645 also contains 1.96, so the two thresholds cannot be separated; that cell is not reported.
Placebo thresholds 1.30 and 2.30 (bins of 0.10), which have no conventional meaning, calibrate the test.
Sensitivity ("excl. Z at a threshold", bins of 0.10): estimates whose |Z| lies in [c, c + 0.002) for c = 1.645, 1.96, 2.326 or
2.576 are dropped. Most of them have standard errors derived from significance stars (the largest SE consistent with the reported
significance), which places Z at the threshold by construction.
Writes output/table_WA3.csv. Run from this folder: python3 table_WA3.py   (about 10 seconds)
"""
import numpy as np, pandas as pd
from scipy import stats

NB, PER_SIDE = 500, 8
CONV = [1.645, 1.96]
AT = [1.645, 1.96, 2.326, 2.576]


def excess(z, c, h):
    edges = c + h * np.arange(-PER_SIDE, PER_SIDE + 1)
    cnt, _ = np.histogram(z, edges); mid = (edges[:-1] + edges[1:]) / 2
    bun = (mid >= c) & (mid < c + h)
    excl = bun.copy()
    for o in CONV:                                         # the bin covering the other conventional threshold's bunching region
        if o != c: excl |= (mid >= o) & (mid < o + h)
    x = mid - c
    coef = np.polyfit(x[~excl], cnt[~excl], 3); cf = np.polyval(coef, x)
    cf = cf + (cnt.sum() - cf.sum()) / len(cf)            # upward shift: same mass over the window
    return cnt[bun].sum(), cf[bun].sum()


def test(z, cl, c, h):
    obs, exp = excess(z, c, h)
    rng = np.random.default_rng(1); ids = np.unique(cl); idx = {k: np.where(cl == k)[0] for k in ids}; boot = []
    for _ in range(NB):
        s = rng.choice(ids, len(ids), replace=True); zz = np.concatenate([z[idx[k]] for k in s])
        o, e = excess(zz, c, h); boot.append(o - e)
    se = np.std(boot, ddof=1)
    return dict(observed=int(obs), counterfactual=exp, excess_estimates=obs - exp, se_boot=se, excess_mass=(obs - exp) / exp,
                p_one_sided=stats.norm.sf((obs - exp) / se))


rows = []
for f, dom in [('C_full', 'Behavior'), ('C_know', 'Knowledge')]:
    d = pd.read_csv(f'data/{f}.csv')
    za = (d.yi / d.sei).abs(); at = np.zeros(len(d), bool)
    for c0 in AT: at |= ((za >= c0) & (za < c0 + 0.002)).to_numpy()
    for lab, x in [('Combined', d), ('Combined, excl. Z at a threshold', d[~at])]:
        z = (x.yi / x.sei).to_numpy(); cl = x.cluster.astype(str).to_numpy()
        specs = [(c, h) for c in CONV for h in (0.05, 0.10, 0.40)] if lab == 'Combined' else [(c, 0.10) for c in CONV]
        if lab == 'Combined': specs += [(1.30, 0.10), (2.30, 0.10)]
        for c, h in specs:
            if c == 1.645 and h == 0.40: continue          # the bin [1.645, 2.045) also contains 1.96
            rows.append(dict(outcome=dom, data=lab, threshold=c, bin_width=h, window=PER_SIDE * h,
                             bunching_region=f'[{c}, {c + h:.3f})', **test(z, cl, c, h)))
    print(f'{dom}: {at.sum()} estimates with |Z| within 0.002 above a conventional threshold dropped in the sensitivity analysis')
res = pd.DataFrame(rows)
res.to_csv('output/table_WA3.csv', index=False)
pd.set_option('display.width', 200)
print(res.assign(excess_mass=(100 * res.excess_mass).round(1), counterfactual=res.counterfactual.round(1),
                 p_one_sided=res.p_one_sided.round(3)).drop(columns=['excess_estimates', 'se_boot']).to_string(index=False))
