"""Web appendix C: the excess mass of Z statistics just above 1.96 (behavior, combined data) without significance-derived standard
errors, over 80 specifications of the test.

Data: data/C_full_se_source.csv, the behavior sample (C_full) with se_from_significance = 1 for the 56 estimates of the newer
trials whose standard error was derived from a reported p-value or from significance stars (from our coding notes: the
statistical information and input values used for each effect size, and the October 2026 change log). Kaiser et al.'s file does not record how each standard error was obtained; an estimate whose |Z| lies within 0.002 of
1.645, 1.96, 2.326, 2.576, or 3.291 carries the mark of a standard error derived from a significance level and is also dropped.
Samples: all estimates; without |Z| in [c, c + 0.002) for c = 1.645, 1.96, 2.326, 2.576 (the sensitivity analysis of Table WA3);
without significance-derived standard errors (both rules above).
Test: as in table_WA3.py (Allen et al. 2017): bins with an edge at 1.96, bunching region [1.96, 1.96 + h), polynomial fitted to the
other bins of the window (leaving out the bin of the other threshold's bunching region), counterfactual shifted to the observed
mass, standard error from 500 bootstrap samples of trials (seed 1), one-sided p. Grid: h = 0.05, 0.10, 0.15, 0.20, 0.25; degree 1-4;
6, 8, 10, or 12 bins on each side.
Writes output/excess_mass_robustness_grid.csv and output/excess_mass_robustness_summary.csv. Run: python3 excess_mass_robustness.py
"""
import numpy as np, pandas as pd
from scipy import stats

NB, CONV, C0 = 500, [1.645, 1.96], 1.96


def excess(z, c, h, per_side, deg):
    edges = c + h * np.arange(-per_side, per_side + 1)
    cnt, _ = np.histogram(z, edges); mid = (edges[:-1] + edges[1:]) / 2
    bun = (mid >= c) & (mid < c + h); excl = bun.copy()
    for o in CONV:
        if o != c: excl |= (mid >= o) & (mid < o + h)
    coef = np.polyfit(mid[~excl] - c, cnt[~excl], deg); cf = np.polyval(coef, mid - c)
    cf = cf + (cnt.sum() - cf.sum()) / len(cf)
    return cnt[bun].sum(), cf[bun].sum()


def test(z, cl, c, h, per_side, deg):
    o, e = excess(z, c, h, per_side, deg)
    rng = np.random.default_rng(1); ids = np.unique(cl); idx = {k: np.where(cl == k)[0] for k in ids}; boot = []
    for _ in range(NB):
        s = rng.choice(ids, len(ids), replace=True); oo, ee = excess(np.concatenate([z[idx[k]] for k in s]), c, h, per_side, deg)
        boot.append(oo - ee)
    se = np.std(boot, ddof=1)
    return dict(observed=int(o), counterfactual=e, excess_mass=(o - e) / e, p_one_sided=stats.norm.sf((o - e) / se))


d = pd.read_csv('data/C_full_se_source.csv'); az = (d.yi / d.sei).abs()
at_paper = np.any([((az >= c) & (az < c + 0.002)) for c in [1.645, 1.96, 2.326, 2.576]], axis=0)
at_crit = np.any([np.abs(az - c) < 0.002 for c in [1.645, 1.96, 2.326, 2.576, 3.291]], axis=0)
samples = {'All estimates': np.zeros(len(d), bool), 'Without |Z| at a threshold (Table WA3)': at_paper,
           'Without significance-derived standard errors': d.se_from_significance.astype(bool).to_numpy() | at_crit}
rows = []
for lab, drop in samples.items():
    x = d[~drop]; z = (x.yi / x.sei).to_numpy(); cl = x.cluster.astype(str).to_numpy()
    print(f'{lab}: {len(x)} estimates ({int(drop.sum())} dropped)')
    for h in (0.05, 0.10, 0.15, 0.20, 0.25):
        for deg in (1, 2, 3, 4):
            for ps in (6, 8, 10, 12):
                rows.append(dict(sample=lab, estimates=len(x), bin_width=h, degree=deg, bins_per_side=ps, **test(z, cl, C0, h, ps, deg)))
g = pd.DataFrame(rows); g.to_csv('output/excess_mass_robustness_grid.csv', index=False)
s = g.groupby('sample', sort=False).apply(lambda x: pd.Series(dict(
    estimates=x.estimates.iloc[0], specifications=len(x), positive=int((x.excess_mass > 0).sum()), p_below_10=int((x.p_one_sided < .10).sum()),
    p_below_05=int((x.p_one_sided < .05).sum()), min_excess=x.excess_mass.min(), median_excess=x.excess_mass.median(), max_excess=x.excess_mass.max(),
    default_excess=x[(x.bin_width == 0.10) & (x.degree == 3) & (x.bins_per_side == 8)].excess_mass.iloc[0],
    default_p=x[(x.bin_width == 0.10) & (x.degree == 3) & (x.bins_per_side == 8)].p_one_sided.iloc[0]))).reset_index()
for c in ('estimates', 'specifications', 'positive', 'p_below_10', 'p_below_05'): s[c] = s[c].astype(int)
s.to_csv('output/excess_mass_robustness_summary.csv', index=False)
pd.set_option('display.width', 200); print(s.round(3).to_string(index=False))
