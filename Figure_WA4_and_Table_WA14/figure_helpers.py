"""Plotting helpers shared by the figure folders (matplotlib)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.transforms as mt, matplotlib.ticker
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
INK, MUTED, RED, BLUE, GRID = '#222222', '#6B6B6B', '#C0392B', '#1F5AA6', '#E6E6E6'
z95, z50 = 1.959964, 0.674490


def z_histogram(d, out, w, edge=2.0):
    """Stacked histogram of Z = g / SE: Kaiser et al. (gray) and newer trials (blue). Bins of width w, one of whose edges is at Z =
    edge, closed on the right as in R and JASP; the edges run from the last one below -4 to the first one above 8, and values beyond
    them go to the end bins.
    Figures 2 and WA5: w = 0.40 and an edge at 2.00, so that one bin runs from 1.60 to 2.00 and contains both positive critical
    values; with bins of 0.40, that placement is the one at which the bin containing Z = 1.96 stands out most from the two bins
    next to it, in Kaiser et al.'s data and in the combined data alike (estimates of effects on behavior; starting points 0.01
    apart).
    Returns a line with the counts in the bin that contains 1.96 and in the two bins next to it."""
    z = (d.yi / d.sei).values; kai = d.source.str.startswith('Kaiser').values
    k0, k1 = int(np.floor((-4 - edge) / w + 1e-9)), int(np.ceil((8 - edge) / w - 1e-9))
    e = edge + w * np.arange(k0, k1 + 1)                     # bin i is (e[i], e[i + 1]]
    def counts(x):                                           # the tolerance keeps a value that sits on an edge (up to rounding) in the bin below it
        i = np.ceil((x - e[0]) / w - 1e-9).astype(int) - 1
        return np.bincount(np.clip(i, 0, len(e) - 2), minlength=len(e) - 1)
    ck, cu = counts(z[kai]), counts(z[~kai])
    n_new_trials = d[~kai].cluster.str.startswith('U-').groupby(d[~kai].cluster).first().sum()
    lw = 0.2 if w < 0.1 else 0.3 if w < 0.2 else 0.5
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.bar(e[:-1], ck, width=w, align='edge', color=MUTED, edgecolor='white', linewidth=lw,
           label=f'Kaiser et al. (2022) ({kai.sum()} estimates, {d[kai].cluster.nunique()} studies)')
    ax.bar(e[:-1], cu, width=w, bottom=ck, align='edge', color=BLUE, edgecolor='white', linewidth=lw,
           label=f'New RCTs, 2019–2026 ({(~kai).sum()} estimates, {n_new_trials} new trials)')
    for x in (-1.96, -1.645, 1.645, 1.96): ax.axvline(x, color=INK, lw=0.8, ls='--' if abs(x) < 1.9 else '-')
    ax.set_xlabel('Z statistic (estimate / standard error)'); ax.set_ylabel('Number of estimates')
    ax.set_xlim(e[0], e[-1]); ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=8.5)
    fig.tight_layout(); fig.savefig(out, dpi=300); plt.close(fig)
    s, t = int(np.searchsorted(e, 1.96) - 1), ck + cu       # s: the bin that contains 1.96
    f = lambda x: f'{x:.3f}'[:-1] if f'{x:.3f}'.endswith('0') else f'{x:.3f}'
    return (f'bins of {w:.2f}: the bin ({f(e[s])}, {f(e[s + 1])}] holds {t[s]} estimates (Kaiser et al. {ck[s]}); the bins next to it '
            f'hold {t[s - 1]} and {t[s + 1]} (Kaiser et al. {ck[s - 1]} and {ck[s + 1]}); the tallest other bin holds {np.delete(t, s).max()} '
            f'(Kaiser et al. {np.delete(ck, s).max()})')


STY = {'unc': (INK, 'o'), 'sym': (RED, 's'), 'asym': (BLUE, '^'), 'rob': (BLUE, 'D')}


def forest(rows, heads, path, xlim, xticks, h=4.0, legend=False):
    """rows: (label, estimate, lo95, hi95, lo50, hi50, style); style 'na' = estimate not available."""
    fig, ax = plt.subplots(figsize=(6.8, h)); y = 0; pos = []
    tr = lambda: mt.blended_transform_factory(ax.transAxes, ax.transData)
    for i, r in enumerate(rows):
        if i in heads:
            y -= 0.4; ax.text(-0.02, y, heads[i], fontsize=8.5, fontweight='bold', color=INK, va='center', ha='right', transform=tr()); y -= 0.85
        pos.append(y); y -= 1
    for (lab, e, l95, h95, l50, h50, s), yy in zip(rows, pos):
        if s == 'na':
            ax.text(0.01, yy, 'corrected estimate not estimable', fontsize=7.5, style='italic', va='center', ha='left', color=MUTED,
                    bbox=dict(facecolor='white', edgecolor='none', pad=1), zorder=4, transform=tr())
            ax.text(1.01, yy, '—', fontsize=8, va='center', ha='left', color=MUTED, transform=tr()); continue
        col, m = STY[s]
        ax.plot([l95, h95], [yy, yy], color=col, lw=1.2, solid_capstyle='butt'); ax.plot([l50, h50], [yy, yy], color=col, lw=3.6, solid_capstyle='butt')
        ax.plot(e, yy, marker=m, color=col, ms=7, mec='white', mew=0.8, zorder=3)
        ax.text(1.01, yy, f'{e:.3f} [{l95:.3f}, {h95:.3f}]'.replace('-', '−'), fontsize=8, va='center', ha='left', color=INK, transform=tr())
    ax.set_yticks(pos); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.axvline(0, color=MUTED, lw=0.8); ax.set_xlim(*xlim); ax.set_xticks(xticks)
    ax.xaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    ax.set_xlabel("Average effect (Hedges' g)"); ax.set_ylim(min(pos) - 0.8, 0.2)
    ax.text(1.01, 1.0, 'Estimate [95% interval]', transform=ax.transAxes, fontsize=8, ha='left', va='bottom', color=MUTED)
    if legend:
        ax.legend(handles=[Line2D([], [], marker='o', color=INK, ls='', label='Uncorrected random effects'),
                           Line2D([], [], marker='^', color=BLUE, ls='', label='Corrected (Andrews–Kasy Eq. 2, Z ≥ 1.645)')],
                  loc='lower center', bbox_to_anchor=(0.5, 1.01), frameon=False, fontsize=8, ncol=1)
    fig.tight_layout(); fig.subplots_adjust(right=0.74); fig.savefig(path, dpi=300); plt.close(fig)


def forest_by_method(re_csv, ak_csv, robma_csv, out, xlim, xticks):
    """Figure 3 / Figure WA6: uncorrected random effects, Andrews-Kasy equations 1 and 2 (t-distributed true effects), multilevel RoBMA-PSMA."""
    r = pd.read_csv(re_csv).iloc[0]; ak = pd.read_csv(ak_csv).set_index('spec'); rb = pd.read_csv(robma_csv).iloc[0]
    rows = [('Uncorrected random effects', r.mu, r.lo, r.hi, r.mu - z50 * r.se_CR2, r.mu + z50 * r.se_CR2, 'unc')]
    for sp, lab, st in [('sym1645', 'Eq. 1 symmetric, |Z| ≥ 1.645', 'sym'), ('sym196', 'Eq. 1 symmetric, |Z| ≥ 1.96', 'sym'),
                        ('symboth', 'Eq. 1 symmetric, both cutoffs', 'sym'), ('asym1645', 'Eq. 2 asymmetric, Z ≥ 1.645', 'asym'),
                        ('asym196', 'Eq. 2 asymmetric, Z ≥ 1.96', 'asym'), ('asymboth', 'Eq. 2 asymmetric, both cutoffs', 'asym')]:
        x = ak.loc[sp]; rows.append((lab, x.mu, x.mu - z95 * x.se_mu, x.mu + z95 * x.se_mu, x.mu - z50 * x.se_mu, x.mu + z50 * x.se_mu, st))
    rows.append(('Multilevel RoBMA-PSMA', rb.mu, rb.mu_lo, rb.mu_hi, rb.mu_q25, rb.mu_q75, 'rob'))
    forest(rows, {0: 'Uncorrected', 1: 'Andrews–Kasy (2019), t true effects', 7: 'RoBMA-PSMA'}, out, xlim, xticks)


def subgroup_forest(csv, out, domain):
    """Figure WA4: uncorrected and corrected (equation 2) estimates by subgroup."""
    d = pd.read_csv(csv); rows = []; heads = {}; last = None
    names = {'Outcome': 'Outcome', 'Intervention': 'Intervention type', 'Sample': 'Sample', 'Journal': 'Journal', 'All': 'All'}
    for x in d.itertuples():
        cat = x.group.split(':')[0]
        lab = (x.group.split(': ')[-1].capitalize() if ':' in x.group else 'All estimates').replace('Other behaviors and downstream outcomes', 'Other outcomes')
        if cat != last: heads[len(rows)] = names[cat]; last = cat
        rows.append((f'{lab} (k = {x.k:,})', x.re_mu, x.re_lo, x.re_hi, x.re_mu - z50 * x.re_se, x.re_mu + z50 * x.re_se, 'unc'))
        if x.converged and not np.isnan(x.ak_mu):
            rows.append(('', x.ak_mu, x.ak_mu - z95 * x.ak_se, x.ak_mu + z95 * x.ak_se, x.ak_mu - z50 * x.ak_se, x.ak_mu + z50 * x.ak_se, 'asym'))
        else:
            rows.append(('', np.nan, np.nan, np.nan, np.nan, np.nan, 'na'))
    xl, xt = ((-0.12, 0.26), [-0.1, 0, 0.1, 0.2]) if domain == 'behavior' else ((-0.1, 0.48), [0, 0.2, 0.4])
    forest(rows, heads, out, xl, xt, h=0.2 * len(rows) + 0.3 * len(heads) + 1.2, legend=True)
