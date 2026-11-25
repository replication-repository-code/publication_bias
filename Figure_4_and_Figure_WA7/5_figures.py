"""Figures 4 and WA7: specification curves of the corrected average effect on financial behavior and downstream outcomes (figure 4)
and on financial knowledge (figure WA7, web appendix R), combined data and Kaiser et al.'s data alone.
Top panel: every corrected estimate, ranked, with its 95% interval (credible interval for RoBMA-PSMA); dashed and dotted lines:
uncorrected random-effects estimates of the combined data and of Kaiser et al.'s data alone. Bottom panel: the data and the
correction model of each specification (one row per model, in blocks of selection models, precision-based models, and other
corrections). Run from this folder: python3 5_figures.py   (writes output/figure_4.png and output/figure_WA7.png)"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False})
INK, MUTED, GRID, CI = '#222222', '#6B6B6B', '#E6E6E6', '#B8B8B8'
AK, GMM = 'Andrews and Kasy (maximum likelihood)', 'Andrews and Kasy (GMM)'

d = pd.read_csv('output/spec_curve_estimates.csv')
MAIN = {'behavior': 'C_full', 'knowledge': 'C_know'}
KAISER = {'behavior': 'K_only', 'knowledge': 'K_know'}
DATA = {'behavior': [('C_full', 'Combined data, all outcomes'), ('K_only', "Kaiser et al.'s data alone")],
        'knowledge': [('C_know', 'Combined data, all knowledge'), ('K_know', "Kaiser et al.'s data alone")]}
YL = {'behavior': ((-0.15, 0.25), [-0.1, 0, 0.1, 0.2]), 'knowledge': ((-0.35, 0.40), [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4])}


CUT = [('1645', '1.645'), ('196', '1.96'), ('both', 'both')]


def model_rows(x):
    """One row per correction model: selection models, precision-based models, other corrections."""
    sel, prec, other = [], [], []
    for eq, pre in (('Eq. 2', 'asym'), ('Eq. 1', 'sym')):
        for te, lab in (('t', 't'), ('Normal', 'normal')):
            for c, cl in CUT:
                sel.append(((x.method == AK) & (x.selection == pre + c) & (x.true_effects == te) & (x.unit == 'estimates'), f'{eq}, {cl}, {lab}'))
            if eq == 'Eq. 2' and te == 't':
                for c, cl in CUT:
                    sel.append(((x.method == AK) & (x.selection == pre + c) & (x.unit == 'trials'), f'{eq}, {cl}, t, each trial counted once'))
        for c, cl in CUT:
            sel.append(((x.method == GMM) & (x.selection == pre + c), f'{eq}, {cl}, GMM'))
    sel += [((x.method == '3PSM'), '3PSM'), ((x.method == '4PSM'), '4PSM')]
    prec = [((x.method == m), m) for m in ('PET', 'PEESE', 'UWLS', 'WAAP')]
    other = [((x.method == 'Trim-and-fill'), 'Trim-and-fill')]
    for v, lab in (('multilevel', 'RoBMA-PSMA, multilevel'), ('trial averages', 'RoBMA-PSMA, trial averages'),
                   ('trial averages, near-zero null prior', 'RoBMA-PSMA, trial averages, null prior'),
                   ('trial averages, effect prior centered on a positive effect', 'RoBMA-PSMA, trial averages, effect prior')):
        other.append(((x.method == 'RoBMA-PSMA') & (x.variant == v), lab))
    keep = lambda rows: [(m, lab) for m, lab in rows if m.any()]
    out = [('Selection models', keep(sel)), ('Precision-based models', keep(prec)), ('Other corrections', keep(other))]
    assert sum(int(m.sum()) for _, rows in out for m, _ in rows) == len(x), 'every specification in exactly one model row'
    return out


def curve(outcome, path):
    x = d[(d.outcome == outcome) & d.corrected].sort_values(['estimate', 'method', 'variant']).reset_index(drop=True)
    n = len(x); pos = np.arange(1, n + 1)
    unc = d[(d.outcome == outcome) & (d.method == 'Random effects, estimates nested in trials')].set_index('sample').estimate
    groups = [('Data', [(x['sample'] == s, lab) for s, lab in DATA[outcome]])] + model_rows(x)
    nrow = sum(len(g[1]) + 1 for g in groups)
    RH = 0.118                                                   # inches per indicator row
    top_h, gap, bot_h = 2.75, 0.12, RH * nrow
    H = 0.15 + top_h + gap + bot_h + 0.55
    fig = plt.figure(figsize=(6.5, H))
    L, R = 0.355, 0.985
    ax = fig.add_axes([L, (0.55 + bot_h + gap) / H, R - L, top_h / H])
    ab = fig.add_axes([L, 0.55 / H, R - L, bot_h / H], sharex=ax)
    (y0, y1), yt = YL[outcome]
    lo = x.ci_lo.to_numpy(); hi = x.ci_hi.to_numpy(); loc, hic = np.clip(lo, y0, y1), np.clip(hi, y0, y1)
    ax.vlines(pos, loc, hic, color=CI, lw=0.9, zorder=1)
    cut_hi, cut_lo = hi > y1, lo < y0                            # intervals cut at the edge of the panel
    ax.scatter(pos[cut_hi], np.full(cut_hi.sum(), y1), marker=6, s=12, color=CI, linewidths=0, zorder=1, clip_on=False)
    ax.scatter(pos[cut_lo], np.full(cut_lo.sum(), y0), marker=7, s=12, color=CI, linewidths=0, zorder=1, clip_on=False)
    ax.axhline(0, color=MUTED, lw=0.8, zorder=0)
    lines = [(MAIN[outcome], (0, (4, 3)), 'Uncorrected, combined data'), (KAISER[outcome], (0, (1, 2)), "Uncorrected, Kaiser et al.'s data")]
    for s, ls, lab in lines:
        ax.axhline(unc[s], color=INK, lw=0.8, ls=ls, zorder=0)
        ax.text(1, unc[s] + 0.008 * (y1 - y0), f'{lab}: g = {unc[s]:.3f}', fontsize=7.3, color=INK, ha='left', va='bottom', zorder=2,
                bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
    ax.scatter(pos, x.estimate, s=15, color=INK, edgecolor='white', linewidths=0.5, zorder=3)
    ax.set_ylim(y0, y1); ax.set_yticks(yt); ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_yticklabels([f'{v:.1f}'.replace('-', '−') for v in yt], fontsize=7.5)
    ax.set_ylabel("Corrected average effect (Hedges' g)", fontsize=8)
    ax.tick_params(axis='x', labelbottom=False, length=0); ax.spines['bottom'].set_visible(False)
    ax.set_xlim(0.2, n + 0.8)
    # indicator panel: a bold heading row for each group, then one row per choice
    yy = 0; ticks, labels, bold = [], [], []
    for gname, rows in groups:
        ticks.append(yy); labels.append(gname); bold.append(True); yy -= 1
        first = yy
        for hit, lab in rows:
            hit = hit.to_numpy()
            ab.scatter(pos[hit], np.full(hit.sum(), yy), s=12, marker='o', color=INK, linewidths=0, zorder=2)
            ticks.append(yy); labels.append(lab); bold.append(False); yy -= 1
        ab.axhspan(yy + 0.5, first + 0.5, color='#F2F2F2', zorder=0, lw=0)
    ab.set_yticks(ticks); ab.set_yticklabels(labels, fontsize=7.1); ab.tick_params(axis='y', length=0, pad=3)
    for t, b in zip(ab.get_yticklabels(), bold):
        if b: t.set_fontweight('bold')
    ab.set_ylim(yy + 0.5, 0.6); ab.spines['left'].set_visible(False)
    step = 10 if n > 40 else 5
    ticks_x = [1] + [k for k in range(step, n + 1, step) if n - k >= step / 2.5] + [n]
    ab.set_xticks(ticks_x); ab.tick_params(axis='x', labelsize=7.5)
    ab.set_xlabel('Specifications, ranked by the corrected estimate', fontsize=8)
    fig.savefig(path, dpi=300); plt.close(fig)
    return n


for o, name in (('behavior', 'figure_4.png'), ('knowledge', 'figure_WA7.png')):
    print(o, curve(o, os.path.join('output', name)), 'specifications')
