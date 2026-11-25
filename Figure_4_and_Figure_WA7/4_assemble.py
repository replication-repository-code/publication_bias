"""Collect every estimate of the specification curve into one table (one row per specification).

Data: the combined data with all outcomes (behavior: C_full; knowledge: C_know) and Kaiser et al.'s data alone (K_only, K_know).
Corrections: Andrews and Kasy's equations 1 and 2 (maximum likelihood with t-distributed or normal true effects, and GMM), the 3PSM
and 4PSM with t-distributed true effects, PET, PEESE, UWLS, WAAP, trim-and-fill, and RoBMA-PSMA.
Sources: copies of outputs of other folders of the package (data/; see README) and the fits of this folder (output/ak_fits.csv,
output/gmm_fits.csv, output/precision_and_trimfill.csv). Run from this folder: python3 4_assemble.py
Writes output/spec_curve_estimates.csv."""
import os
import numpy as np, pandas as pd

OUT = 'output/spec_curve_estimates.csv'
Z = 1.959963984540054
P = lambda f: os.path.join('data', f)        # copies of the package's outputs
N = lambda f: os.path.join('output', f)      # fits of this folder

SAMPLES = {  # code: outcome, label, estimates, trials
    'C_full': ('behavior', 'Combined data, all outcomes', 1092, 130),
    'K_only': ('behavior', "Kaiser et al.'s data alone", 454, 64),
    'C_know': ('knowledge', 'Combined data, all knowledge outcomes', 449, 115),
    'K_know': ('knowledge', "Kaiser et al.'s data alone", 214, 50)}
SEL = {'sym1645': 'Eq. 1, |Z| = 1.645', 'sym196': 'Eq. 1, |Z| = 1.96', 'symboth': 'Eq. 1, |Z| = 1.645 and 1.96',
       'asym1645': 'Eq. 2, Z = −1.645, 0, 1.645', 'asym196': 'Eq. 2, Z = −1.96, 0, 1.96', 'asymboth': 'Eq. 2, Z = ±1.645, ±1.96, 0'}
AK, GMM = 'Andrews and Kasy (maximum likelihood)', 'Andrews and Kasy (GMM)'
keep = lambda df: df[df['sample'].isin(SAMPLES)]
rows = []


def add(sample, method, variant, g, se, lo, hi, source, selection='', true_effects='', corrected=True, interval='95% CI', note='',
        unit='estimates'):
    o, lab, k, t = SAMPLES[sample]
    rows.append(dict(outcome=o, sample=sample, sample_label=lab, k=k, trials=t, unit=unit, corrected=corrected, method=method,
                     variant=variant, selection=selection, true_effects=true_effects, estimate=g, se=se, ci_lo=lo, ci_hi=hi,
                     interval=interval, source=source, note=note))


# ---------------------------------------------------------------- uncorrected random-effects model (estimates nested in trials)
re = keep(pd.concat([pd.read_csv(P('random_effects_behavior.csv')), pd.read_csv(P('random_effects_knowledge.csv')),
                     pd.read_csv(P('random_effects_other_samples.csv'))]))
for x in re.itertuples():
    add(x.sample, 'Random effects, estimates nested in trials', 'uncorrected', x.mu, x.se_CR2, x.lo, x.hi,
        'package: random_effects_*.csv', corrected=False, note='CR2 cluster-robust CI (Satterthwaite df)')

# ---------------------------------------------------------------- Andrews and Kasy, maximum likelihood
pk = keep(pd.concat([pd.read_csv(P('andrews_kasy_t_and_normal_behavior.csv')), pd.read_csv(P('andrews_kasy_t_and_normal_knowledge.csv')),
                     pd.read_csv(P('andrews_kasy_t_other_samples.csv'))])).assign(src='package: andrews_kasy_*.csv')
new = pd.read_csv(N('ak_fits.csv')).assign(src='new fit (1_andrews_kasy.R)')
assert not new.duplicated(['sample', 'model', 'spec']).any() and (new.conv == 0).all()
# normal models: the package has them for C_full and C_know; the new fits must reproduce them
chk = pk[pk.model == 'normal'].merge(new[new.model == 'normal'], on=['sample', 'model', 'spec'], suffixes=('', '_new'))
assert len(chk) == 14 and np.allclose(chk.mu, chk.mu_new, atol=5e-5) and np.allclose(chk.se_mu, chk.se_mu_new, atol=5e-5), chk
ml = pd.concat([pk, new[~((new.model == 'normal') & new['sample'].isin(['C_full', 'C_know']))]])
# 3PSM and 4PSM (t): the package has them for C_full and C_know (selection_models.csv); the new fits must reproduce them
sm = keep(pd.read_csv(P('selection_models.csv')))
for x in sm[sm.true_effects.str.startswith('t')].itertuples():
    spec = {'3PSM': 'psm3', '4PSM': 'psm4'}[x.selection]
    m = (ml['sample'] == x.sample) & (ml.model == 't') & (ml.spec == spec)
    assert m.sum() == 1 and abs(ml[m].mu.iloc[0] - x.g) < 5e-5 and abs(ml[m].se_mu.iloc[0] - x.se) < 5e-5, (x.sample, spec)
    ml = pd.concat([ml[~m], pd.DataFrame([dict(sample=x.sample, model='t', spec=spec, mu=x.g, se_mu=x.se, src='package: selection_models.csv')])])
assert not ml.duplicated(['sample', 'model', 'spec']).any() and len(ml) == 4 * (7 + 7 + 2), len(ml)
for x in ml.itertuples():
    te = 't' if x.model == 't' else 'Normal'
    lo, hi = x.mu - Z * x.se_mu, x.mu + Z * x.se_mu
    if x.spec == 'none':
        add(x.sample, 'Andrews and Kasy, no selection', f'uncorrected; {te} true effects', x.mu, x.se_mu, lo, hi, x.src, true_effects=te,
            corrected=False, note='Benchmark of the selection models: same likelihood with all publication probabilities = 1')
    elif x.spec in ('psm3', 'psm4'):
        lab = {'psm3': ('3PSM', 'one step at one-tailed p = .025 (Z = 1.96)'), 'psm4': ('4PSM', 'steps at one-tailed p = .025 and .50 (Z = 1.96 and 0)')}[x.spec]
        add(x.sample, lab[0], f'{lab[1]}; t true effects', x.mu, x.se_mu, lo, hi, x.src, true_effects='t',
            note="Andrews and Kasy's t specification with the 3PSM or 4PSM steps; standard errors clustered by trial")
    else:
        add(x.sample, AK, f'{SEL[x.spec]}; {te} true effects', x.mu, x.se_mu, lo, hi, x.src, selection=x.spec, true_effects=te,
            note='Standard errors clustered by trial')
# each trial counted once: equation 2 (t) with each estimate's contribution to the likelihood weighted by one over the number of
# estimates from its trial (body and web appendix J; Results_reported_in_text/weighted_eq2.R, outputs copied to data/). The
# estimates and confidence limits enter rounded to three decimals, as reported in the text (the GMM estimates are rounded too, as
# in table WA9); standard error = width of the 95% interval / 3.92
for dom, s0 in (('behavior', 'C_full'), ('knowledge', 'C_know')):
    for x in pd.read_csv(P(f'weighted_eq2_{dom}.csv')).itertuples():
        src = f'package: weighted_eq2_{dom}.csv (rounded to 3 decimals)'
        mu, lo, hi = round(x.mu, 3), round(x.lo, 3), round(x.hi, 3); se = (hi - lo) / (2 * Z)
        if x.spec == 'none':
            add(s0, 'Andrews and Kasy, no selection', 'uncorrected; t true effects; each trial counted once', mu, se, lo, hi,
                src, true_effects='t', corrected=False, unit='trials', note='Estimates weighted by 1 / estimates in their trial')
        else:
            add(s0, AK, f'{SEL[x.spec]}; t true effects; each trial counted once', mu, se, lo, hi, src,
                selection=x.spec, true_effects='t', unit='trials',
                note='Estimates weighted by 1 / estimates in their trial; standard errors clustered by trial')

# ---------------------------------------------------------------- Andrews and Kasy, nonparametric GMM
for x in keep(pd.read_csv(P('gmm.csv'))).itertuples():
    add(x.sample, GMM, SEL[x.spec.replace('_3', '')] + (' (three regions)' if x.spec.endswith('_3') else ''), x.gmm_g, x.gmm_se,
        x.gmm_ci_lo, x.gmm_ci_hi, 'package: gmm.csv (table WA9, rounded to 3 decimals)', selection=x.spec.replace('_3', ''),
        note='Nonparametric in the distribution of true effects; standard errors clustered by trial')
g2 = pd.read_csv(N('gmm_fits.csv'))
assert not g2.duplicated(['sample', 'spec']).any() and set(g2['sample']) == {'K_know'} and len(g2) == 6
for x in g2.itertuples():
    add(x.sample, GMM, SEL[x.spec.replace('_3', '')] + (' (three regions)' if x.spec.endswith('_3') else ''), x.gmm_g, x.gmm_se,
        x.gmm_g - Z * x.gmm_se, x.gmm_g + Z * x.gmm_se, 'new fit (2_gmm.R)', selection=x.spec.replace('_3', ''),
        note='Nonparametric in the distribution of true effects; standard errors clustered by trial')

# ---------------------------------------------------------------- precision-based corrections and trim-and-fill
pt = pd.read_csv(N('precision_and_trimfill.csv'))
for x in pt.itertuples():
    if x.method in ('PET-PEESE', 'Two-level random effects (REML)') or np.isnan(x.g): continue   # PET-PEESE is PET or PEESE (note)
    note = 'CR2 cluster-robust CI (Satterthwaite df)'
    if x.method == 'WAAP':
        note = f'{int(x.n_estimates)} adequately powered estimates from {int(x.n_trials)} trials; CR2 CI'
    if x.method == 'Trim-and-fill':
        note = f'{int(x.k0)} estimates imputed ({x.side} side): equals the two-level random-effects model; CI treats estimates as independent'
    if x.method in ('PET', 'PEESE') and pt[(pt['sample'] == x.sample) & (pt.method == 'PET-PEESE')].rule.iloc[0] == x.method:
        note += '; selected by conditional PET-PEESE'
    add(x.sample, x.method, '', x.g, x.se, x.lo, x.hi, 'new computation (3_corrections.py; reproduces Table_WA10/output/corrections.csv)', note=note)

# ---------------------------------------------------------------- RoBMA-PSMA (posterior mean, 95% credible interval)
rb = pd.read_csv(P('robma_summaries.csv')); rb['sample'] = rb.model.str.split(' ').str[0]
for x in keep(rb).itertuples():
    v = x.model.split(' ')[1]
    add(x.sample, 'RoBMA-PSMA', {'ml': 'multilevel', 'avg': 'trial averages'}[v], x.mu, np.nan, x.mu_lo, x.mu_hi,
        'package: robma_summaries.csv', interval='95% CrI', note='Model-averaged posterior mean (includes the models without an effect)',
        unit='trials' if v == 'avg' else 'estimates')
pr = pd.read_csv(P('robma.csv'))
for x in pr[pr.model.str.contains('nullnormal|effprior')].itertuples():
    v = 'trial averages, near-zero null prior' if 'nullnormal' in x.model else 'trial averages, effect prior centered on a positive effect'
    add(x.model.split('_avg')[0], 'RoBMA-PSMA', v, x.mu, np.nan, x.mu_lo, x.mu_hi, 'package: robma.csv (table WA13)', interval='95% CrI',
        note='Alternative priors (table WA13)', unit='trials')

d = pd.DataFrame(rows)
d['sig_positive'] = d.ci_lo > 0
d = d.sort_values(['outcome', 'corrected', 'estimate'], ascending=[True, False, True]).reset_index(drop=True)
d.to_csv(OUT, index=False)
print(d.groupby(['outcome', 'corrected', 'method']).size().to_string())
print('total', len(d), '| corrected:', d[d.corrected].groupby('outcome').size().to_dict())

