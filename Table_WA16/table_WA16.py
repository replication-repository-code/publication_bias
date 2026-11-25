"""Table WA16: costs and dollar benefits of financial education programs, as reported and corrected for publication bias (web
appendix P).
Programs: every program that reports both its cost per participant and an effect on behavior in dollars, or one we can convert
to dollars, 11 of the 20 studies in Kaiser et al. (2022) that report costs and 4 of the newer trials; one row per program arm. Per
participant, in U.S. dollars: for studies in Kaiser et al., costs as reported in the study; for the newer trials, costs and effects
as reported or converted at the exchange rate of the year of the intervention. The dollar effects, their standard errors, and the
costs are entered below (P), with, for each effect, whether we count it in the program's benefit (rules in web appendix P and the
note to table WA16). Gine and Mansuri (2014) report savings in logs: we convert the effect to dollars as they convert their effect
on expenditure, from the control mean of the log, exp(mean) x (exp(b) - 1), with its standard error by the delta method, at their
exchange rate of 60 rupees per dollar; their cost is the recurring cost of the training, which they use in their benefit-cost
analysis. Three more studies in Kaiser et al. report costs and dollar effects only on outcomes we do not count as benefits, the
loan amount (Bruhn and Zia 2013) and the amount remitted (Gibson, McKenzie, and Zia 2014, two samples); the last line of the log
gives the counts and medians if they were counted.
Benefit: the sum of the counted (nonredundant) dollar effects. g: that sum divided by its standard deviation, where each effect's
standard deviation is its standard error divided by SE_g, the median standard error of g of the trial's estimates in our data
(data/C_full.csv), treating the effects of a program as uncorrelated; for the two effects reported without a standard error
(Bhutoria and Vignoles 2018; Collins and Urban 2016), the standard deviation is the effect divided by the effect size of that
outcome in our data. Corrected g and benefit-cost ratio: divided by six (an uncorrected average effect of about g = 0.12 against a
corrected one of about 0.02). Breakeven g: the effect at which the benefit would just equal the cost, g x cost / benefit. The log
also reports the counts and medians when benefits are divided by five or by four instead of six.
Run from this folder: python3 table_WA16.py   (writes output/table_WA16.csv, formatted as in the paper, and
output/table_WA16_programs.csv, with full precision and the effects counted)"""
import math
import numpy as np, pandas as pd


def from_log(b, se, mean_log, rate):
    """Dollar effect of a coefficient b (standard error se) on the log of an outcome whose control mean of the log is mean_log,
    at rate units of local currency per dollar, as Gine and Mansuri (2014) convert their effect on expenditure."""
    level = math.exp(mean_log)
    return level * (math.exp(b) - 1) / rate, level * math.exp(b) * se / rate


c = pd.read_csv('data/C_full.csv')
seg = lambda key, col='cluster': float(c[c[col] == key].sei.median())
# trial (cluster) or paper identifier of each program in the combined data
SEG = {'Attanasio': seg('K6'), 'Barua': seg('K8'), 'Berry': seg('K15'), 'Bhutoria': seg('K17'), 'Bruhn': seg('P017', 'paper'),
       'Calderone': seg('K28'), 'Collins2013': seg('K38'), 'Collins': seg('K39'), 'Drexler': seg('K42'), 'Gine': seg('K53'),
       'Grohmann': seg('P044', 'paper'), 'BruhnZia': seg('K25'), 'Gibson1': seg('K50'), 'Gibson2': seg('K51'),
       'Horn': seg('P047', 'paper'), 'Sayinzoga': seg('K70'), 'Supanantaroek': seg('K76'), 'Theodos': seg('P083', 'paper')}
# Gine and Mansuri (2014), table A5: log of savings, business training 0.140 (0.153), control mean of the log 9.42 (60 Rs per USD)
GM_SAVINGS, GM_SE = from_log(0.140, 0.153, 9.42, 60.0)
# (study, arm label, cost, [(outcome, dollar effect, SE or None, coded g if no SE, counted?)])
P = [
 ('Attanasio et al. (2019)', 'Attanasio', '', 22.00, [('Account balance at partner banks', 1.41, 0.63, None, True)]),
 ('Barua et al. (2020)', 'Barua', '', 43.50, [('Savings held in the Philippines', 258.15, 1813.01, None, True),
                                               ('Savings held in Singapore', -52.12, 59.24, None, True)]),
 ('Berry et al. (2018)', 'Berry', '1', 0.62, [('Total savings', -0.23, 1.05, None, True)]),
 ('Berry et al. (2018)', 'Berry', '2', 0.62, [('Total savings', -0.88, 1.03, None, True)]),
 ('Bhutoria and Vignoles (2018)', 'Bhutoria', '', 0.76, [('Total savings', 1.52, None, 0.289963, True)]),
 ('Bruhn et al. (2025)', 'Bruhn', '', 74.57, [('Interest saved on credit card debt, per year', 12.75, 5.76, None, True),
                                              ('Interest saved on overdrafts, per year', 3.83, 1.62, None, True)]),
 ('Calderone et al. (2018)', 'Calderone', '1', 25.00, [('Savings account balance at partner company', 0.97, 0.46, None, False),
      ('Other formal savings', 38.54, 13.32, None, False), ('Informal savings', -0.33, 1.02, None, False),
      ('Total savings', 39.00, 13.46, None, False), ('Total savings including assets', 69.08, 20.08, None, True),
      ('Total amount of outstanding loans', 8.46, 23.45, None, False), ('Spending on temptation goods', -0.37, 0.18, None, False)]),
 ('Calderone et al. (2018)', 'Calderone', '2', 25.00, [('Savings account balance at partner company', 1.80, 0.71, None, False),
      ('Other formal savings', 16.88, 11.97, None, False), ('Informal savings', 0.28, 1.14, None, False),
      ('Total savings', 19.65, 12.20, None, False), ('Total savings including assets', 56.87, 22.01, None, True),
      ('Total amount of outstanding loans', -23.48, 25.27, None, False), ('Spending on temptation goods', -0.43, 0.19, None, False)]),
 # Collins (2013), table 4 (ITT): savings (bank account balances); debt excludes mortgages; net worth = savings - debt, which
 # leaves out the assets the debt may finance (cars), so we count savings
 ('Collins (2013)', 'Collins2013', '', 100.00, [('Savings', -353.96, 363.63, None, True),
      ('Total outstanding debt', 4596.13, 1966.05, None, False), ('Net worth (savings minus debt)', -4361.00, 1906.91, None, False)]),
 ('Collins and Urban (2016)', 'Collins', '', 210.00, [('Retirement contributions, per year (12 x $26 per month)', 312.00, None, 0.783564, True)]),
 ('Drexler et al. (2014)', 'Drexler', '1', 19.60, [('Total savings', 33.24, 66.72, None, True), ('Savings in the past month', 12.40, 16.24, None, False)]),
 ('Drexler et al. (2014)', 'Drexler', '2', 19.60, [('Total savings', -82.66, 81.45, None, True), ('Savings in the past month', 8.43, 18.09, None, False)]),
 ('Gine and Mansuri (2014)', 'Gine', '', 20.44, [('Savings, converted from logs', GM_SAVINGS, GM_SE, None, True)]),
 ('Grohmann, Menkhoff, and Seitz (2022)', 'Grohmann', '1', 29.70, [('Investment in physical capital', -2.89, 25.27, None, True),
      ('Inventory', 245.33, 254.92, None, True), ('Formal savings', 75.64, 68.00, None, True), ('Informal savings', 7.42, 14.90, None, True),
      ('Business loan amount', 17650.08, 57070.46, None, False)]),
 ('Grohmann, Menkhoff, and Seitz (2022)', 'Grohmann', '2', 35.60, [('Investment in physical capital', 20.48, 27.91, None, True),
      ('Inventory', 254.13, 233.03, None, True), ('Formal savings', 71.22, 57.33, None, True), ('Informal savings', 43.90, 15.31, None, True),
      ('Business loan amount', -2842.02, 55400.33, None, False)]),
 ('Horn et al. (2023)', 'Horn', '1', 29.00, [('Total savings', 33.45, 26.58, None, True), ('Business investment in the past year', 10.09, 24.65, None, True)]),
 ('Horn et al. (2023)', 'Horn', '2', 63.00, [('Total savings', 41.58, 30.67, None, True), ('Business investment in the past year', 54.78, 24.04, None, True)]),
 ('Horn et al. (2023)', 'Horn', '3', 92.00, [('Total savings', 63.40, 28.33, None, True), ('Business investment in the past year', 28.20, 28.23, None, True)]),
 ('Sayinzoga et al. (2016)', 'Sayinzoga', '', 3.50, [('Savings in cooperative, past 12 months', 25.56, 8.12, None, True)]),
 ('Supanantaroek et al. (2016)', 'Supanantaroek', '', 8.00, [('Total savings', 6.07, 7.18, None, True)]),
 ('Theodos et al. (2020)', 'Theodos', '1', 0.53, [('Reduction in revolving debt at partner credit union', 104.20, 51.74, None, False),
      ('Reduction in total revolving debt', 144.50, 154.90, None, True), ('Purchases on partner credit union card', 1.88, 8.74, None, False),
      ('Total credit for open revolving trades', -595.40, 198.60, None, False)]),
 ('Theodos et al. (2020)', 'Theodos', '2', 0.53, [('Reduction in revolving debt at partner credit union', 57.58, 51.11, None, False),
      ('Reduction in total revolving debt', 74.83, 154.80, None, True), ('Purchases on partner credit union card', 2.30, 8.74, None, False),
      ('Total credit for open revolving trades', -679.30, 199.70, None, False)]),
]
# not in the table: dollar effects only on outcomes we do not count as benefits (used only in the last line of the log)
NOT_BENEFITS = [
 ('Bruhn and Zia (2013)', 'BruhnZia', '', 245.00, [('Loan amount (0.429 KM, SE 849.06, at 0.70 USD per KM)', 0.30, 594.34, None, True)]),
 ('Gibson, McKenzie, and Zia (2014), sample 1', 'Gibson1', '', 22.90, [('Total amount remitted', 3.06, 26.42, None, True)]),
 ('Gibson, McKenzie, and Zia (2014), sample 2', 'Gibson2', '', 22.90, [('Total amount remitted', 62.94, 243.66, None, True)]),
]


def table(programs):
    rows = []
    for study, key, arm, cost, outs in programs:
        sg = SEG[key]; B = 0.0; var = 0.0
        for name, b, se, g0, counted in outs:
            if not counted: continue
            sd = se / sg if se is not None else b / g0
            B += b; var += sd ** 2
        g = B / math.sqrt(var)
        rows.append(dict(study=study, arm=arm, cost=cost, benefit=B, g=g, g_corrected=g / 6, benefit_cost=B / cost,
                         benefit_cost_corrected=B / cost / 6, breakeven_g=(g * cost / B) if B > 0 else np.nan, se_g=sg,
                         effects_counted='; '.join(o[0] for o in outs if o[4]),
                         effects_not_counted='; '.join(o[0] for o in outs if not o[4])))
    return pd.DataFrame(rows)


R = table(P)
R.to_csv('output/table_WA16_programs.csv', index=False)


def num(v, d):
    s = f'{v:.{d}f}'
    return '−' + s[1:] if s.startswith('-') else s


def usd(v):
    return ('−$' if v < 0 else '$') + f'{abs(v):,.2f}'


def breakeven(v):
    return '—' if np.isnan(v) else '< 0.001' if v < 0.001 else num(v, 3)


T = pd.DataFrame(dict(study=R.study, arm=R.arm, g_reported=R.g.map(lambda v: num(v, 3)),
                      g_corrected=R.g_corrected.map(lambda v: num(v, 3)), cost=R.cost.map(usd), benefit=R.benefit.map(usd),
                      benefit_cost_reported=R.benefit_cost.map(lambda v: num(v, 2)),
                      benefit_cost_corrected=R.benefit_cost_corrected.map(lambda v: num(v, 2)), breakeven_g=R.breakeven_g.map(breakeven)))
T.to_csv('output/table_WA16.csv', index=False)
pd.set_option('display.width', 250)
print(T.to_string(index=False))
print(f'\nGine and Mansuri (2014): savings effect ${GM_SAVINGS:.2f} (SE ${GM_SE:.2f}), converted from logs.')
best = R.loc[R.groupby('study').benefit_cost.idxmax()]
print(f"\n{len(best)} programs, {len(R)} arms. Programs with at least one arm whose benefit exceeds its cost: {int((best.benefit_cost > 1).sum())}; "
      f"more than six times its cost: {int((best.benefit_cost > 6).sum())} ({'; '.join(best[best.benefit_cost > 6].study)}).")
print(f'Median benefit-cost ratio across the {len(R)} arms: {R.benefit_cost.median():.2f} reported, {R.benefit_cost_corrected.median():.2f} corrected.')
for k in (5, 4):
    print(f"Dividing by {k} instead of 6: programs whose benefit exceeds {k} times its cost: {int((best.benefit_cost > k).sum())} "
          f"({'; '.join(best[best.benefit_cost > k].study)}); median corrected ratio {R.benefit_cost.median() / k:.2f}.")
A = table(P + NOT_BENEFITS); ba = A.loc[A.groupby('study').benefit_cost.idxmax()]
print(f"Counting also the loan amount (Bruhn and Zia 2013) and the amounts remitted (Gibson, McKenzie, and Zia 2014): {len(ba)} programs, "
      f"{len(A)} arms; benefit exceeds cost: {int((ba.benefit_cost > 1).sum())}; more than six times: {int((ba.benefit_cost > 6).sum())}; "
      f"median {A.benefit_cost.median():.2f} reported, {A.benefit_cost_corrected.median():.2f} corrected.")
