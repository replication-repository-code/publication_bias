# Table WA1: estimates (trials) by outcome class and source.  Run: python3 table_WA1.py
import pandas as pd
b = pd.read_csv('data/C_full.csv'); k = pd.read_csv('data/C_know.csv')
def cnt(d): return f"{len(d)} ({d.cluster.nunique()})"
kb, ub = b[b.source.str.startswith('Kaiser')], b[~b.source.str.startswith('Kaiser')]
kk, uk = k[k.source.str.startswith('Kaiser')], k[~k.source.str.startswith('Kaiser')]
rows = [('A: Financial behavior', cnt(kb), cnt(ub[ub.role.str[:2] == 'R3'])),
        ('A: Decision quality', '-', cnt(ub[ub.role.str[:2] == 'DQ'])),
        ('A: Financial well-being', '-', cnt(ub[ub.role.str[:2] == 'R4'])),
        ('A: Economic outcomes', '-', cnt(ub[ub.role.str[:2] == 'R5'])),
        ('A: Mechanisms', '-', cnt(ub[ub.role.str[:2] == 'R6'])),
        ('A: Total', cnt(kb), cnt(ub) + f", of which {ub.cluster.str.startswith('U-').groupby(ub.cluster).first().sum()} not in Kaiser et al."),
        ('B: Objective financial knowledge', cnt(kk), cnt(uk[uk.category == 'A'])),
        ('B: Subjective financial knowledge', '-', cnt(uk[uk.category == 'sk'])),
        ('B: Total', cnt(kk), cnt(uk))]
out = pd.DataFrame(rows, columns=['outcome_class', 'Kaiser et al. (2022): estimates (studies)', 'Newer trials: estimates (trials)'])
out.to_csv('output/table_WA1.csv', index=False); print(out.to_string(index=False))
