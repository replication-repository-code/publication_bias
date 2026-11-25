# Figures 2 and WA5: distribution of Z statistics (estimate / standard error), Kaiser et al. (2022) versus newer trials, bins of
# 0.40: effects on financial behavior and downstream outcomes (figure 2) and on financial knowledge (figure WA5, web appendix R).
# Also prints the counts in the bin that contains 1.96 and the shares of estimates with Z >= 1.96 and Z <= -1.96.
# Run from this folder: python3 figure_2_and_WA5.py
from figure_helpers import *
for fig, f, out in [('Figure 2', 'C_full', 'output/figure_2.png'), ('Figure WA5', 'C_know', 'output/figure_WA5.png')]:
    d = pd.read_csv(f'data/{f}.csv'); z = d.yi / d.sei
    print(f'{fig}, ' + z_histogram(d, out, w=0.40))
    k = d.source.str.startswith('Kaiser')
    print(f"{fig}: {len(d)} estimates ({k.sum()} Kaiser et al., {(~k).sum()} newer); "
          f"Z >= 1.96: {(z >= 1.96).mean():.1%} (Kaiser {(z[k] >= 1.96).mean():.1%}, newer {(z[~k] >= 1.96).mean():.1%}); "
          f"Z <= -1.96: {(z <= -1.96).mean():.1%} (Kaiser {(z[k] <= -1.96).mean():.1%}, newer {(z[~k] <= -1.96).mean():.1%})")
