# Figure 3, step 3: draw the figure from the outputs of steps 1 and 2.  Run: python3 3_plot.py
from figure_helpers import *
forest_by_method('output/random_effects.csv', 'output/andrews_kasy_t.csv', 'output/robma_C_full_ml_summary.csv',
                 'output/figure_3.png', *((-0.03, 0.16), [0, 0.05, 0.10, 0.15]))
