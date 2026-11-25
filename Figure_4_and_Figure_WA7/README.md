# Figures 4 and WA7: specification curves of the corrected effects

Every corrected estimate of the average effect from the corrections used in the paper, for the combined data and for Kaiser et
al.'s data alone: 59 specifications for financial behavior and downstream outcomes (figure 4, body of the paper and web appendix
N) and 59 for financial knowledge (figure WA7, web appendix R). The upper panel of each figure ranks the estimates; the lower panel
shows the data and the correction behind each one (Simonsohn, Simmons, and Nelson 2020).

## Specifications (per outcome)

| Choice | Options |
|---|---|
| Data, behavior | Combined data (`C_full`: 1,092 estimates, 130 trials); Kaiser et al.'s data alone (`K_only`: 454 estimates, 64 studies) |
| Data, knowledge | Combined data (`C_know`: 449 estimates, 115 trials); Kaiser et al.'s data alone (`K_know`: 214 estimates, 50 studies) |
| Andrews and Kasy (2019), maximum likelihood | Equation 1 or equation 2, cutoffs at 1.645, 1.96, or both, × *t*-distributed or normal true effects: 12; standard errors clustered by trial |
| Andrews and Kasy (2019), GMM (nonparametric) | Equation 1 or equation 2, cutoffs at 1.645, 1.96, or both: 6; for knowledge, equation 2 in its three-region version, as in table WA20 |
| 3PSM and 4PSM | One step at one-tailed *p* = .025 (3PSM) or steps at .025 and .50 (4PSM), *t*-distributed true effects (Andrews and Kasy's *t* specification with these steps): 2; standard errors clustered by trial |
| Precision-based corrections | PET, PEESE, UWLS, WAAP: weights 1/SE², CR2 cluster-robust standard errors by trial (Satterthwaite degrees of freedom): 4 |
| Funnel-plot imputation | Trim-and-fill (L0) on the two-level random-effects model: 1 |
| Bayesian model averaging | RoBMA-PSMA, multilevel and on trial averages: 2; for the combined data, also the two alternative priors of table WA13 (trial averages): 2 |
| Each trial counted once (combined data) | Equation 2 with *t*-distributed true effects, each estimate's contribution to the likelihood weighted by one over the number of estimates from its trial (web appendix J): 3 |

This gives 32 specifications for the combined data and 27 for Kaiser et al.'s data alone. PET-PEESE is not counted separately,
because it equals PET or PEESE; the `note` column of `output/spec_curve_estimates.csv` says which one it selects. The uncorrected
random-effects estimates and the Andrews and Kasy models without selection are in that table as benchmarks (`corrected = False`)
but not in the curves. The intervals are 95% confidence intervals, except for RoBMA-PSMA (95% credible intervals).

## Sources of the estimates

| Estimates | Source |
|---|---|
| Uncorrected random effects; Andrews and Kasy *t* models (equations 1 and 2); normal models, 3PSM, and 4PSM for the combined data; GMM for `C_full`, `K_only`, and `C_know`; RoBMA-PSMA; equation 2 with each trial counted once | Outputs of other folders, copied to `data/`: `random_effects_behavior.csv` and `random_effects_knowledge.csv` (`Figure_3`, `Figure_WA6`), `andrews_kasy_t_and_normal_behavior.csv` (`Table_WA8`), `andrews_kasy_t_and_normal_knowledge.csv` (`Table_WA18`), `selection_models.csv` (`Table_WA10`), `gmm.csv` (`Table_WA9`), `robma.csv` (`Table_WA13`), and `random_effects_other_samples.csv`, `andrews_kasy_t_other_samples.csv`, `robma_summaries.csv`, `weighted_eq2_behavior.csv`, and `weighted_eq2_knowledge.csv` (`Results_reported_in_text`) |
| Normal models, 3PSM, and 4PSM for Kaiser et al.'s data alone | `1_andrews_kasy.R` → `output/ak_fits.csv` |
| GMM for `K_know` | `2_gmm.R` → `output/gmm_fits.csv` |
| PET, PEESE, UWLS, WAAP, trim-and-fill | `3_corrections.py` → `output/precision_and_trimfill.csv` |

`1_andrews_kasy.R` and `2_gmm.R` use the package's Andrews and Kasy and GMM code unchanged (`ak_functions.R`, `ak_likelihood_KH.R`,
`ak_t.R`, `gmm_functions.R`, `gmm_PublicationbiasPackage.R`). `1_andrews_kasy.R` also refits the normal models, the 3PSM, and the
4PSM for the combined data, and `4_assemble.py` checks that these fits reproduce the package's (to 5 × 10⁻⁵). `3_corrections.py`
implements the estimators of `Table_WA10/1_corrections.R` without R packages (`corrections_functions.py`) and checks that it
reproduces UWLS, WAAP, PET, PEESE, and trim-and-fill in `Table_WA10/output/corrections.csv` (to six decimals). RoBMA-PSMA is not
refitted here. The GMM estimates (`gmm.csv`) and the estimates that count each trial once enter rounded to three decimals, as
reported in table WA9 and in the text of web appendix J.

## Run

From this folder, in order:

    Rscript 1_andrews_kasy.R normal      # under a minute
    Rscript 1_andrews_kasy.R t           # about 10 minutes
    Rscript 2_gmm.R                      # about a minute
    python3 3_corrections.py             # seconds
    python3 4_assemble.py                # writes output/spec_curve_estimates.csv
    python3 5_figures.py                 # writes output/figure_4.png and output/figure_WA7.png
    python3 6_summaries.py               # writes the summaries

`1_andrews_kasy.R` and `2_gmm.R` append each fit to their output file as it finishes and skip the fits already there, so an
interrupted run can be resumed; delete `output/ak_fits.csv`, `output/gmm_fits.csv`, and `output/cache/` to refit everything.
Steps 3 to 6 can be run alone, from the files provided.

## Outputs

- `output/figure_4.png`, `output/figure_WA7.png`: the figures.
- `output/spec_curve_estimates.csv`: one row per specification (outcome, data, correction, variant, selection function, true
  effects, estimate, standard error, interval, source).
- `output/spec_curve_summary.csv`: median, quartiles, range, and shares of specifications whose interval excludes zero, whose
  estimate is below the uncorrected estimate of the same data (or below half of it), and whose whole interval lies below that
  uncorrected estimate or below 0.10; for all corrected specifications and by data and correction (the numbers in the body of the
  paper and in web appendix N).
- `output/spec_curve_paper_estimates_rank.csv`: the rank of the estimates of figure 3 (equation 2, *t*-distributed true effects,
  combined data) among the 59 corrected specifications.
- `output/spec_curve_variance_decomposition.csv`: share of the variance of the corrected estimates explained by the choice of
  data and by the choice of correction, and, for the Andrews and Kasy maximum-likelihood models, by the data, the selection
  function, and the distribution of true effects (Shapley decomposition of the R² of a regression on indicators).
- `output/log_*.txt`: the console output of each script.
