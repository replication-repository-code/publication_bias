# Replication package

**How Much Does Financial Education Change Behavior? Publication Bias and an Updated Meta-Analysis**
Daniel Fernandes, John G. Lynch, Jr., and Lena Kim

This package reproduces every table and figure of the paper and its web appendix that uses the meta-analytic data, except figure
WA1 (section 7), and the other numbers reported in the text. Folders, scripts, and output files are named after the tables and
figures of the version of the paper of 23 November 2026 (section 2 lists where to find each one). Every folder is self-contained:
it holds the data it needs (`data/`), the code (R or Python scripts and the helper files they source), and the results we
obtained (`output/`). A script run from inside its folder writes its results to that folder's `output/`; `log_*.txt` files hold
the console output of our runs.

## 1. Quick start

1. Install the software in section 4.
2. Open a terminal **in the folder of the table or figure** you want (the scripts use paths relative to the folder).
3. Run the scripts in the order of their numbers, for example in `Figure_3/`:
   ```
   Rscript 1_estimates.R
   Rscript 2_robma.R
   python3 3_plot.py
   ```
   In RStudio, set the working directory to the folder (Session > Set Working Directory > To Source File Location) and source the
   scripts.
4. Compare the new files in `output/` with the ones we provide, and with the paper.

## 2. Folders

| Folder | Paper | Scripts (run in this order) | Run time* |
|---|---|---|---|
| `00_Data_and_analysis_samples` | Data for all analyses (section 3) | `build_samples.py` | seconds |
| `Figure_2_and_Figure_WA5` | Figure 2: distribution of *Z* statistics, bins of 0.40, behavior; figure WA5: the same for knowledge (web appendix R); the shares of significant estimates in the text | `figure_2_and_WA5.py` | seconds |
| `Table_1` | Table 1: relative publication probabilities, equation 2, behavior; the row "Estimates (trials)" of each panel and the counts in its note | `table_1.R`, `table_1_counts.py` | 10–20 min; seconds |
| `Table_2` | Table 2: log-likelihood fit of equations 1 and 2, behavior | `table_2.R` | 10–30 min |
| `Table_3_and_Table_WA19` | Table 3: publication probabilities and likelihood-ratio tests of the symmetry restrictions, behavior; table WA19: the same for knowledge (web appendix R) | `table_3.R`, `table_WA19.R` | 10–30 min; 1–2 min |
| `Figure_3` | Figure 3: uncorrected, equation 1, equation 2, and multilevel RoBMA-PSMA estimates, behavior; most estimates quoted in the text | `1_estimates.R`, `2_robma.R`, `3_plot.py` | 1–3 h (RoBMA) |
| `Table_4_and_Table_WA20` | Table 4 (behavior) and table WA20 (knowledge, web appendix R): types of corrections for publication bias, what each assumes, and the corrected effects. The script collects the outputs of other folders, copied into `data/` (section 5) | `table_4_and_WA20.py` | seconds |
| `Figure_4_and_Figure_WA7` | Figure 4 (behavior) and figure WA7 (knowledge, web appendix R): specification curves of the corrected estimates; the summaries in the body and in web appendix N. See the folder's README | `1_andrews_kasy.R`, `2_gmm.R`, `3_corrections.py`, `4_assemble.py`, `5_figures.py`, `6_summaries.py` | about 15 min |
| `Table_WA1` | Table WA1: estimates by outcome class and source (web appendix A) | `table_WA1.py` | seconds |
| `Table_WA2` | Table WA2: uncorrected random-effects models (web appendix B); the models with a level for papers quoted in the text | `table_WA2.R` | seconds |
| `Table_WA3` | Table WA3: excess mass of *Z* statistics just above 1.645 and 1.96, bins of 0.05, 0.10, and 0.40, behavior and knowledge (web appendix C); the placebo thresholds and the tests without the estimates whose *Z* lies at a threshold; the excess just above 1.96 for behavior without significance-derived standard errors, over 80 specifications | `table_WA3.py`, `excess_mass_robustness.py` | seconds; under a minute |
| `Table_WA5` | Table WA5: newer trials, journal articles versus working papers and reports (web appendix E) | `1_counts.py`, `2_random_effects.R` | seconds |
| `Table_WA6` | Table WA6: Kaiser et al.'s behaviors plus the newer trials' financial behaviors (web appendix F) | `table_WA6.R` | 10–30 min |
| `Table_WA7` | Table WA7: Kaiser et al.'s behavior data alone (web appendix G) | `table_WA7.R` | 10–30 min |
| `Table_WA8` | Table WA8: models with *t*-distributed and with normal true effects (web appendix H), except its last column (next folder); also the models without selection (the *t* model without selection is quoted in the body) and the implied distribution of true effects quoted in the text | `table_WA8.R`, `implied_distribution.R` | 20–60 min |
| `Likelihood_ratio_tests_clustered` | Likelihood-ratio tests that allow for dependence among estimates from the same trial: the last column of table WA8, the notes of tables 2 and 3, the body, and web appendices G and H. See the folder's README | `1_eq2_vs_eq1.R`, `2_table_WA8.R` | about 30 min |
| `Table_WA9` | Table WA9: parametric versus nonparametric (GMM) estimates (web appendix I); the GMM estimates and 95% confidence intervals quoted in the body | `table_WA9.R` | 1–3 h |
| `Table_WA10` | Table WA10: UWLS, WAAP, PET, PEESE, PET-PEESE, 3PSM and 4PSM with *t*-distributed true effects, trim-and-fill, and leave-one-trial-out (web appendix J); the statistical power and the shares of the inverse-variance weight quoted in the text; the rows of tables 4 and WA20 for these corrections | `1_corrections.R`, `2_selection_models.R`, `3_table_WA10.py` | 15–30 min |
| `Figure_WA2` | Figure WA2: funnel plot of trial-averaged effect sizes, behavior (web appendix K); the regression tests for funnel-plot asymmetry quoted in web appendices J and K | `funnel.R` | seconds |
| `Table_WA11` | Table WA11: multilevel RoBMA-PSMA, behavior (web appendix M) | `table_WA11.R` | 1–3 h |
| `Table_WA12` | Table WA12: most probable RoBMA-PSMA models (web appendix M) | `table_WA12.R` | 1–3 h |
| `Table_WA13` | Table WA13: RoBMA-PSMA under alternative priors (web appendix M); also the RoBMA-PSMA rows of tables 4 and WA20 | `table_WA13.R` | several hours |
| `Figure_WA4_and_Table_WA14` | Figure WA4 and table WA14: subgroups, behavior (web appendix O); the subgroup estimates in the body | `1_subgroups.R`, `2_plot.py` | 30–60 min |
| `Table_WA15` | Table WA15: subgroup estimates with one selection function per dimension, behavior (web appendix O); the same for knowledge, quoted in web appendix R (both in `output/table_WA15.csv`, column `domain`) | `table_WA15.R` | 30–60 min |
| `Table_WA16` | Table WA16: costs and dollar benefits of financial education programs, as reported and corrected for publication bias (web appendix P); the counts and medians in the body and in web appendix P | `table_WA16.py` | seconds |
| `Web_appendix_Q_Fernandes_2014` | Web appendix Q: Andrews and Kasy corrections with normal true effects for the 15 randomized trials in Fernandes, Lynch, and Netemeyer (2014): equation 2 with regions split at *Z* = 0 and 1.645, equation 1 with a cutoff at 1.645, the uncorrected estimate, and the log likelihoods quoted in the text. The data (`data/fernandes_2014_rcts.csv`) give each trial's effect as a correlation *r* with its standard error; the script converts the estimates to *g* = 2*r*/√(1 − *r*²) | `fernandes_2014.R` | seconds |
| `Table_WA17` | Table WA17: relative publication probabilities, knowledge (web appendix R) | `table_WA17.R` | 10–20 min |
| `Table_WA18` | Table WA18: fit of equations 1 and 2, knowledge (web appendix R); the normal models of table WA20 | `table_WA18.R` | 20–40 min |
| `Figure_WA6` | Figure WA6: estimates by method, knowledge (web appendix R) | `1_estimates.R`, `2_robma.R`, `3_plot.py` | 1 h (RoBMA) |
| `Web_appendix_R_knowledge_subgroups` | Subgroup estimates for knowledge, each subgroup with its own publication probabilities, quoted in web appendix R | `1_subgroups.R` | 10–30 min |
| `RoBMA_single_level` | Single-level RoBMA-PSMA without clustering by trial, Kaiser et al.'s data and combined data (response to Reviewer 3 in the revision notes; not in the paper) | `robma_single_level.R` | 30–60 min |
| `Results_reported_in_text` | Numbers quoted in the text but not shown in a table or figure (section 6) | `1_andrews_kasy_and_random_effects.R`, `2_robma.R`, `3_gmm.R`, `loo_eq2.R`, `weighted_eq2.R` | many hours |

\*On one core of a laptop. The Andrews and Kasy models are slow because each is estimated from several starting values and its
likelihood integrates over the distribution of true effects. RoBMA-PSMA runs three MCMC chains in parallel.

**Where to find each table and figure** (output files, in the folders above):

| Paper | Output |
|---|---|
| Figure 2, figure WA5 | `Figure_2_and_Figure_WA5/output/figure_2.png`, `figure_WA5.png` |
| Table 1 | `Table_1/output/table_1.csv`, `table_1_counts.csv` |
| Table 2 | `Table_2/output/table_2.csv` |
| Table 3, table WA19 | `Table_3_and_Table_WA19/output/table_3.csv`, `table_WA19.csv` |
| Figure 3 | `Figure_3/output/figure_3.png` |
| Table 4, table WA20 | `Table_4_and_Table_WA20/output/table_4.csv`, `table_WA20.csv` |
| Figure 4, figure WA7 | `Figure_4_and_Figure_WA7/output/figure_4.png`, `figure_WA7.png` |
| Table WA1 | `Table_WA1/output/table_WA1.csv` |
| Table WA2 | `Table_WA2/output/table_WA2.csv` |
| Table WA3 | `Table_WA3/output/table_WA3.csv` |
| Table WA5 | `Table_WA5/output/table_WA5_counts.csv`, `table_WA5_random_effects.csv` |
| Table WA6 | `Table_WA6/output/table_WA6.csv` |
| Table WA7 | `Table_WA7/output/table_WA7.csv` |
| Table WA8 | `Table_WA8/output/table_WA8.csv`; last column: `Likelihood_ratio_tests_clustered/output/robust_lr_vs_best_C_full.csv` |
| Table WA9 | `Table_WA9/output/table_WA9.csv` |
| Table WA10 | `Table_WA10/output/table_WA10.csv` |
| Figure WA2 | `Figure_WA2/output/figure_WA2.png` |
| Table WA11 | `Table_WA11/output/table_WA11_summary.csv`, `table_WA11_models.csv`, and the predictive density in `log_table_WA11.R.txt` |
| Table WA12 | `Table_WA12/output/table_WA12_ml_models.csv`, `table_WA12_avg_models.csv` |
| Table WA13 | `Table_WA13/output/table_WA13.csv` |
| Figure WA4, table WA14 | `Figure_WA4_and_Table_WA14/output/figure_WA4.png`, `table_WA14.csv` |
| Table WA15 | `Table_WA15/output/table_WA15.csv` |
| Table WA16 | `Table_WA16/output/table_WA16.csv` |
| Table WA17 | `Table_WA17/output/table_WA17.csv` |
| Table WA18 | `Table_WA18/output/table_WA18.csv` |
| Figure WA6 | `Figure_WA6/output/figure_WA6.png` |

Figure 1, figure WA1, figure WA3, and table WA4 are not in this package (section 7).

**Large files.** The RoBMA folders do not include the fitted RoBMA objects (`.rds`, up to 180 MB each). Their `output/` folders
contain the summaries the paper uses (`*_summary.csv`, `*_models.csv`, `table_WA13.csv`). The plotting scripts read these
summaries, so the figures can be redrawn without refitting. Running a RoBMA script refits the models and saves the `.rds` files;
if an `.rds` file is already in `output/`, the script reuses it.

## 3. Data (`00_Data_and_analysis_samples`)

- `FinEd_meta_update_effect_sizes_OSF.xlsx`: our coding of the newer randomized trials of financial education, for the 873
  estimates that enter the analyses (98 papers, 88 trials). The sheet *Effect sizes* has one row per estimate (treatment arm,
  outcome, and time of measurement) and the variables the analyses use: `Row ID` (estimate), `Paper ID`, `Trial ID`, `Authors`
  and `Year` of the paper, `Publication status` (as coded), `Top journal (your list)` (Yes = top economics journal), `Intervention
  type`, `Low-income sample`, `Outcome category (A/sk/B/C/D/E/F/Other)` (A objective knowledge, sk subjective knowledge, B–F
  Kaiser et al.'s behavior codes, Other outcomes outside them), `Outcome role (decision quality / positive outcomes)` (R1 to R6
  and DQ, below), `Hedges' g aligned incl. authors' view` (Hedges' *g*, positive = beneficial), and `SE(g)` (its standard
  error).
- `Kaiser_et_al_2022_data_OSF.xlsx`: the effect sizes of Kaiser, Lusardi, Menkhoff, and Urban (2022, *Journal of Financial
  Economics*), from their replication files: all 673 estimates (458 behavior, codes B–G, and 215 knowledge, code A), with the
  variables the analyses use: `Study_ID`, `Study_Name`, `outcome_code`, `outcome_domain`, `hedges_g`, `se_g`,
  `intervention_channel`, `low_income_dummy`, and `top_econ_journal`.
- `build_samples.py` builds `analysis_samples/*.csv` from these two files. It assigns each estimate to the samples below, drops
  four behavior estimates and one knowledge estimate that Kaiser et al.'s file lists twice (same study, estimate, and standard
  error), assigns each follow-up paper to the Kaiser et al. trial it follows up, checks that every trial has a single low-income
  code, and groups Kaiser et al.'s two G1 estimates (their "Saving & Investing" domain) with saving and investing. It also
  writes `kaiser_raw.csv`, Kaiser et al.'s 458 behavior estimates before duplicates are removed, for the folder
  `RoBMA_single_level`. The CSV files in every other folder are copies of these files.

| File | Content | Estimates (trials) |
|---|---|---|
| `C_full.csv` | Kaiser et al.'s behavior estimates (codes B–G) + the newer trials' financial behavior, decision quality, financial well-being, economic outcomes, and mechanisms | 1,092 (130) |
| `C_noR6.csv`, `C_noR5R6.csv` | As C_full without mechanisms; without mechanisms and economic outcomes | 1,018 (128); 982 (128) |
| `C_R3.csv` | Kaiser et al.'s behavior estimates + the newer trials' financial behaviors only | 909 (125) |
| `K_only.csv` | Kaiser et al.'s behavior estimates alone | 454 (64) |
| `C_know.csv` | Kaiser et al.'s knowledge estimates (code A) + the newer trials' objective (A) and subjective (sk) knowledge | 449 (115) |
| `C_knowobj.csv` | As C_know with objective knowledge only | 364 (105) |
| `K_know.csv` | Kaiser et al.'s knowledge estimates alone | 214 (50) |

`paper_status.csv` lists the 98 papers of the update that contribute estimates to the analyses, with their authors, year,
publication status as coded, and `published` (1 = journal article or book chapter by October 2026, 0 = working paper, report, or
thesis; two working papers that have since appeared in journals count as published, see `build_samples.py`): 73 are published and
25 are not. It is used for table WA5.

Variables in the CSV files: `source` (Kaiser et al. or the update), `study`, `paper`, `cluster` (trial; used for clustering and as
the random-effects level), `yi` (Hedges' *g*, positive = beneficial), `sei` (standard error), `category` (outcome code), `role` (R1
objective knowledge, R2 subjective knowledge, R3 financial behavior, DQ decision quality, R4 financial well-being, R5 economic
outcomes, R6 mechanisms), and the subgroup codes `interv`, `lowinc`, `topj`, `otype`.

## 4. Software

We ran the analyses with:

- R 4.3.3 with metafor 4.4.0, clubSandwich 0.7.0, RoBMA 4.0.0 (requires JAGS 4.3; we used JAGS 4.3.2), MASS 7.3-60, loo 2.6.0,
  and coda 0.19-4. The scripts that use only the Andrews and Kasy and GMM code need only base R and MASS: those of `Table_1`,
  `Table_2`, `Table_3_and_Table_WA19`, `Table_WA6` to `Table_WA9`, `Table_WA15`, `Table_WA17`, `Table_WA18`,
  `Likelihood_ratio_tests_clustered`, and `Web_appendix_Q_Fernandes_2014`, the two R scripts of `Figure_4_and_Figure_WA7`, and
  `3_gmm.R`, `loo_eq2.R`, and `weighted_eq2.R` in `Results_reported_in_text`.
- Python 3 with pandas, numpy, matplotlib (3.11), openpyxl, and scipy.

Install the R packages with `install.packages(c("metafor", "clubSandwich", "RoBMA", "MASS", "loo", "coda"))` after installing
JAGS (https://mcmc-jags.sourceforge.io). Install the Python packages with `pip install pandas numpy matplotlib openpyxl scipy`.

## 5. Methods implemented in the code

- **Uncorrected random-effects model** (`re_functions.R`): estimates nested in trials (metafor `rma.mv`, REML), with CR2
  cluster-robust standard errors and Satterthwaite degrees of freedom by trial (clubSandwich).
- **Andrews and Kasy (2019) maximum likelihood, *t*-distributed true effects** (`ak_t.R`, used through `ak_functions.R` and
  `ak_general.R`; body of the paper and every web appendix): Andrews and Kasy's own specification, as in their minimum-wage
  application. True effects follow a *t* distribution, θ = μ + τ·*T* with *T* ~ *t*(ν), and each estimate is normal around its
  true effect, *X* | θ ~ N(θ, SE²). Publication probabilities are step functions of *Z* = *g*/SE: equation 1 is symmetric (regions
  of |*Z*|), equation 2 is asymmetric (regions of signed *Z*), with cutoffs at 1.645, at 1.96, or at both; the publication
  probability of significant positive estimates is fixed at 1. The *t* is written as a scale mixture of normals, θ | *V* ~ N(μ,
  τ²/*V*) with *V* ~ Gamma(ν/2, rate ν/2), so the latent density of each estimate and the probability of each region of *Z* are
  integrals over *V*; they are computed with the trapezoid rule in log *V* (relative error below 10^-8 against adaptive
  quadrature). The log-likelihood has an analytic gradient. Each model is maximized with `nlminb` from several starting values,
  including the estimates of the normal model, with ν between 0.05 and 1,000, and restarted from its solution with `nlminb`'s
  default tolerances, which confirms the maximum and records the convergence code. Standard errors are clustered by trial
  (sandwich estimator). The three-category version of equation 2 (one probability for all negative estimates) is used for
  subgroups with fewer than five significantly negative estimates and for knowledge in table WA9 and in the GMM row of table WA20.
- **Normal true effects** (`ak_functions.R`, `ak_likelihood_KH.R`; table WA8, web appendix H; web appendix Q; and the normal rows
  of tables 4 and WA20): true effects are normal with standard deviation τ, so the latent distribution of each estimate is N(μ,
  SE² + τ²); the *t* model converges to it as ν grows. The likelihood and the robust variance come from `ak_likelihood_KH.R`, the
  code used for the earlier version of the paper (`AndrewsKasyFunctions_KH.R`) without its plotting libraries. The scripts report
  likelihood-ratio tests with *p*-values that treat estimates as independent: chi-square for restrictions on the publication
  probabilities, and a 50:50 mixture of chi-square(0) and chi-square(1) for *t* versus normal true effects (the normal is the
  boundary case of the *t*).
- **Likelihood-ratio tests that allow for dependence within trials** (`Likelihood_ratio_tests_clustered`; notes of tables 2 and
  3, last column of table WA8, body, and web appendices G and H): for nested models, the likelihood-ratio statistic is compared with
  a weighted sum of chi-square(1) variables whose weights are the eigenvalues of the ratio of the cluster-robust to the model-based
  covariance matrix of the restricted parameters (Vuong 1989; Rotnitzky and Jewell 1990). For *t* versus normal true effects, each
  statistic is divided by the design effect of the degrees-of-freedom parameter and compared with the 50:50 mixture.
- **General Andrews and Kasy model** (`ak_general.R`; tables 3, WA19, and WA15): the same likelihood with any pattern of equality
  restrictions on the publication probabilities (for example *p*1 = 1, or *p*2 = *p*3, tested against equation 2 with
  likelihood-ratio tests), and with a design matrix for the location of the distribution of true effects. With subgroup
  indicators, it estimates subgroup-specific means under one common set of publication probabilities, Andrews and Kasy's (2019)
  "conditioning on covariates."
- **Three- and four-parameter selection models** (`Table_WA10`): the 3PSM has one step at a one-tailed *p* = .025 (*Z* = 1.96),
  so it separates significant positive estimates from all others; the 4PSM adds a step at *p* = .50 (*Z* = 0). With
  *t*-distributed true effects (`2_selection_models.R`; tables 4, WA10, and WA20, and the text), they are Andrews and Kasy's *t*
  model with these steps, with standard errors clustered by trial. With normal true effects, as they are usually estimated
  (`1_corrections.R`; the 3PSM and 4PSM rows of `output/corrections.csv`, which the paper does not report), they are metafor's
  `selmodel` (step function, `alternative = "greater"`) on the random-effects model, which treats estimates as independent.
  `output/selection_models.csv` reports all eight fits (two models, normal or *t* true effects, behavior or knowledge) with their
  log-likelihoods, which are comparable across specifications.
- **Clustered Vuong test** (`Table_WA8/output/table_WA8.csv`, columns `vuong_z` and `vuong_p_one_sided`; not reported in the
  paper, whose table WA8 reports the likelihood-ratio tests above): per-estimate log-likelihood differences are summed within
  trials; *z* = mean / standard error across trials.
- **Implied distribution of true effects** (`Table_WA8/implied_distribution.R`; web appendix H): the shares of true effects above
  0.10 and 0.05 and below 0 and −0.05 under each fitted *t* model, θ = μ + τ·*T* with *T* ~ *t*(ν).
- **Influence of single trials on equation 2** (`Results_reported_in_text/loo_eq2.R`; body and web appendices J and R): equation 2
  with *t*-distributed true effects, at each choice of cutoffs, refitted without one trial at a time: for behavior, each of the 30
  trials most likely to drive the result (the 19 with significantly negative estimates at the 10% level, the 16 with 20 or more
  estimates, and the 3 with the largest inverse-variance weight); for knowledge, each of the three trials with significantly
  negative estimates; and, for behavior, Sconti, Caserta, and Ferrante (2024) (trial U-P081), which randomized one group of
  classes per arm (web appendix A). Each fit starts `nlminb` from the full-data solution (and from a second nearby point), keeps
  the better maximum, and restarts from it with `nlminb`'s default tolerances, as `fit_t` does; standard errors are not computed.
- **Each trial counted once** (`Results_reported_in_text/weighted_eq2.R`; body and web appendices J and R): the log-likelihood
  contribution of each estimate is weighted by one over the number of estimates from its trial (weights rescaled to sum to the
  number of estimates), for equation 2 and for the model without selection, with *t*-distributed true effects; standard errors are
  clustered by trial (sandwich with the weighted scores).
- **Nonparametric GMM estimator** (`gmm_functions.R`, `gmm_PublicationbiasPackage.R`; table WA9): Andrews and Kasy's online
  supplement, section C.1, using the authors' R code as distributed with Kaiser et al.'s replication files. We made four changes:
  we fixed a typo that stops the symmetric model with two cutoffs from running; we solve the second-stage moments for the mean and
  variance in closed form; we start the first-stage minimization from the authors' starting values and from 0.2 and keep the
  smaller objective; and we use a generalized inverse (`MASS::ginv`) of the moment covariance matrix in the continuously updated
  objective, which matters only when that matrix is nearly singular. On Kaiser et al.'s data as published, the code reproduces the
  GMM estimates in the earlier version of the paper.
- **RoBMA-PSMA** (`robma_functions.R`): RoBMA 4.0.0 with its default 36-model ensemble. The multilevel version nests estimates in
  trials. The trial-averaged version uses one estimate per trial (the mean of the trial's estimates, with a standard error that
  assumes a correlation of 0.5 between estimates of the same trial). Three chains with 5,000 adaptation and 10,000 burn-in
  iterations, then 60,000 (multilevel) or 50,000 (trial-averaged) iterations per chain; seed 1. Alternative priors (table WA13): a
  near-zero normal(0, 0.05) null prior on the effect in place of the spike at zero (95% of its mass between −0.1 and 0.1), or an
  effect prior centered at 0.1 (behavior) or 0.2 (knowledge) with the default standard deviation 1/√2 (RoBMA 4.0.0's default
  alternative is normal(0, 0.71)).
- **Funnel-plot test** (`Figure_WA2/funnel.R`; web appendices J and K): regression test (metafor `regtest`) on trial averages of
  *g* and of its standard error (the simple mean of the trial's standard errors, unlike RoBMA's trial averages), for the combined
  data and Kaiser et al.'s data alone (`output/funnel_tests.csv`, which also has the knowledge data, not reported).
- **Additional corrections** (`Table_WA10/1_corrections.R`; web appendix J): UWLS, WAAP, PET, and PEESE are weighted least squares
  with weights 1/SE² and CR2 standard errors by trial; PET-PEESE uses PEESE when the PET intercept is significantly positive
  (one-tailed .05). Trim-and-fill is metafor's `trimfill` on the two-level random-effects model. Power is computed from each
  estimate's standard error with the normal approximation (`output/power.csv`); `data/benchmarks.csv` holds the corrected
  estimates used as power benchmarks, copied from the outputs of `Figure_3` and `Figure_WA6`.
  `output/inverse_variance_weight_shares.csv` gives the share of the inverse-variance weight carried by the largest trial and the
  three largest trials (body and web appendix J). `3_table_WA10.py` assembles table WA10 from `output/corrections.csv` and, for the
  3PSM and 4PSM, `output/selection_models.csv`.
- **Excess mass** (`Table_WA3`, web appendix C): *Z* statistics are grouped in bins of width *h* = 0.05, 0.10, or 0.40 with an edge
  at the threshold *c* and 8 bins on each side of it, so that the window extends 0.4, 0.8, or 3.2 on each side. The bunching region
  is the first bin above the threshold, [*c*, *c* + *h*). A cubic polynomial fitted to the counts of the other bins of the window,
  leaving out also the bin that covers the bunching region of the other conventional threshold (1.645 or 1.96), gives the
  counterfactual, which is shifted to match the number of estimates in the window (Allen, Dechow, Pope, and Wu 2017, following
  Chetty et al. 2011). Excess mass is the observed count in the bunching region minus the counterfactual count, relative to the
  counterfactual count; its standard error comes from 500 bootstrap samples of trials (seed 1, drawn anew for each test), and
  *p*-values are one-sided. With bins of 0.40, the first bin above 1.645 also contains 1.96, so that test is not run. The script
  also runs placebo thresholds (1.30 and 2.30, bins of 0.10) and the tests without the estimates whose |*Z*| lies within 0.002
  above 1.645, 1.96, 2.326, or 2.576; most of them have standard errors derived from significance stars, which place them at the
  threshold by construction. `excess_mass_robustness.py` repeats the test at 1.96 for behavior without significance-derived
  standard errors, over 80 specifications (bins of 0.05, 0.10, 0.15, 0.20, or 0.25; polynomials of degree 1 to 4; 6, 8, 10, or 12
  bins on each side), and, for comparison, with all estimates and without the estimates at a threshold.
  `data/C_full_se_source.csv` is the behavior sample with `se_from_significance` = 1 for the 56 estimates of the newer trials
  whose standard error was derived from a reported *p*-value or from significance stars, as recorded in our coding notes;
  estimates whose |*Z*| lies within 0.002 of 1.645, 1.96, 2.326, 2.576, or 3.291 are dropped as well, because Kaiser et al.'s file
  does not record how each standard error was obtained and such a value marks one derived from a significance level.
- **Histograms of *Z*** (figures 2 and WA5; `z_histogram` in `figure_helpers.py`): bins of 0.40 closed on the right, as in R and
  JASP, with an edge at *Z* = 2.00, so that one bin runs from 1.60 to 2.00 and contains both positive critical values; with bins of
  0.40, that placement is the one at which the bin containing 1.96 stands out most from the two bins next to it, in Kaiser et al.'s
  data and in the combined data alike (estimates of effects on behavior). The script writes the counts in that bin and the bins
  next to it to its log. The excess-mass tests above use their own bins.
- **Summary of the corrections** (tables 4 and WA20): `Table_4_and_Table_WA20/data/` holds copies of outputs of other folders:
  `random_effects_behavior.csv` and `andrews_kasy_t_behavior.csv` (`Figure_3`), `random_effects_knowledge.csv` and
  `andrews_kasy_t_knowledge.csv` (`Figure_WA6`), `andrews_kasy_t_and_normal_behavior.csv` (`Table_WA8`,
  `table_WA8_model_details.csv`), `andrews_kasy_t_and_normal_knowledge.csv` (`Table_WA18`, `table_WA18_all_models.csv`), `gmm.csv`
  (`Table_WA9`), `corrections.csv` and `selection_models.csv` (`Table_WA10`), and `robma.csv` (`Table_WA13`). For each Andrews and
  Kasy specification, the script gives the estimate and 95% confidence interval at each choice of cutoffs (1.645, 1.96, and both);
  for knowledge, the GMM row uses the three-region version of equation 2, as in table WA9.
- **Specification curves** (figures 4 and WA7): every corrected estimate of the corrections above, for the combined data and for
  Kaiser et al.'s data alone, 59 per outcome; the folder's README lists the specifications and the source of each estimate.
- **Cost-effectiveness** (`Table_WA16/table_WA16.py`; web appendix P): the costs and the dollar effects, with their standard
  errors, of the 15 programs that report both are entered in the script, with, for each effect, whether we count it in the
  program's benefit (rules in web appendix P). Gine and Mansuri (2014) report savings in logs; the script converts the effect to
  dollars as the authors convert their effect on expenditure, from the control mean of the log. The benefit is the sum of the
  counted effects, and *g* is that sum divided by its standard deviation, from each effect's standard error and the median
  standard error of *g* of the trial's estimates in the combined data (`data/C_full.csv`). The corrected *g* and benefit-cost
  ratio divide the reported ones by six; the breakeven *g* is *g* × cost / benefit. The log also gives the counts and medians when
  benefits are divided by five or by four, and when the loan amount (Bruhn and Zia 2013) and the amounts remitted (Gibson,
  McKenzie, and Zia 2014), which we do not count as benefits, are added.

## 6. Results reported in the text (`Results_reported_in_text`)

- `output/andrews_kasy_t_other_samples.csv`: Andrews and Kasy models with *t*-distributed true effects for the data without
  mechanisms, and without mechanisms and economic outcomes (web appendix A); Kaiser et al.'s behavior data alone (body and web
  appendix G); objective knowledge alone and Kaiser et al.'s knowledge data alone (web appendix R); financial behaviors only (body
  and web appendix F).
- `output/random_effects_other_samples.csv`: the uncorrected estimates for the same data.
- `output/robma_summaries.csv` and `output/robma_*_models.csv`: multilevel and trial-averaged RoBMA-PSMA for these data sets.
- `output/gmm_other_samples.csv`: nonparametric GMM estimates for the financial behaviors only.
- `output/loo_eq2_behavior_targeted.csv` and `output/loo_eq2_knowledge_signeg.csv` (with `_full.csv`, the full-data fits):
  equation 2 leaving out one trial at a time (body and web appendices J and R); `output/loo_eq2_behavior_sconti.csv` (with
  `_full.csv`): equation 2 without Sconti, Caserta, and Ferrante (2024) (web appendix A).
- `output/weighted_eq2_behavior.csv` and `output/weighted_eq2_knowledge.csv`: equation 2 and the model without selection with each
  trial counted once (body and web appendices J and R).
- In other folders: `Table_WA2/output/table_WA2_paper_level.csv` (random-effects models with a level for papers, web appendix B);
  `Table_WA8/output/implied_distribution.csv` (web appendix H); `Table_WA10/output/power.csv` and
  `Table_WA10/output/inverse_variance_weight_shares.csv` (body and web appendix J); `Figure_WA2/output/funnel_tests.csv` (web
  appendices J and K); `Web_appendix_R_knowledge_subgroups/output/knowledge_subgroups.csv` and the knowledge rows of
  `Table_WA15/output/table_WA15.csv` (web appendix R); `Figure_4_and_Figure_WA7/output/spec_curve_summary.csv` and
  `spec_curve_paper_estimates_rank.csv` (body and web appendix N); the logs of `Figure_2_and_Figure_WA5` (shares of significant
  estimates, body) and `Table_WA16` (body and web appendix P).

## 7. Tables and figures not in this package

These items do not use the data above. Their materials are in our earlier OSF project:

- Figure 1 (replications of Thaler's "beer on the beach" study; data shared by Zoey Chen, Derek Rucker, and Blake McShane);
- Figure WA3 (funnel plot of the randomized trials in Fernandes, Lynch, and Netemeyer 2014, web appendix K);
- Table WA4 (simulation of the selection models, web appendix D).

Figure WA1 (search, screening, and inclusion of the newer trials, web appendix A) counts records from our search and screening,
including the papers and reports that do not enter the analyses; these records are not part of this package.

## 8. Notes on numerical reproducibility

- The Andrews and Kasy and random-effects point estimates are deterministic, up to the convergence tolerance of the optimizer
  (for example, a model fitted without the cached fits of the other models, see below, can differ beyond the sixth decimal).
  Standard errors that rely on numerical second derivatives (the normal models of `ak_likelihood_KH.R`, and the information matrix
  of `ak_t.R`, which differentiates the analytic gradient numerically) can differ across machines (BLAS libraries) in the fourth
  decimal, which occasionally changes the last reported digit. When the optimizer of a normal model tries parameters at which the
  likelihood cannot be evaluated, `ak_likelihood_KH.R` prints those parameters (lines such as `[1] -3.68`); these lines do not
  affect the results.
- RoBMA uses a fixed seed. With the same versions of RoBMA and JAGS, the results are the same; with other versions, Monte Carlo
  error can change the last digit, and Bayes factors reported as "> 179,999" or "> 149,999" mean that no posterior draw came from
  the models without that component.
- The GMM first stage is a nonconvex minimization; small differences across platforms in the third decimal are possible.
- The Andrews and Kasy scripts and the GMM scripts store each fit in `output/cache/`, under a name that depends on the data and the
  specification, and reuse it when a script is run again, so an interrupted run can be resumed; delete that folder to refit
  everything. The package is distributed without these folders, so the first run fits every model.
- We gave the folders, scripts, and output files the numbering of this version of the paper after running the R scripts; the logs
  of the R scripts are those of our runs, under the new names. The Python scripts were run again after the renaming.

## 9. Contact

Daniel Fernandes, Católica Lisbon School of Business and Economics, daniel.fernandes@ucp.pt
