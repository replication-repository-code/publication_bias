# Likelihood-ratio tests that allow for dependence within trials

The likelihood of Andrews and Kasy's (2019) models treats the estimates as independent, but estimates from the same trial are
dependent. For nested models, the likelihood-ratio statistic then follows a weighted sum of χ²(1) variables, with weights equal to
the eigenvalues of (R H⁻¹ R′)⁻¹ (R V R′), where H⁻¹ is the model-based covariance of the larger model, V its cluster-robust
(sandwich) covariance, and R the restrictions (Vuong 1989; Rotnitzky and Jewell 1990). The scripts compare each statistic with that
distribution (saddlepoint approximation of the tail probability; Kuonen 1999). The weights are the design effects of the
restrictions; with independent estimates they would all be 1.

| Script | Results | Reported in |
|---|---|---|
| `1_eq2_vs_eq1.R` | Equation 2 against equation 1 (*t* model) at 1.645, 1.96, and both, for the combined behavior data (`C_full`), Kaiser et al.'s behavior data (`K_only`), and knowledge (`C_know`); the two restrictions of table 3 (p1 = 1; p2 = p3) for the combined data. Also the cluster-robust Wald test of the same restrictions and the clustered Vuong z, for comparison only (it is a test for non-nested models) | Body ("Comparing Fit of Equations 1 and 2") and notes of tables 2 and 3; web appendix G (Kaiser et al.'s data alone) |
| `2_table_WA8.R` | Each *t* model against the best model (equation 2 at both cutoffs) and each normal model against the *t* model with the same selection (statistic divided by the design effect of the degrees-of-freedom parameter, evaluated at the *t* estimate, and compared with a 50:50 mixture of χ²(0) and χ²(1), because the normal model lies on the boundary) | Last column of table WA8 and web appendix H |

Run from this folder: `Rscript 1_eq2_vs_eq1.R` (about 20 minutes), then `Rscript 2_table_WA8.R` (about 10 minutes). The scripts
use the package's Andrews and Kasy code unchanged (`ak_functions.R`, `ak_likelihood_KH.R`, `ak_t.R`), whose fits reproduce tables
2, 3, WA7, WA8, and WA18; `robust_helpers.R` recomputes the model-based and clustered covariance at each estimate exactly as
`ak_t.R` does (the script checks that the standard errors match). Outputs: `output/robust_lr_eq2_vs_eq1.csv`,
`output/robust_lr_table3.csv`, `output/robust_lr_vs_best_C_full.csv` (the last column of table WA8), and the logs.

Main results (p allowing for dependence within trials; in parentheses, p treating estimates as independent):

- Combined behavior data, equation 2 against equation 1: 1.645, LR = 126.70, p = 2 × 10⁻¹⁶ (3 × 10⁻²⁸); 1.96, LR = 112.73,
  p = 6 × 10⁻⁹ (3 × 10⁻²⁵); both, LR = 144.88, p = 5 × 10⁻¹¹ (3 × 10⁻³¹).
- Table 3 restrictions: p1 = 1, p = 3 × 10⁻¹⁷ (1.645) and 2 × 10⁻¹² (1.96); p2 = p3, p = 3 × 10⁻⁷ and 7 × 10⁻⁵.
- Kaiser et al.'s data alone: p = .0003 (1.645), .0007 (1.96), .0002 (both).
- Knowledge: p < 10⁻¹¹ at every cutoff.
- Table WA8: every *t* model fits worse than the best model (largest p = .018, equation 2 at 1.645); every normal model fits worse
  than the *t* model with the same selection (design effects 2.0 to 4.5; p < 10⁻¹⁶).
