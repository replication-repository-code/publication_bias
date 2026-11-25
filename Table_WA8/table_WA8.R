# Table WA8: Andrews and Kasy (2019) models with t-distributed and with normal true effects, combined data, with clustered
# Vuong tests against the minimum-AIC model and LR tests of t versus normal (p from a 50:50 mixture of chi2(0) and chi2(1)).  Run from this folder: Rscript table_WA8.R   (about 20-40 minutes)
# The last column of table WA8 (likelihood-ratio tests that allow for dependence within trials) comes from the folder
# Likelihood_ratio_tests_clustered (2_table_WA8.R).
source("ak_functions.R")
d <- read.csv("data/C_full.csv")
tab <- ak_table(d, models = c("t", "normal"), sample_name = "C_full")
out <- data.frame(distribution = tab$model, selection = spec_label[tab$spec], g = round(tab$mu, 3), SE = round(tab$se_mu, 3),
                  tau = round(tab$tau, 3), df = round(tab$df_t, 2), log_lik = round(tab$logL, 2), parameters = tab$npar,
                  AIC = round(tab$AIC, 1), vuong_z = round(tab$vuong_z, 2), vuong_p_one_sided = signif(tab$vuong_p, 3), best = tab$best,
                  LR_t_vs_normal = round(tab$LR_t_vs_normal, 1), p_t_vs_normal = signif(tab$p_t_vs_normal, 3))
write.csv(out, "output/table_WA8.csv", row.names = FALSE); write.csv(tab, "output/table_WA8_model_details.csv", row.names = FALSE); print(out)
