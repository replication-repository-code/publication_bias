# Table WA7: Andrews and Kasy (2019) t models, Kaiser et al.'s behavior data alone.  Run from this folder: Rscript table_WA7.R   (about 10-20 minutes)
source("ak_functions.R")
d <- read.csv("data/K_only.csv")
tab <- ak_table(d, models = "t", sample_name = "K_only")
out <- data.frame(model = spec_label[tab$spec], g = round(tab$mu, 3), SE = round(tab$se_mu, 3),
                  ci_lo = round(tab$mu - 1.959964 * tab$se_mu, 3), ci_hi = round(tab$mu + 1.959964 * tab$se_mu, 3),
                  log_lik = round(tab$logL, 2), publication_probabilities = tab$betas, se_publication_probabilities = tab$se_betas)
write.csv(out, "output/table_WA7.csv", row.names = FALSE); print(out)
lr <- lr_eq2_vs_eq1(tab); write.csv(lr, "output/table_WA7_LR_eq2_vs_eq1.csv", row.names = FALSE); print(lr, digits = 4)
