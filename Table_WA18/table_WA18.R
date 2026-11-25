# Table WA18: log-likelihood fit of equations 1 and 2 (t), financial knowledge, combined data; also the normal models
# of table WA20 (web appendix R).  Run: Rscript table_WA18.R   (about 15-30 minutes)
source("ak_functions.R")
d <- read.csv("data/C_know.csv")
tab <- ak_table(d, models = c("t", "normal"), sample_name = "C_know")
t <- tab[tab$model == "t", ]
out <- data.frame(model = spec_label[t$spec], log_likelihood = round(t$logL, 2), AIC = round(t$AIC, 2), g = round(t$mu, 3), SE = round(t$se_mu, 3))
write.csv(out, "output/table_WA18.csv", row.names = FALSE); write.csv(tab, "output/table_WA18_all_models.csv", row.names = FALSE)
print(out)
lr <- rbind(lr_eq2_vs_eq1(tab, "t"), lr_eq2_vs_eq1(tab, "normal")); write.csv(lr, "output/table_WA18_LR_eq2_vs_eq1.csv", row.names = FALSE); print(lr, digits = 4)
