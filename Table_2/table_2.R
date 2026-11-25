# Table 2: log-likelihood fit of equations 1 and 2 (t-distributed true effects), combined data.
# Run from this folder: Rscript table_2.R   (about 10-20 minutes)
source("ak_functions.R")
d <- read.csv("data/C_full.csv")
tab <- ak_table(d, models = "t", sample_name = "C_full")
g <- function(s) tab[tab$spec == s, ]
out <- data.frame(model = spec_label[tab$spec], publication_probability_parameters = tab$npar - 3,
                  log_likelihood = round(tab$logL, 2), AIC = round(tab$AIC, 2), LR_vs_eq1 = NA, df = NA, p = NA)
for (p in list(c("asym1645", "sym1645"), c("asym196", "sym196"), c("asymboth", "symboth"))) {
  i <- tab$spec == p[1]; LR <- 2 * (g(p[1])$logL - g(p[2])$logL); ddf <- g(p[1])$npar - g(p[2])$npar
  out$LR_vs_eq1[i] <- round(LR, 2); out$df[i] <- ddf; out$p[i] <- signif(lr_p(LR, ddf, boundary = FALSE), 3) }
write.csv(out, "output/table_2.csv", row.names = FALSE); write.csv(tab, "output/table_2_model_details.csv", row.names = FALSE); print(out)
