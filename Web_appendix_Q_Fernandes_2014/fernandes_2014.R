# Web appendix Q: Andrews and Kasy (2019) corrections for the 15 randomized trials in Fernandes, Lynch, and Netemeyer (2014),
# with normally distributed true effects.
# Data: data/fernandes_2014_rcts.csv, one row per trial: eff = the effect as a correlation r, se = its standard error,
#   int = trial number. The models are fitted in r and the estimates converted to g = 2r / sqrt(1 - r^2).
# output/equation_2.csv: asymmetric equation 2 with three regions, Z < 0, 0 <= Z < 1.645, and Z >= 1.645 (publication
#   probability fixed at 1). No estimate has |Z| >= 1.96 or Z < -1.645, so only this three-region version can be estimated.
# output/equation_1.csv: symmetric equation 1 with one cutoff, |Z| = 1.645.
# output/text_estimates.csv: the estimates in g with their 95% confidence intervals and the log likelihoods of these two models
# and of the model without selection, and the uncorrected (fixed-effect) estimate, as quoted in the text.
# No estimate lies exactly at a cutoff, so the regions match the intervals in the output files.
# The variance of true effects is estimated at zero with these data, so t-distributed true effects give the same estimates.
# Run from this folder: Rscript fernandes_2014.R   (seconds)
source("ak_functions.R")
d0 <- read.csv("data/fernandes_2014_rcts.csv", fileEncoding = "UTF-8-BOM")
d <- data.frame(yi = d0$eff, sei = d0$se, cluster = d0$int)
tab <- ak_table(d, models = "normal", spec_names = c("none", "sym1645", "asym1645_3"), sample_name = "Fernandes et al. (2014) RCTs")
to_g <- function(r) 2 * r / sqrt(1 - r^2)
se_g <- function(r, se) 2 * se / (1 - r^2)^1.5                     # delta method
one <- function(s) tab[tab$spec == s, ]
table_of <- function(s, regions) {
  r <- one(s); b <- as.numeric(strsplit(r$betas, "; ")[[1]]); sb <- as.numeric(strsplit(r$se_betas, "; ")[[1]])
  data.frame(quantity = c("Corrected r", "Corrected g", paste("Publication probability,", regions)),
             estimate = round(c(r$mu, to_g(r$mu), b), 3), se = round(c(r$se_mu, se_g(r$mu, r$se_mu), sb), 3))
}
eq2 <- table_of("asym1645_3", c("Z < 0", "0 <= Z < 1.645"))
eq1 <- table_of("sym1645", "|Z| < 1.645")
w <- 1 / d$sei^2; fe <- sum(w * d$yi) / sum(w); se_fe <- sqrt(1 / sum(w))
text <- rbind(
  data.frame(model = c("Equation 2, Z = 0 and 1.645", "Equation 1, |Z| = 1.645", "No selection"),
             spec = c("asym1645_3", "sym1645", "none"),
             g = sapply(c("asym1645_3", "sym1645", "none"), function(s) to_g(one(s)$mu)),
             ci_lo = sapply(c("asym1645_3", "sym1645", "none"), function(s) to_g(one(s)$mu) - 1.96 * se_g(one(s)$mu, one(s)$se_mu)),
             ci_hi = sapply(c("asym1645_3", "sym1645", "none"), function(s) to_g(one(s)$mu) + 1.96 * se_g(one(s)$mu, one(s)$se_mu)),
             logL = sapply(c("asym1645_3", "sym1645", "none"), function(s) one(s)$logL)),
  data.frame(model = "Uncorrected fixed-effect estimate", spec = "inverse-variance weighted mean", g = to_g(fe),
             ci_lo = to_g(fe) - 1.96 * se_g(fe, se_fe), ci_hi = to_g(fe) + 1.96 * se_g(fe, se_fe), logL = NA))
rownames(text) <- NULL
write.csv(eq2, "output/equation_2.csv", row.names = FALSE)
write.csv(eq1, "output/equation_1.csv", row.names = FALSE)
write.csv(text, "output/text_estimates.csv", row.names = FALSE)
write.csv(tab, "output/model_details.csv", row.names = FALSE)
cat("\nEquation 2, three regions (output/equation_2.csv)\n"); print(eq2); cat("\nEquation 1, |Z| = 1.645 (output/equation_1.csv)\n"); print(eq1)
cat("\nNumbers quoted in the text (g with 95% CI; log likelihood)\n"); print(text, digits = 4)
