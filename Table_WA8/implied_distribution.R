# Web appendix H: how true effects vary under the fitted t distributions (Andrews and Kasy's t specification, table WA8).
# For each t model, the share of true effects above 0.10 and 0.05, below 0, and below -0.05, from the fitted location (mu),
# scale (tau), and degrees of freedom: theta = mu + tau * T, T ~ t(df).
# Run from this folder after table_WA8.R: Rscript implied_distribution.R   (seconds; writes output/implied_distribution.csv)
d <- read.csv("output/table_WA8_model_details.csv"); d <- d[d$model == "t", ]
F <- function(x) pt((x - d$mu) / d$tau, d$df_t)
out <- data.frame(spec = d$spec, mu = d$mu, tau = d$tau, df_t = d$df_t, share_above_0.10 = 1 - F(0.10), share_above_0.05 = 1 - F(0.05),
                  share_below_0 = F(0), share_below_minus_0.05 = F(-0.05))
write.csv(out, "output/implied_distribution.csv", row.names = FALSE); print(out, digits = 3)
