# Uncorrected random-effects model: estimates nested in trials (REML, metafor), CR2 cluster-robust
# standard errors with Satterthwaite degrees of freedom (clubSandwich), clustered by trial.
suppressPackageStartupMessages({library(metafor); library(clubSandwich)})
re_model <- function(d, sample_name = "") {
  d$es <- factor(seq_len(nrow(d)))
  f <- rma.mv(yi, sei^2, random = ~ 1 | cluster/es, data = d, method = "REML", sparse = TRUE)
  ci <- conf_int(f, vcov = "CR2", cluster = d$cluster); ct <- coef_test(f, vcov = "CR2", cluster = d$cluster)
  W <- diag(1 / f$vi); X <- model.matrix(f); P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- 100 * f$sigma2 / (sum(f$sigma2) + (f$k - f$p) / sum(diag(P)))
  data.frame(sample = sample_name, k = nrow(d), trials = length(unique(d$cluster)), mu = as.numeric(f$b), se_CR2 = ci$SE,
             lo = ci$CI_L, hi = ci$CI_U, df = ct$df_Satt, sigma2_trial = f$sigma2[1], sigma2_est = f$sigma2[2],
             I2_trial = I2[1], I2_est = I2[2], logLik = as.numeric(logLik(f)))
}
# Funnel-plot asymmetry test (Egger-type regression test, metafor::regtest) on trial averages of g and SE
funnel_test <- function(d, sample_name = "") {
  a <- aggregate(cbind(yi, sei) ~ cluster, d, mean); ra <- rma(yi, sei^2, data = a); rt <- regtest(ra)
  list(row = data.frame(sample = sample_name, trials = nrow(a), z = rt$zval, p = rt$pval, est_avg = as.numeric(ra$b)), fit = ra, agg = a)
}
