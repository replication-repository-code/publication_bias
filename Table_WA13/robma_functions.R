# RoBMA-PSMA (RoBMA package 4.0.0; Bartos et al. 2023; multilevel version: Bartos, Maier, and Wagenmakers 2026)
suppressPackageStartupMessages(library(RoBMA))
# One estimate per trial: mean of the trial's estimates; its SE assumes a correlation of 0.5 between estimates of the same trial
trial_average <- function(d, r_within = 0.5) {
  do.call(rbind, lapply(split(d, d$cluster), function(x) {
    k <- nrow(x); v <- x$sei^2; V <- outer(sqrt(v), sqrt(v)) * r_within; diag(V) <- v
    data.frame(trial = x$cluster[1], yi = mean(x$yi), sei = sqrt(sum(V)) / k, k = k)
  }))
}
# type = "ml" (estimates nested in trials) or "avg" (trial averages). prior: "default", "nullnormal", or "effprior" (with effm).
# RoBMA 4.0.0 defaults for the effect: null = spike at 0; alternative = Normal(0, 1/sqrt(2)) (printed as Normal(0, 0.71)).
# "nullnormal" replaces the spike by a near-zero (perinull) prior, Normal(0, 0.05): 95% of its mass on effects between -0.1 and 0.1.
# "effprior" moves the centre of the alternative to effm and keeps its default standard deviation, 1/sqrt(2).
fit_robma <- function(d, type = "ml", prior = "default", effm = NA) {
  if (type == "avg") {
    a <- trial_average(d)
    args <- list(yi = a$yi, sei = a$sei, measure = "SMD", adapt = 5000, burnin = 10000, sample = 50000, parallel = TRUE, seed = 1, silent = TRUE)
  } else {
    args <- list(yi = d$yi, sei = d$sei, measure = "SMD", cluster = d$cluster, adapt = 5000, burnin = 10000, sample = 60000, parallel = TRUE, seed = 1, silent = TRUE)
  }
  if (prior == "nullnormal") args$prior_effect_null <- prior("normal", list(mean = 0, sd = 0.05))
  if (prior == "effprior") args$prior_effect <- prior("normal", list(mean = effm, sd = 1 / sqrt(2)))
  do.call(RoBMA, args)
}
# Compact summary of a fit: inclusion Bayes factors, model-averaged estimates, interquartile range of mu, convergence
robma_summary <- function(fit, label = "") {
  s <- summary(fit); comp <- s$inclusion_components; est <- s$estimates
  bf <- as.numeric(comp$inclusion_BF); bound <- attr(comp$inclusion_BF, "bound_operator"); if (is.null(bound)) bound <- rep(NA, 3)
  fmt_bf <- function(i) if (!is.na(bound[i])) paste0(">", round(bf[i])) else if (bf[i] >= 1000) as.character(round(bf[i])) else as.character(signif(bf[i], 4))
  mu_draws <- unlist(lapply(fit$fit$mcmc, function(m) as.matrix(m)[, "mu"]))   # only the column needed (saves memory)
  g <- function(par, col) if (par %in% rownames(est)) as.numeric(est[par, col]) else NA
  h <- if (isTRUE(attr(s, "multilevel"))) NULL else try(summary_heterogeneity(fit)$estimates, silent = TRUE)   # I2 for trial averages only
  data.frame(model = label,
    BF_effect = fmt_bf(1), post_prob_effect = as.numeric(comp$post_prob[1]), BF_heterogeneity = fmt_bf(2), BF_bias = fmt_bf(3),
    mu = g("mu", "Mean"), mu_lo = g("mu", "0.025"), mu_hi = g("mu", "0.975"),
    mu_q25 = unname(quantile(mu_draws, 0.25)), mu_q75 = unname(quantile(mu_draws, 0.75)),
    tau = g("tau", "Mean"), tau_lo = g("tau", "0.025"), tau_hi = g("tau", "0.975"),
    rho = g("rho", "Mean"), rho_lo = g("rho", "0.025"), rho_hi = g("rho", "0.975"),
    I2 = if (!is.null(h) && !inherits(h, "try-error") && "I2" %in% rownames(h)) as.numeric(h["I2", "Mean"]) else NA,
    mu_ESS = g("mu", "ESS"), mu_Rhat = g("mu", "R_hat"), n_draws = length(mu_draws))
}
# Posterior probability of each of the 36 models, Bayes factor against the most probable model (the model with the
# highest posterior probability), and mean log-likelihood
robma_models <- function(fit, type = "ml") {
  ind <- c("mu_indicator", "tau_indicator", "bias_indicator")
  mc <- do.call(rbind, lapply(fit$fit$mcmc, function(m) as.matrix(m)[, ind]))
  unit <- if (type == "ml") "cluster" else "estimate"; th <- if (type == "ml") 10 else 1
  f <- fit
  if (th > 1) for (j in seq_along(f$fit$mcmc)) { m <- f$fit$mcmc[[j]]; f$fit$mcmc[[j]] <- coda::mcmc(as.matrix(m)[seq(1, nrow(m), by = th), , drop = FALSE]) }
  ll <- RoBMA:::.log_lik.brma(f, unit = unit); mcl <- do.call(rbind, lapply(f$fit$mcmc, function(m) as.matrix(m)[, ind])); rm(f); gc()
  tot <- rowSums(ll); keyl <- paste(mcl[, "mu_indicator"], mcl[, "tau_indicator"], mcl[, "bias_indicator"])
  bias_lab <- c("None", "omega[two-sided: .05]", "omega[two-sided: .05, .1]", "omega[one-sided: .05]", "omega[one-sided: .025, .05]",
                "omega[one-sided: .05, .5]", "omega[one-sided: .025, .05, .5]", "PET", "PEESE")
  bias_prior <- c(0.5, rep(0.5 / 12, 6), 0.125, 0.125)
  key <- paste(mc[, "mu_indicator"], mc[, "tau_indicator"], mc[, "bias_indicator"])
  grid <- expand.grid(mu = c(1, 2), tau = c(1, 2), bias = 1:9)
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]; s <- key == paste(g$mu, g$tau, g$bias); sl <- keyl == paste(g$mu, g$tau, g$bias)
    data.frame(effect = ifelse(g$mu == 1, "No", "Yes"), heterogeneity = ifelse(g$tau == 1, "No", "Yes"), bias = bias_lab[g$bias],
               prior_prob = 0.25 * bias_prior[g$bias], n_draws = sum(s), post_prob = mean(s),
               mean_loglik = if (any(sl)) mean(tot[sl]) else NA)
  }))
  res$post_odds_ratio <- res$post_prob / res$prior_prob; b <- which.max(res$post_prob)
  res$BF_vs_best <- res$post_odds_ratio / res$post_odds_ratio[b]
  lo <- loo::loo(ll, r_eff = NA)
  attr(res, "elpd_loo") <- lo$estimates["elpd_loo", 1:2]
  res[order(-res$post_prob), ]
}
