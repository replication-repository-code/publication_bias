# Additional publication-bias corrections and diagnostics (web appendix J): UWLS, WAAP, PET, PEESE, PET-PEESE, three- and
# four-parameter selection models (3PSM, 4PSM) with normal true effects, trim-and-fill, leave-one-trial-out sensitivity, and
# study-level power, for behavior (C_full) and knowledge (C_know). Writes output/corrections.csv (rows of tables 4, WA10, and
# WA20; the 3PSM and 4PSM of these tables have t-distributed true effects and come from 2_selection_models.R), output/power.csv
# (statistical power, text of web appendix J), and output/inverse_variance_weight_shares.csv (text of the body and of web
# appendix J). Run from this folder:
#   Rscript 1_corrections.R
source("re_functions.R")   # loads metafor and clubSandwich; re_model()

cr2 <- function(fit, cl, coef = 1) {   # CR2 cluster-robust test of one coefficient (Satterthwaite df)
  ct <- coef_test(fit, vcov = "CR2", cluster = cl)[coef, ]
  est <- as.numeric(coef(fit)[coef]); se <- ct$SE; df <- ct$df_Satt
  c(g = est, se = se, lo = est - qt(.975, df) * se, hi = est + qt(.975, df) * se, p = ct$p_Satt, df = df)
}
res <- list(); power <- list(); wshare <- list()
for (s in c("C_full", "C_know")) {
  d <- read.csv(paste0("data/", s, ".csv")); d$w <- 1 / d$sei^2; cl <- d$cluster
  add <- function(method, v, note = "") res[[length(res) + 1]] <<- data.frame(sample = s, method = method, g = v[["g"]], se = v[["se"]], lo = v[["lo"]], hi = v[["hi"]], note = note)

  # Share of the inverse-variance weight (1/SE^2, used by UWLS, WAAP, PET and PEESE) carried by the largest trials
  tw <- sort(tapply(d$w, cl, sum) / sum(d$w), decreasing = TRUE)
  wshare[[s]] <- data.frame(sample = s, trials = length(tw), largest_trial = names(tw)[1], share_largest = tw[1],
                            three_largest = paste(names(tw)[1:3], collapse = "; "), share_three_largest = sum(tw[1:3]))
  # Uncorrected benchmark: random-effects model with estimates nested in trials (as in the body of the paper)
  re <- re_model(d, s); add("Uncorrected random-effects model", c(g = re$mu, se = re$se_CR2, lo = re$lo, hi = re$hi))

  # UWLS (Stanley and Doucouliagos 2015): weighted least squares with weights 1/SE^2 and a multiplicative error
  uwls <- lm(yi ~ 1, data = d, weights = w); u <- cr2(uwls, cl); add("UWLS", u, "CR2 standard errors by trial")

  # WAAP (Stanley, Doucouliagos, and Ioannidis 2017): UWLS on adequately powered estimates (SE <= |UWLS| / 2.8)
  ap <- d$sei <= abs(u[["g"]]) / 2.8
  if (sum(ap) >= 2 && length(unique(cl[ap])) >= 2) {
    wa <- cr2(lm(yi ~ 1, data = d[ap, ], weights = w), cl[ap])
    add("WAAP", wa, sprintf("%d adequately powered estimates from %d trials", sum(ap), length(unique(cl[ap]))))
  } else add("WAAP", c(g = NA, se = NA, lo = NA, hi = NA), sprintf("%d adequately powered estimates: WAAP not defined; UWLS applies", sum(ap)))

  # PET and PEESE (Stanley and Doucouliagos 2014): fixed-effect WLS of g on SE or SE^2, weights 1/SE^2
  pet <- lm(yi ~ sei, data = d, weights = w); pt <- cr2(pet, cl); ptb <- cr2(pet, cl, 2)
  add("PET", pt, sprintf("Slope on SE = %.2f (p %s)", ptb[["g"]], format.pval(ptb[["p"]], digits = 2, eps = .001)))
  peese <- lm(yi ~ I(sei^2), data = d, weights = w); pe <- cr2(peese, cl); add("PEESE", pe)
  # Conditional PET-PEESE: PEESE if the PET intercept is significantly positive (one-tailed .05), otherwise PET
  pet_pos <- pt[["g"]] > 0 && pt[["p"]] / 2 < .05
  add("PET-PEESE", if (pet_pos) pe else pt, if (pet_pos) "PET intercept significantly positive: PEESE" else "PET intercept not significantly positive: PET")

  # Three-parameter selection model (Vevea and Hedges 1995): random effects, one step at one-tailed p = .025
  m0 <- rma(yi, sei^2, data = d, method = "ML")
  sm <- selmodel(m0, type = "stepfun", steps = 0.025, alternative = "greater")
  add("Three-parameter selection model (3PSM)", c(g = as.numeric(sm$beta), se = sm$se, lo = sm$ci.lb, hi = sm$ci.ub),
      sprintf("Relative probability of nonsignificant or negative estimates = %.3f; LR test vs. no selection: chi2(1) = %.2f, p %s",
              sm$delta[2], sm$LRT, format.pval(sm$LRTp, digits = 2, eps = .001)))

  # Four-parameter selection model (McShane, Boeckenholt, and Hansen 2016): random effects, steps at one-tailed p = .025 and .50
  sm4 <- selmodel(m0, type = "stepfun", steps = c(0.025, 0.5), alternative = "greater")
  add("Four-parameter selection model (4PSM)", c(g = as.numeric(sm4$beta), se = sm4$se, lo = sm4$ci.lb, hi = sm4$ci.ub),
      sprintf("Relative probability of nonsignificant positive estimates = %.3f, of negative estimates = %.3f; LR test vs. no selection: chi2(2) = %.2f, p %s",
              sm4$delta[2], sm4$delta[3], sm4$LRT, format.pval(sm4$LRTp, digits = 2, eps = .001)))

  # Trim-and-fill (Duval and Tweedie 2000) on the two-level random-effects model
  m2 <- rma(yi, sei^2, data = d, method = "REML"); tf <- trimfill(m2)
  add("Trim-and-fill", c(g = as.numeric(tf$beta), se = tf$se, lo = tf$ci.lb, hi = tf$ci.ub),
      sprintf("%d estimates imputed on the %s side; two-level model before imputation g = %.3f", tf$k0, ifelse(tf$side == "left", "left", "right"), as.numeric(m2$beta)))

  # Leave-one-trial-out: uncorrected random-effects model without each trial in turn
  loo <- sapply(unique(cl), function(t) re_model(d[cl != t, ])$mu)
  add("Uncorrected, leave one trial out (range)", c(g = NA, se = NA, lo = min(loo), hi = max(loo)),
      sprintf("Range of %d estimates, each without one trial", length(loo)))

  # Study-level power: probability that a two-sided test at .05 is significant if the true effect equals the benchmark
  bm <- read.csv("data/benchmarks.csv"); bm <- rbind(data.frame(sample = s, benchmark = "Uncorrected random-effects model", g = re$mu), bm[bm$sample == s, ])
  for (i in seq_len(nrow(bm))) {
    pw <- pnorm(bm$g[i] / d$sei - qnorm(.975)) + pnorm(-bm$g[i] / d$sei - qnorm(.975))
    power[[length(power) + 1]] <- data.frame(sample = s, benchmark = bm$benchmark[i], g = bm$g[i], median_power = median(pw),
                                             mean_power = mean(pw), share_adequately_powered = mean(pw >= .80))
  }
}
out <- do.call(rbind, res); write.csv(out, "output/corrections.csv", row.names = FALSE)
write.csv(do.call(rbind, power), "output/power.csv", row.names = FALSE)
write.csv(do.call(rbind, wshare), "output/inverse_variance_weight_shares.csv", row.names = FALSE); print(do.call(rbind, wshare), digits = 3)
print(out, digits = 3); print(do.call(rbind, power), digits = 3)
