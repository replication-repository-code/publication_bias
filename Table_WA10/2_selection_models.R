# Selection functions of the three- and four-parameter selection models combined with normal or with t-distributed true
# effects, behavior (C_full) and knowledge (C_know). 3PSM: one step at one-tailed p = .025 (Z = 1.96); 4PSM: steps at
# one-tailed p = .025 and .50 (Z = 1.96 and Z = 0). With normal true effects, these are metafor's selmodel step-function
# models (Vevea and Hedges 1995), estimates treated as independent; with t-distributed true effects, they are Andrews and
# Kasy's (2019) t specification (ak_t.R) with the same steps, standard errors clustered by trial.
# Run from this folder: Rscript 2_selection_models.R   (5-15 minutes)
source("ak_functions.R"); suppressMessages(library(metafor))
steps <- list("3PSM" = list(p = 0.025, cu = c(1.96)), "4PSM" = list(p = c(0.025, 0.5), cu = c(0, 1.96)))
out <- list()
for (s in c("C_full", "C_know")) {
  d <- read.csv(paste0("data/", s, ".csv")); m0 <- rma(yi, sei^2, data = d, method = "ML")
  for (k in names(steps)) {
    sm <- selmodel(m0, type = "stepfun", steps = steps[[k]]$p, alternative = "greater")
    pr <- rev(sm$delta[-1])   # probabilities relative to significant positive estimates, from the most negative region up
    out[[paste(s, k, "normal")]] <- data.frame(sample = s, selection = k, true_effects = "Normal", g = as.numeric(sm$beta), se = sm$se,
      lo = sm$ci.lb, hi = sm$ci.ub, tau = sqrt(sm$tau2), df = NA, publication_probabilities = paste(sprintf("%.3f", pr), collapse = "; "),
      logL = as.numeric(logLik(sm)), note = "metafor selmodel; estimates treated as independent")
    f <- fit_one(d$yi, d$sei, as.character(d$cluster), list(cu = steps[[k]]$cu, sym = FALSE, sel = TRUE), "t")
    out[[paste(s, k, "t")]] <- data.frame(sample = s, selection = k, true_effects = "t (Andrews and Kasy)", g = f$P[1], se = f$se[1],
      lo = f$P[1] - qnorm(0.975) * f$se[1], hi = f$P[1] + qnorm(0.975) * f$se[1], tau = f$P[2], df = f$P[3],
      publication_probabilities = paste(sprintf("%.3f", f$betas), collapse = "; "), logL = f$logL, note = "standard errors clustered by trial")
  }
}
r <- do.call(rbind, out); rownames(r) <- NULL; write.csv(r, "output/selection_models.csv", row.names = FALSE); print(r, digits = 3)
