# Influence of single trials on Andrews and Kasy's equation 2 (t-distributed true effects).
# For each cutoff specification, refit equation 2 leaving out one trial at a time, starting from the full-data solution.
# Run from this folder (writes output/loo_eq2_<label>.csv and output/loo_eq2_<label>_full.csv; see below):
#   Rscript loo_eq2.R data/C_full.csv behavior_targeted K1,K11,K26,K28,K29,K38,K47,K50,K6,K60,K66,K7,K8,U-P007,U-P013,U-P014,U-P021,U-P038,U-P042,U-P044,U-P054,U-P056,U-P058,U-P083,U-P088,U-P094,U-P122,U-P129,U-P135,U-P138
#     (the 30 trials most likely to drive the result: the 19 with significantly negative estimates at the 10% level, the 16 with
#     20 or more estimates, and the 3 with the largest inverse-variance weight; about 40 minutes; leave out the list to refit
#     without every trial, which takes several hours)
#   Rscript loo_eq2.R data/C_know.csv knowledge_signeg K13,U-P026,U-P129   (the three trials with significantly negative knowledge estimates)
#   Rscript loo_eq2.R data/C_full.csv behavior_sconti U-P081   (Sconti, Caserta, and Ferrante 2024, one group of classes per arm; about 1 minute)
suppressMessages(library(parallel))
source("ak_functions.R")
args <- commandArgs(TRUE); dfile <- args[1]; label <- args[2]
d <- read.csv(dfile); d$cluster <- as.character(d$cluster)
trials <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else sort(unique(d$cluster))
spec_names <- c("asym1645", "asym196", "asymboth")
warm_fit <- function(X, s, cu, sym, map, P0) {
  Tp <- Tpowers_fun(X / s, cu, sym); M <- matrix(1, length(X), 1)
  f <- function(P) { v <- try(t_fg(P, cu, sym, X, s, Tp, map, M, grad = FALSE)$LLH, silent = TRUE); if (inherits(v, "try-error") || !is.finite(v)) 1e10 else v }
  gr <- function(P) { v <- try(t_fg(P, cu, sym, X, s, Tp, map, M)$grad, silent = TRUE); if (inherits(v, "try-error") || any(!is.finite(v))) rep(0, length(P)) else v }
  nb <- max(c(0, map)); lower <- c(-Inf, 1e-6, 0.05, rep(0, nb)); upper <- c(Inf, Inf, 1000, rep(Inf, nb))
  best <- NULL
  for (st in list(P0, c(P0[1], P0[2] * 1.5, min(P0[3] + 1, 999), P0[-(1:3)]))) {
    o <- try(nlminb(st, f, gr, lower = lower, upper = upper, control = list(eval.max = 5000, iter.max = 3000, abs.tol = 1e-10, rel.tol = 1e-12)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$objective < best$objective)) best <- o
  }
  # as in fit_t: restart from the solution with nlminb's default tolerances, which confirms the maximum and records convergence
  P <- best$par; logL <- -best$objective; conv <- best$convergence; msg <- best$message
  o <- try(nlminb(P, f, gr, lower = lower, upper = upper), silent = TRUE)
  if (!inherits(o, "try-error") && -o$objective >= logL - 1e-8) { P <- o$par; logL <- -o$objective; conv <- o$convergence; msg <- o$message }
  list(P = P, logL = logL, conv = conv, msg = msg)
}
full <- lapply(spec_names, function(k) { sp <- specs[[k]]; fit_one(d$yi, d$sei, d$cluster, sp, "t") }); names(full) <- spec_names
jobs <- expand.grid(spec = spec_names, trial = trials, stringsAsFactors = FALSE)
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(j) {
  k <- jobs$spec[j]; tr <- jobs$trial[j]; sp <- specs[[k]]; keep <- d$cluster != tr
  map <- c(seq_along(sp$cu), 0); z <- d$yi / d$sei
  o <- warm_fit(d$yi[keep], d$sei[keep], sp$cu, sp$sym, map, full[[k]]$P)
  if (o$conv != 0) {   # no convergence code from the warm start: refit with fit_t's starting values plus the warm start
    f <- fit_t(d$yi[keep], d$sei[keep], d$cluster[keep], sp$cu, sp$sym, map, starts_extra = list(full[[k]]$P))
    o <- list(P = f$P, logL = f$logL, conv = f$conv, msg = "refit with fit_t (no convergence from the warm start)")
  }
  data.frame(sample = label, spec = k, left_out = tr, n_left_out = sum(!keep), sig_neg_10_left_out = sum(z[!keep] < -1.645),
             mu = o$P[1], tau = o$P[2], df_t = o$P[3], betas = paste(sprintf("%.3f", o$P[-(1:3)]), collapse = "; "), logL = o$logL, conv = o$conv, message = o$msg)
}, mc.cores = 2)
out <- do.call(rbind, res)
write.csv(out, sprintf("output/loo_eq2_%s.csv", label), row.names = FALSE)
full_tab <- do.call(rbind, lapply(spec_names, function(k) data.frame(sample = label, spec = k, mu = full[[k]]$P[1], se_mu = full[[k]]$se[1], logL = full[[k]]$logL)))
write.csv(full_tab, sprintf("output/loo_eq2_%s_full.csv", label), row.names = FALSE)
cat(sprintf("%s: %d fits in %.1f min\n", label, nrow(out), as.numeric(Sys.time() - t0, units = "mins")))
for (k in spec_names) { r <- out[out$spec == k, ]; cat(sprintf("%-9s full mu = %.4f | leave-one-trial-out range %.4f to %.4f (min when leaving out %s, max when leaving out %s) | nonconverged %d\n",
  k, full[[k]]$P[1], min(r$mu), max(r$mu), r$left_out[which.min(r$mu)], r$left_out[which.max(r$mu)], sum(r$conv != 0))) }
