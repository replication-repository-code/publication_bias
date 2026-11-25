# Andrews and Kasy's equation 2 (t-distributed true effects) with each trial counted once: each estimate's log-likelihood
# contribution is weighted by 1 / (number of estimates from its trial), rescaled so that the weights sum to the number of
# estimates. Standard errors: sandwich clustered by trial (weighted scores). Also the weighted model without selection.
# Run from this folder (writes output/weighted_eq2_<label>.csv; about 10 minutes each):
#   Rscript weighted_eq2.R data/C_full.csv behavior
#   Rscript weighted_eq2.R data/C_know.csv knowledge
source("ak_functions.R")
args <- commandArgs(TRUE); d <- read.csv(args[1]); label <- args[2]
X <- d$yi; s <- d$sei; cl <- as.character(d$cluster); n <- length(X)
nc <- table(cl); w <- as.numeric(1 / nc[cl]); w <- w * n / sum(w)
fit_t_w <- function(cu, sym, map, starts) {
  M <- matrix(1, n, 1); nb <- max(c(0, map)); Tp <- Tpowers_fun(X / s, cu, sym)
  f <- function(P) { v <- try(t_fg(P, cu, sym, X, s, Tp, map, M, grad = FALSE), silent = TRUE)
    if (inherits(v, "try-error")) return(1e10); val <- -sum(w * v$logL); if (!is.finite(val)) 1e10 else val }
  gr <- function(P) { v <- try(t_fg(P, cu, sym, X, s, Tp, map, M), silent = TRUE)
    if (inherits(v, "try-error")) return(rep(0, length(P))); g <- -colSums(w * v$score); if (any(!is.finite(g))) rep(0, length(P)) else g }
  lower <- c(-Inf, 1e-6, 0.05, rep(0, nb)); upper <- c(Inf, Inf, 1000, rep(Inf, nb)); best <- NULL
  for (st in starts) { o <- try(nlminb(st, f, gr, lower = lower, upper = upper, control = list(eval.max = 5000, iter.max = 3000, abs.tol = 1e-10, rel.tol = 1e-12)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$objective < best$objective)) best <- o }
  P <- best$par; o <- try(nlminb(P, f, gr, lower = lower, upper = upper), silent = TRUE)
  if (!inherits(o, "try-error") && o$objective <= best$objective + 1e-8) { P <- o$par; best$objective <- o$objective; best$convergence <- o$convergence }
  h <- 1e-6; Info <- sapply(seq_along(P), function(j) { pp <- P; pm <- P; pp[j] <- pp[j] + h; pm[j] <- pm[j] - h; (gr(pp) - gr(pm)) / (2 * h) })
  Info <- (Info + t(Info)) / 2; ev <- t_fg(P, cu, sym, X, s, Tp, map, M)
  Vr <- try(n * solve(Info) %*% Clustered_covariance_estimate(w * ev$score, cl) %*% solve(Info), silent = TRUE)
  se <- if (inherits(Vr, "try-error")) rep(NA, length(P)) else sqrt(pmax(diag(Vr), 0))
  list(P = P, se = se, objective = best$objective, conv = best$convergence)
}
rows <- list()
for (k in c("none", "asym1645", "asym196", "asymboth")) {
  sp <- specs[[k]]; uw <- fit_one(X, s, cl, sp, "t"); map <- if (sp$sel) c(seq_along(sp$cu), 0) else rep(0, length(sp$cu) + 1)
  nbk <- max(c(0, map)); P0 <- uw$P
  starts <- list(P0, c(P0[1], P0[2] * 1.5, 3, rep(0.5, nbk)), c(0.02, 0.05, 2, rep(0.5, nbk)), c(0.003, 0.1, 3, rep(1, nbk)))
  fw <- fit_t_w(sp$cu, sp$sym, map, starts)
  rows[[k]] <- data.frame(sample = label, spec = k, weighting = "each trial counted once", mu = fw$P[1], se_mu = fw$se[1],
    lo = fw$P[1] - 1.96 * fw$se[1], hi = fw$P[1] + 1.96 * fw$se[1], tau = fw$P[2], df_t = fw$P[3],
    betas = paste(sprintf("%.3f", fw$P[-(1:3)]), collapse = "; "), se_betas = paste(sprintf("%.3f", fw$se[-(1:3)]), collapse = "; "),
    conv = fw$conv, unweighted_mu = uw$P[1], unweighted_se_mu = uw$se[1])
  print(rows[[k]])
}
write.csv(do.call(rbind, rows), sprintf("output/weighted_eq2_%s.csv", label), row.names = FALSE)
