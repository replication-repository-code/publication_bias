# General Andrews-Kasy likelihood model used for Tables 3, WA19, and WA15 (and the knowledge rows of table_WA15.csv).
# Requires ak_likelihood_KH.R, ak_t.R, and ak_functions.R to be sourced first.
# Regions are defined by ascending cutoffs on signed Z (asymmetric). `map` gives, for each region (lowest to highest Z),
# the index of its free publication probability, or 0 to fix it at 1 (the reference region). Regions sharing an index share
# one probability, so map = c(1, 2, 3, 0) is equation 2, c(0, 1, 2, 0) imposes p1 = 1, c(1, 2, 2, 0) imposes p2 = p3, and
# c(0, 1, 1, 0) is equation 1. M is a design matrix for the location (mean) of the distribution of true effects: an
# intercept by default, or subgroup indicators to estimate subgroup-specific means with one common selection function
# (Andrews and Kasy 2019, "conditioning on covariates"). The distribution of true effects is Andrews and Kasy's t
# (model = "t", ak_t.R) or normal (model = "normal", ak_likelihood_KH.R). Fits are cached in output/cache/.
fit_gen <- function(X, s, cl, cu, map, model = "t", M = NULL) {
  if (is.null(M)) M <- matrix(1, length(X), 1)
  ak_cached(paste0("gen_", model), list(X = X, s = s, cl = cl, cu = cu, map = map, model = model, M = M), {
    if (model == "t") fit_t(X, s, cl, cu, FALSE, map, M) else {
      Tp <- Tpowers_fun(X / s, cu, FALSE); nm <- ncol(M); nb <- max(map)
      unpack <- function(P) list(b = P[1:nm], tau = P[nm + 1],
                                 beta = { fb <- P[(length(P) - nb + 1):length(P)]; out <- rep(1, length(map)); out[map > 0] <- fb[map[map > 0]]; out })
      full <- function(P) { u <- unpack(P); VariationVarianceLogLikelihood(as.vector(M %*% u$b), u$tau, u$beta, cu, FALSE, X, s, Tp, df = Inf) }
      obj <- function(P) { v <- full(P)$LLH; if (!is.finite(v)) 1e10 else v }
      b0 <- qr.solve(M, X)
      starts <- lapply(list(c(b0, 0.1), c(rep(0, nm), 0.1), c(b0, sd(X))), function(st) c(st, rep(0.5, nb)))
      best <- NULL
      for (st in starts) {
        o <- try(nlminb(st, obj, lower = c(rep(-Inf, nm), rep(0, length(st) - nm)), control = list(eval.max = MaxEval, iter.max = MaxIter, abs.tol = Tol)), silent = TRUE)
        if (!inherits(o, "try-error") && (is.null(best) || o$objective < best$objective)) best <- o
      }
      P <- best$par; V <- try(RobustVariance(stepsize, length(X), P, full, cluster_ID = cl), silent = TRUE)
      se <- if (inherits(V, "try-error")) rep(NA, length(P)) else sqrt(pmax(diag(V), 0))
      list(P = P, se = se, logL = -best$objective, npar = length(P), nm = nm, nb = nb, conv = best$convergence)
    }
  })
}
# publication probabilities with standard errors, e.g. "0.011 (0.006); 0.237 (0.050); 0.396 (0.078)"
fmt_b <- function(f) { i <- (length(f$P) - f$nb + 1):length(f$P); paste(sprintf("%.3f (%.3f)", f$P[i], f$se[i]), collapse = "; ") }
