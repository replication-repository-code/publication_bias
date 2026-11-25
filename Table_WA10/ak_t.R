# Andrews and Kasy's (2019) t specification: t-distributed true effects with normally distributed sampling errors,
#   true effect theta ~ mu + tau * t(nu),   estimate X | theta ~ N(theta, SE^2)   (Andrews and Kasy 2019, section III.C).
# The t distribution is a scale mixture of normals: theta | V ~ N(mu, tau^2 / V) with V ~ Gamma(nu/2, rate nu/2), so
#   X | V ~ N(mu, SE^2 + tau^2 / V).
# The latent (pre-selection) density of each estimate and the probabilities of the selection regions are therefore integrals
# over V. They are computed with the trapezoid rule in u = log V on [-40, 5], with step 0.4 when nu <= 4 and a step that
# shrinks with the concentration of V otherwise (relative error below 1e-8 against adaptive quadrature). As nu grows, the
# model converges to the normal model of ak_functions.R (true effects ~ N(mu, tau^2)).
# Selection is a step function of Z = X / SE (regions from Tpowers_fun in ak_likelihood_KH.R; symmetric = regions of |Z|).
# `map` gives, for each region (lowest to highest Z), the index of its free publication probability, or 0 for a region whose
# probability is fixed at 1 (the reference region); regions sharing an index share one probability. M is a design matrix for
# the location of the distribution of true effects (an intercept by default; subgroup indicators give subgroup-specific means
# under one common selection function, Andrews and Kasy's "conditioning on covariates").
# Parameters: P = (location coefficients, tau, nu, free publication probabilities). The log-likelihood has an analytic gradient;
# standard errors are clustered by trial (sandwich: information from the derivative of the analytic gradient, analytic scores,
# and the clustered score covariance of ak_likelihood_KH.R).

t_grid <- function(nu, lo = -40, hi = 5) {
  a <- nu / 2; h <- if (a <= 2) 0.4 else 0.4 * sqrt(2 / a)
  u <- seq(lo, hi, by = h); V <- exp(u)
  lw <- log(h) + a * log(a) - lgamma(a) + a * u - a * V; keep <- lw > -700
  list(u = u[keep], V = V[keep], lw = lw[keep], a = a)
}

t_fg <- function(P, cu, sym, X, s, Tp, map, M, grad = TRUE) {
  nm <- ncol(M); nb <- max(c(0, map)); b <- P[1:nm]; tau <- P[nm + 1]; nu <- P[nm + 2]
  beta <- rep(1, length(map)); if (nb > 0) beta[map > 0] <- P[nm + 2 + map[map > 0]]
  lambda <- as.vector(M %*% b)
  G0 <- t_grid(nu); u <- G0$u; V <- G0$V; a <- G0$a
  w <- exp(G0$lw); dlw_dnu <- 0.5 * (log(a) + 1 - digamma(a) + u - V)
  n <- length(X); one <- rep(1, length(V))
  S <- sqrt(outer(s^2, one) + outer(rep(1, n), tau^2 / V))                 # n x K standard deviations of X | V
  dS_dtau <- outer(rep(1, n), tau / V) / S
  r <- (X - lambda) / S; g <- dnorm(r) / S
  fX <- as.vector(g %*% w)                                                  # latent density of each estimate
  Fpiece <- function(v) { q <- (v * s - lambda) / S; ph <- dnorm(q)          # P(Z* < v) and its derivatives
    list(F = as.vector(pnorm(q) %*% w), dl = as.vector((-ph / S) %*% w), dt = as.vector((-ph * q / S * dS_dtau) %*% w),
         dn = as.vector(pnorm(q) %*% (w * dlw_dnu))) }
  K <- length(cu); Fm <- Dl <- Dt <- Dn <- matrix(0, n, K)
  for (j in seq_len(K)) {
    if (sym) { p1 <- Fpiece(cu[j]); p2 <- Fpiece(-cu[j]); Fm[, j] <- p1$F - p2$F; Dl[, j] <- p1$dl - p2$dl; Dt[, j] <- p1$dt - p2$dt; Dn[, j] <- p1$dn - p2$dn }
    else { p1 <- Fpiece(cu[j]); Fm[, j] <- p1$F; Dl[, j] <- p1$dl; Dt[, j] <- p1$dt; Dn[, j] <- p1$dn }
  }
  cellF <- cbind(Fm, 1) - cbind(0, Fm)                                      # latent probability of each selection region
  cdiff <- function(Dm) cbind(Dm, 0) - cbind(0, Dm)
  nc <- as.vector(cellF %*% beta); phat <- as.vector(Tp %*% beta)
  lg <- log(phat) + log(fX) - log(nc)                                       # pointwise log-likelihood
  out <- list(LLH = -sum(lg), logL = lg)
  if (grad) {
    dfX_dl <- as.vector((g * r / S) %*% w); dfX_dt <- as.vector((g * (r^2 - 1) / S * dS_dtau) %*% w); dfX_dn <- as.vector(g %*% (w * dlw_dnu))
    gl <- dfX_dl / fX - as.vector(cdiff(Dl) %*% beta) / nc
    G <- cbind(gl * M, dfX_dt / fX - as.vector(cdiff(Dt) %*% beta) / nc, dfX_dn / fX - as.vector(cdiff(Dn) %*% beta) / nc)
    if (nb > 0) G <- cbind(G, sapply(seq_len(nb), function(k) { r_k <- which(map == k)
      rowSums(Tp[, r_k, drop = FALSE]) / phat - rowSums(cellF[, r_k, drop = FALSE]) / nc }))
    out$score <- G; out$grad <- -colSums(G)
  }
  out
}

# Fit one specification from several starting values and keep the maximum of the likelihood; then restart nlminb from the
# solution with its default tolerances, which confirms the maximum (it never accepts a lower likelihood) and records the
# convergence code, and compute clustered standard errors at the final point.
fit_t <- function(X, s, cl, cu, sym, map, M = NULL, starts_extra = list()) {
  if (is.null(M)) M <- matrix(1, length(X), 1)
  nm <- ncol(M); nb <- max(c(0, map)); n <- length(X); Tp <- Tpowers_fun(X / s, cu, sym)
  f <- function(P) { v <- try(t_fg(P, cu, sym, X, s, Tp, map, M, grad = FALSE)$LLH, silent = TRUE); if (inherits(v, "try-error") || !is.finite(v)) 1e10 else v }
  gr <- function(P) { v <- try(t_fg(P, cu, sym, X, s, Tp, map, M)$grad, silent = TRUE); if (inherits(v, "try-error") || any(!is.finite(v))) rep(0, length(P)) else v }
  b0 <- qr.solve(M, X)
  starts <- c(list(c(rep(0.003, nm), 0.1, 3, rep(1, nb)), c(rep(0.02, nm), 0.08, 2, rep(0.5, nb)), c(b0, sd(X), 10, rep(1, nb)),
                   c(rep(0.01, nm), 0.03, 1.5, rep(0.3, nb))), starts_extra)
  lower <- c(rep(-Inf, nm), 1e-6, 0.05, rep(0, nb)); upper <- c(rep(Inf, nm), Inf, 1000, rep(Inf, nb))
  best <- NULL
  for (st in starts) {
    o <- try(nlminb(st, f, gr, lower = lower, upper = upper, control = list(eval.max = 5000, iter.max = 3000, abs.tol = 1e-10, rel.tol = 1e-12)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$objective < best$objective)) best <- o
  }
  P <- best$par; logL <- -best$objective; conv <- best$convergence
  o <- try(nlminb(P, f, gr, lower = lower, upper = upper), silent = TRUE)
  if (!inherits(o, "try-error") && -o$objective >= logL - 1e-8) { P <- o$par; logL <- -o$objective; conv <- o$convergence }
  hstep <- 1e-6
  Info <- sapply(seq_along(P), function(j) { pp <- P; pm <- P; pp[j] <- pp[j] + hstep; pm[j] <- pm[j] - hstep; (gr(pp) - gr(pm)) / (2 * hstep) })
  Info <- (Info + t(Info)) / 2
  ev <- t_fg(P, cu, sym, X, s, Tp, map, M)
  Cov <- Clustered_covariance_estimate(ev$score, cl)
  Vr <- try(n * solve(Info) %*% Cov %*% solve(Info), silent = TRUE)
  se <- if (inherits(Vr, "try-error")) rep(NA, length(P)) else sqrt(pmax(diag(Vr), 0))
  i_b <- if (nb) (nm + 3):(nm + 2 + nb) else integer(0)
  list(P = P, se = se, logL = logL, conv = conv, pw = ev$logL, npar = length(P), nm = nm, nb = nb,
       betas = P[i_b], se_betas = se[i_b])
}
