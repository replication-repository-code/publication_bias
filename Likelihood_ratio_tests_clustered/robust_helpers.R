# Helpers: tail probability of a weighted sum of chi2(1) variables (saddlepoint approximation, Kuonen 1999; it matches
# simulation to two or three digits), and the model-based and clustered sandwich covariance of an Andrews and Kasy t model at its estimate (as in ak_t.R).
weighted_p <- function(x, lambda) {         # P(sum lambda_j chi2_1 > x): saddlepoint approximation (Kuonen 1999), accurate in the tail
  K <- function(t) -0.5 * sum(log(1 - 2 * t * lambda)); K1 <- function(t) sum(lambda / (1 - 2 * t * lambda))
  K2 <- function(t) sum(2 * lambda^2 / (1 - 2 * t * lambda)^2)
  if (abs(x - sum(lambda)) < 1e-8) return(0.5)
  th <- uniroot(function(t) K1(t) - x, c(-1e6, (1 - 1e-12) / (2 * max(lambda))), tol = 1e-14)$root
  w <- sign(th) * sqrt(2 * (th * x - K(th))); v <- th * sqrt(K2(th))
  pnorm(w + log(v / w) / w, lower.tail = FALSE)
}
robust_parts <- function(X, s, cl, sp, P, sel = TRUE) {   # model-based and clustered sandwich covariance at P (as in fit_t)
  nb <- if (sel) length(sp$cu) else 0; map <- if (sel) c(seq_along(sp$cu), 0) else rep(0, length(sp$cu) + 1); M <- matrix(1, length(X), 1); Tp <- Tpowers_fun(X / s, sp$cu, sp$sym)
  gr <- function(Q) t_fg(Q, sp$cu, sp$sym, X, s, Tp, map, M)$grad
  h <- 1e-6
  Info <- sapply(seq_along(P), function(j) { pp <- P; pm <- P; pp[j] <- pp[j] + h; pm[j] <- pm[j] - h; (gr(pp) - gr(pm)) / (2 * h) })
  Info <- (Info + t(Info)) / 2
  ev <- t_fg(P, sp$cu, sp$sym, X, s, Tp, map, M)
  Cov <- Clustered_covariance_estimate(ev$score, cl)
  list(Hinv = solve(Info), V = length(X) * solve(Info) %*% Cov %*% solve(Info))
}
