# Andrews-Kasy (2019) GMM estimator (nonparametric in the distribution of true effects),
# using the moment functions of the authors' PublicationbiasPackage.R (Kaiser et al. 2022 replication files).
src <- readLines("gmm_PublicationbiasPackage.R")
# fix a typo in the original code for symmetric models with several cutoffs (ifelse without yes/no values)
src <- gsub("(ifelse(abs(T)>=cutoffs[m-1]))", "(ifelse(abs(T)>=cutoffs[m-1],1,0))", src, fixed = TRUE)
src <- gsub("solve(Sigma_hat)", "MASS::ginv(Sigma_hat)", src, fixed = TRUE)  # robust to near-singular moment covariance
keep <- c("Step_function_normal_cdf","MetastudyMoments","ComputingMetastudyMoments","MetastudyGMMObjective","ComputingMetastudyGMMObjective","Clustered_covariance_estimate")
starts <- grep("^[A-Za-z_]+ <- function", src); names(starts) <- sub(" <- function.*", "", src[starts])
ends <- c(starts[-1] - 1, length(src))
for (k in keep) eval(parse(text = src[starts[k]:ends[which(names(starts) == k)]]), envir = globalenv())
gmm_fit <- function(X, sigma, cl, cutoffs, symmetric) {
  X <- as.matrix(X); sigma <- as.matrix(sigma); cl <- as.matrix(as.numeric(factor(cl)))
  nb <- length(cutoffs)
  obj1 <- function(P) MetastudyGMMObjective(P, cutoffs, symmetric, X, sigma, cl) + max(-min(P), 0) * 1e5
  # the continuously updating objective can have several local minima (e.g. a publication probability collapsing to zero);
  # start from all-ones (the authors' choice) and from 0.2 and keep the solution with the smaller objective
  m1 <- NULL
  for (st in list(rep(1, nb), rep(0.2, nb))) {
    m <- optim(par = st, fn = obj1, method = "BFGS", control = list(abstol = 1e-8, maxit = 1e5))
    if (is.null(m1) || m$value < m1$value) m1 <- m }
  b <- m1$par
  # Second stage: the mean and variance moments are exactly identified given the publication probabilities,
  # so solve them in closed form (identical to minimizing the authors' second-stage objective, but numerically stable)
  Tz <- as.vector(X / sigma); pz <- rep(1, length(Tz)); bb <- c(b, 1)
  if (symmetric == 1) { idx <- findInterval(abs(Tz), c(0, cutoffs, Inf)) } else { idx <- findInterval(Tz, c(-Inf, cutoffs, Inf)) }
  pz <- bb[idx]; w <- 1 / pz
  theta <- sum(w * X) / sum(w); tau2 <- sum(w * ((X - theta)^2 - sigma^2)) / sum(w)
  th <- c(theta, sqrt(max(tau2, 0))); m2 <- list(convergence = 0); all <- c(th, b)
  mom <- ComputingMetastudyMoments(c(all, 1), cutoffs, symmetric, X, sigma)
  S <- Clustered_covariance_estimate(mom$raw_moments, cl)
  G <- matrix(0, length(mom$moment_mean), length(all)); h <- 1e-3
  for (j in seq_along(all)) { p <- all; p[j] <- p[j] + h; m <- all; m[j] <- m[j] - h
    G[, j] <- (ComputingMetastudyMoments(c(p, 1), cutoffs, symmetric, X, sigma)$moment_mean - ComputingMetastudyMoments(c(m, 1), cutoffs, symmetric, X, sigma)$moment_mean) / (2 * h) }
  V <- MASS::ginv(t(G) %*% solve(S) %*% G) / length(X)
  list(mu = all[1], se_mu = sqrt(V[1, 1]), tau = all[2], se_tau = sqrt(V[2, 2]), betas = b, se_betas = sqrt(diag(V))[-(1:2)], J1 = m1$value, conv = c(m1$convergence, m2$convergence))
}
