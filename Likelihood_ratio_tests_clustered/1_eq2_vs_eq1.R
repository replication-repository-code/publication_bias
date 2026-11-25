# Cluster-robust versions of the likelihood-ratio tests of equation 2 against equation 1 (Andrews and Kasy's t specification).
# For nested models, Vuong (1989, theorem 3.3) shows that under misspecification the LR statistic converges to a weighted sum of
# chi2(1) variables, with weights equal to the eigenvalues of (R H^-1 R')^-1 (R V R'), where H^-1 is the model-based covariance
# and V the sandwich covariance; with the score covariance clustered by trial, the same result covers dependence within trials.
# Reported for each comparison: the LR statistic, its usual chi2 p-value, the weights, the p-value from the weighted sum (Imhof),
# the cluster-robust Wald test of the same restrictions, and the clustered Vuong z of the package (trial sums of the pointwise
# log-likelihood differences), shown only for comparison: it is a test for non-nested models.
# Writes output/robust_lr_eq2_vs_eq1.csv and output/robust_lr_table3.csv. Run from this folder: Rscript 1_eq2_vs_eq1.R (about 20 minutes)
source("ak_functions.R")
source("robust_helpers.R")
# restrictions of equation 1 on equation 2: parameters (mu, tau, df, betas); R theta = r
restr <- list(
  "1.645" = list(full = "asym1645", sub = "sym1645", R = rbind(c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 0, 1, -1)), r = c(1, 0)),
  "1.96"  = list(full = "asym196", sub = "sym196", R = rbind(c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 0, 1, -1)), r = c(1, 0)),
  "both"  = list(full = "asymboth", sub = "symboth", R = rbind(c(0, 0, 0, 1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, 0, 0, -1), c(0, 0, 0, 0, 0, 1, -1, 0)), r = c(1, 0, 0)))
out <- list()
for (smp in c("C_full", "K_only", "C_know")) {
  d <- read.csv(paste0("data/", smp, ".csv")); X <- d$yi; s <- d$sei; cl <- as.character(d$cluster)
  for (k in names(restr)) {
    a <- restr[[k]]
    ff <- fit_one(X, s, cl, specs[[a$full]], "t"); fs <- fit_one(X, s, cl, specs[[a$sub]], "t")
    LR <- 2 * (ff$logL - fs$logL); q <- nrow(a$R)
    rp <- robust_parts(X, s, cl, specs[[a$full]], ff$P)
    stopifnot(max(abs(sqrt(diag(rp$V)) - ff$se)) < 1e-6)                       # same sandwich as the package's standard errors
    g <- as.vector(a$R %*% ff$P - a$r)
    W <- as.numeric(t(g) %*% solve(a$R %*% rp$V %*% t(a$R)) %*% g)
    lam <- Re(eigen(solve(a$R %*% rp$Hinv %*% t(a$R)) %*% (a$R %*% rp$V %*% t(a$R)), only.values = TRUE)$values)
    dd <- tapply(ff$pw - fs$pw, cl, sum); vz <- mean(dd) / (sd(dd) / sqrt(length(dd)))
    out[[paste(smp, k)]] <- data.frame(sample = smp, cutoffs = k, q = q, LR = LR, p_chi2 = pchisq(LR, q, lower.tail = FALSE),
      weights = paste(sprintf("%.2f", sort(lam, decreasing = TRUE)), collapse = "; "), mean_weight = mean(lam),
      LR_over_mean_weight = LR / mean(lam), p_weighted = weighted_p(LR, lam), p_rao_scott = pchisq(LR / mean(lam), q, lower.tail = FALSE),
      Wald_robust = W, p_Wald_robust = pchisq(W, q, lower.tail = FALSE), vuong_z = vz, vuong_p = pnorm(vz, lower.tail = FALSE))
    cat(sprintf("%-7s %-6s LR = %7.2f (q = %d)  weights = %-18s  p(weighted) = %.2e  robust Wald = %6.2f (p = %.2e)  Vuong z = %.2f\n",
                smp, k, LR, q, out[[paste(smp, k)]]$weights, out[[paste(smp, k)]]$p_weighted, W, pchisq(W, q, lower.tail = FALSE), vz))
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, "output/robust_lr_eq2_vs_eq1.csv", row.names = FALSE)
# the separate restrictions of table 3 (p1 = 1; p2 = p3) for the combined data: robust Wald and Rao-Scott-adjusted LR
d <- read.csv("data/C_full.csv"); X <- d$yi; s <- d$sei; cl <- as.character(d$cluster)
t3 <- list()
for (k in c("1.645", "1.96")) {
  a <- restr[[k]]; ff <- fit_one(X, s, cl, specs[[a$full]], "t"); rp <- robust_parts(X, s, cl, specs[[a$full]], ff$P)
  LRtab <- if (k == "1.645") c(119.83, 44.11) else c(103.14, 47.05)       # table 3 of the paper
  for (j in 1:2) {
    R <- a$R[j, , drop = FALSE]; g <- as.vector(R %*% ff$P - a$r[j])
    W <- as.numeric(g^2 / (R %*% rp$V %*% t(R))); lam <- as.numeric((R %*% rp$V %*% t(R)) / (R %*% rp$Hinv %*% t(R)))
    t3[[paste(k, j)]] <- data.frame(cutoffs = k, restriction = c("p1 = 1", "p2 = p3")[j], LR_table3 = LRtab[j], weight = lam,
      LR_adjusted = LRtab[j] / lam, p_adjusted = pchisq(LRtab[j] / lam, 1, lower.tail = FALSE), Wald_robust = W,
      p_Wald_robust = pchisq(W, 1, lower.tail = FALSE))
  }
}
t3 <- do.call(rbind, t3); rownames(t3) <- NULL; write.csv(t3, "output/robust_lr_table3.csv", row.names = FALSE)
print(res[, c("sample", "cutoffs", "q", "LR", "p_chi2", "weights", "p_weighted", "Wald_robust", "p_Wald_robust", "vuong_z", "vuong_p")])
print(t3)
