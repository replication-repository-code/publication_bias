# Table WA8 without the Vuong test: each t model against the best (minimum-AIC) model, the t model with equation 2 at both
# cutoffs, with the LR statistic compared with the weighted sum of chi2(1) that allows for dependence within trials (weights from
# the clustered sandwich of the best model); and the t-versus-normal LR tests with the boundary mixture, scaled by the design
# effect of the degrees-of-freedom parameter (evaluated at the t fit, an approximation).
# Writes output/robust_lr_vs_best_C_full.csv (the last column of table WA8). Run from this folder: Rscript 2_table_WA8.R (about 10 minutes)
source("ak_functions.R")
source("robust_helpers.R")
d <- read.csv("data/C_full.csv"); X <- d$yi; s <- d$sei; cl <- as.character(d$cluster)
best <- fit_one(X, s, cl, specs$asymboth, "t"); rp <- robust_parts(X, s, cl, specs$asymboth, best$P)
# parameters of the best model: mu, tau, df, b1..b5 for Z < -1.96, [-1.96, -1.645), [-1.645, 0), [0, 1.645), [1.645, 1.96)
e <- function(j) { v <- rep(0, 8); v[3 + j] <- 1; v }
R <- list(none = list(rbind(e(1), e(2), e(3), e(4), e(5)), rep(1, 5)),
          sym1645 = list(rbind(e(1), e(2), e(5), e(3) - e(4)), c(1, 1, 1, 0)),
          asym1645 = list(rbind(e(1) - e(2), e(5)), c(0, 1)),
          sym196 = list(rbind(e(1), e(2) - e(3), e(3) - e(4), e(4) - e(5)), c(1, 0, 0, 0)),
          asym196 = list(rbind(e(2) - e(3), e(4) - e(5)), c(0, 0)),
          symboth = list(rbind(e(1), e(2) - e(5), e(3) - e(4)), c(1, 0, 0)))
rows <- list()
for (k in names(R)) {
  f <- fit_one(X, s, cl, specs[[k]], "t"); A <- R[[k]][[1]]
  LR <- 2 * (best$logL - f$logL); q <- nrow(A)
  lam <- Re(eigen(solve(A %*% rp$Hinv %*% t(A)) %*% (A %*% rp$V %*% t(A)), only.values = TRUE)$values)
  # t versus normal with the same selection: design effect of df at the t fit
  if (k == "none") { sp <- specs$none } else sp <- specs[[k]]
  rk <- robust_parts(X, s, cl, sp, f$P, sel = sp$sel); de <- rk$V[3, 3] / rk$Hinv[3, 3]
  fn <- fit_one(X, s, cl, sp, "normal"); LRtn <- 2 * (f$logL - fn$logL)
  rows[[k]] <- data.frame(model = k, LR_vs_best = LR, q = q, p_chi2 = pchisq(LR, q, lower.tail = FALSE),
                          weights = paste(sprintf("%.2f", sort(lam, decreasing = TRUE)), collapse = "; "), p_clustered = weighted_p(LR, lam),
                          LR_t_vs_normal = LRtn, design_effect_df = de, p_t_vs_normal_clustered = 0.5 * pchisq(LRtn / de, 1, lower.tail = FALSE))
  cat(sprintf("%-9s LR vs best = %7.2f (q = %d)  weights = %-26s p = %.2e (indep. %.2e) | t vs normal LR = %.1f, design effect %.2f, p = %.1e\n",
              k, LR, q, rows[[k]]$weights, rows[[k]]$p_clustered, rows[[k]]$p_chi2, LRtn, de, rows[[k]]$p_t_vs_normal_clustered))
}
rb <- robust_parts(X, s, cl, specs$asymboth, best$P); de <- rb$V[3, 3] / rb$Hinv[3, 3]
fn <- fit_one(X, s, cl, specs$asymboth, "normal"); LRtn <- 2 * (best$logL - fn$logL)
cat(sprintf("asymboth  t vs normal LR = %.1f, design effect %.2f, p = %.1e\n", LRtn, de, 0.5 * pchisq(LRtn / de, 1, lower.tail = FALSE)))
res <- do.call(rbind, rows); rownames(res) <- NULL
res <- rbind(res, data.frame(model = "asymboth", LR_vs_best = NA, q = NA, p_chi2 = NA, weights = NA, p_clustered = NA, LR_t_vs_normal = LRtn,
                             design_effect_df = de, p_t_vs_normal_clustered = 0.5 * pchisq(LRtn / de, 1, lower.tail = FALSE)))
write.csv(res, "output/robust_lr_vs_best_C_full.csv", row.names = FALSE)
