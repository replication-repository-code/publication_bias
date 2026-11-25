# Table WA9: parametric (maximum likelihood, t) versus nonparametric GMM (Andrews and Kasy 2019, online supplement C.1).
# For knowledge, the asymmetric rows use the three-category version of equation 2 (all negative estimates in one region).
# Run from this folder: Rscript table_WA9.R   (1-2 hours: the GMM first stage is slow for the asymmetric models)
source("ak_functions.R"); source("gmm_functions.R")
# Fits are cached in output/cache/ and reused if present (delete that folder to refit everything).
cache <- function(name, expr) { dir.create("output/cache", showWarnings = FALSE); f <- file.path("output/cache", paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f)); v <- expr; saveRDS(v, f); v }
gspec <- list(sym1645 = list(1.645, 1), asym1645 = list(c(-1.645, 0, 1.645), 0), sym196 = list(1.96, 1), asym196 = list(c(-1.96, 0, 1.96), 0),
              symboth = list(c(1.645, 1.96), 1), asymboth = list(c(-1.96, -1.645, 0, 1.645, 1.96), 0),
              asym1645_3 = list(c(0, 1.645), 0), asym196_3 = list(c(0, 1.96), 0), asymboth_3 = list(c(0, 1.645, 1.96), 0))
rows <- c("sym1645", "asym1645", "sym196", "asym196", "symboth", "asymboth"); out <- list()
for (s in c("C_full", "K_only", "C_know")) {
  d <- read.csv(paste0("data/", s, ".csv"))
  use <- if (s == "C_know") sub("^(asym.*)$", "\\1_3", rows) else rows
  par <- cache(paste0("ak_", s), ak_table(d, models = "t", spec_names = use, sample_name = s))
  for (k in use) {
    g <- cache(paste0("gmm_", s, "_", k), try(gmm_fit(d$yi, d$sei, d$cluster, gspec[[k]][[1]], gspec[[k]][[2]]), silent = TRUE)); p <- par[par$spec == k, ]
    out[[paste(s, k)]] <- data.frame(sample = s, spec = k, parametric_g = round(p$mu, 3), parametric_se = round(p$se_mu, 3),
      gmm_g = if (inherits(g, "try-error")) NA else round(g$mu, 3), gmm_se = if (inherits(g, "try-error")) NA else round(g$se_mu, 3),
      gmm_ci_lo = if (inherits(g, "try-error")) NA else round(g$mu - qnorm(0.975) * g$se_mu, 3), gmm_ci_hi = if (inherits(g, "try-error")) NA else round(g$mu + qnorm(0.975) * g$se_mu, 3),
      gmm_tau = if (inherits(g, "try-error")) NA else round(g$tau, 3), gmm_betas = if (inherits(g, "try-error")) NA else paste(round(g$betas, 3), collapse = "; "))
    print(out[[paste(s, k)]])
  }
}
write.csv(do.call(rbind, out), "output/table_WA9.csv", row.names = FALSE)
