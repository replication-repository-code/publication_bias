# Nonparametric GMM estimates for the other data sets mentioned in the text (financial behaviors only, web appendix F).
# Run from this folder: Rscript 3_gmm.R   (1-2 hours)
source("gmm_functions.R")
# Fits are cached in output/cache/ and reused if present (delete that folder to refit everything).
cache <- function(name, expr) { dir.create("output/cache", showWarnings = FALSE); f <- file.path("output/cache", paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f)); v <- expr; saveRDS(v, f); v }
gspec <- list(sym1645 = list(1.645, 1), asym1645 = list(c(-1.645, 0, 1.645), 0), sym196 = list(1.96, 1), asym196 = list(c(-1.96, 0, 1.96), 0),
              symboth = list(c(1.645, 1.96), 1), asymboth = list(c(-1.96, -1.645, 0, 1.645, 1.96), 0)); out <- list()
for (s in c("C_R3")) { d <- read.csv(paste0("data/", s, ".csv"))
  for (k in names(gspec)) { g <- cache(paste0("gmm_", s, "_", k), try(gmm_fit(d$yi, d$sei, d$cluster, gspec[[k]][[1]], gspec[[k]][[2]]), silent = TRUE))
    out[[paste(s, k)]] <- data.frame(sample = s, spec = k, gmm_g = if (inherits(g, "try-error")) NA else g$mu, gmm_se = if (inherits(g, "try-error")) NA else g$se_mu)
    print(out[[paste(s, k)]]) } }
write.csv(do.call(rbind, out), "output/gmm_other_samples.csv", row.names = FALSE)
