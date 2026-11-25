# Andrews and Kasy's nonparametric GMM estimator (online supplement C.1) for the specification curve, for Kaiser et al.'s knowledge
# data alone (K_know), the one data set of the curve for which the package has no GMM estimates (table WA9 has C_full, K_only, and
# C_know). Uses the package's gmm_functions.R unchanged. Knowledge: three-region equation 2, as in table WA9. Rows are appended to
# output/gmm_fits.csv; fits already there are skipped.
# Usage: Rscript 2_gmm.R   (about a minute)
source("gmm_functions.R")
gspec <- list(sym1645 = list(1.645, 1), asym1645 = list(c(-1.645, 0, 1.645), 0), sym196 = list(1.96, 1), asym196 = list(c(-1.96, 0, 1.96), 0),
              symboth = list(c(1.645, 1.96), 1), asymboth = list(c(-1.96, -1.645, 0, 1.645, 1.96), 0),
              asym1645_3 = list(c(0, 1.645), 0), asym196_3 = list(c(0, 1.96), 0), asymboth_3 = list(c(0, 1.645, 1.96), 0))
a <- commandArgs(TRUE); samples <- if (length(a)) a else c("K_know")
know <- c("C_know", "K_know"); out <- "output/gmm_fits.csv"
for (s in samples) {
  d <- read.csv(paste0("data/", s, ".csv"))
  use <- c("sym1645", "asym1645", "sym196", "asym196", "symboth", "asymboth"); if (s %in% know) use <- sub("^(asym.*)$", "\\1_3", use)
  for (k in use) {
    done <- if (file.exists(out)) read.csv(out) else NULL      # read again before each fit (several runs can share the file)
    if (!is.null(done) && any(done$sample == s & done$spec == k)) next
    t0 <- proc.time()[["elapsed"]]
    g <- try(gmm_fit(d$yi, d$sei, d$cluster, gspec[[k]][[1]], gspec[[k]][[2]]), silent = TRUE)
    ok <- !inherits(g, "try-error")
    row <- data.frame(sample = s, spec = k, gmm_g = if (ok) g$mu else NA, gmm_se = if (ok) g$se_mu else NA, gmm_tau = if (ok) g$tau else NA,
                      gmm_betas = if (ok) paste(sprintf("%.4f", g$betas), collapse = "; ") else NA, J1 = if (ok) g$J1 else NA,
                      conv = if (ok) paste(g$conv, collapse = ";") else as.character(g), seconds = round(proc.time()[["elapsed"]] - t0, 1))
    write.table(row, out, sep = ",", row.names = FALSE, col.names = !file.exists(out), append = file.exists(out))
    cat(sprintf("%-9s %-11s g = %8.5f (SE %.5f)  %6.1fs\n", s, k, row$gmm_g, row$gmm_se, row$seconds))
  }
}
