# Andrews and Kasy (2019) models for the specification curve, combined data (all outcomes) and Kaiser et al.'s data alone, for
# behavior and knowledge. Uses the paper's code (ak_functions.R, ak_likelihood_KH.R, ak_t.R) unchanged; standard errors clustered by
# trial.
#   normal: equations 1 and 2 at 1.645, 1.96, or both with normal true effects, and the model without selection (the package has
#           them for C_full and C_know, and 4_assemble.py checks that these fits reproduce them; K_only and K_know are new);
#   t:      the 3PSM (one step at Z = 1.96) and the 4PSM (steps at Z = 0 and 1.96) with t-distributed true effects (the package
#           has them for C_full and C_know, and 4_assemble.py checks these too; K_only and K_know are new).
# The t models with equations 1 and 2 are in the package's outputs (copied to data/). Rows are appended to output/ak_fits.csv as
# each model finishes; models already in that file are skipped.
# Usage: Rscript 1_andrews_kasy.R normal   (under a minute), then Rscript 1_andrews_kasy.R t   (about 10 minutes)
source("ak_functions.R")
a <- commandArgs(TRUE); model <- a[1]
samples <- if (length(a) > 1) a[-1] else c("C_full", "K_only", "C_know", "K_know")
sel <- if (model == "t") list(psm3 = list(cu = c(1.96), sym = FALSE, sel = TRUE), psm4 = list(cu = c(0, 1.96), sym = FALSE, sel = TRUE)) else specs
out <- "output/ak_fits.csv"
done <- if (file.exists(out)) read.csv(out) else NULL
for (s in samples) {
  d <- read.csv(paste0("data/", s, ".csv")); X <- d$yi; se <- d$sei; cl <- as.character(d$cluster)
  for (k in names(sel)) {
    if (!is.null(done) && any(done$sample == s & done$model == model & done$spec == k)) next
    t0 <- proc.time()[["elapsed"]]
    f <- fit_one(X, se, cl, sel[[k]], model)
    row <- data.frame(sample = s, k = length(X), clusters = length(unique(cl)), model = model, spec = k,
      mu = f$P[1], se_mu = f$se[1], tau = f$P[2], se_tau = f$se[2], df_t = if (model == "t") f$P[3] else NA,
      betas = paste(sprintf("%.4f", f$betas), collapse = "; "), se_betas = paste(sprintf("%.4f", f$se_betas), collapse = "; "),
      logL = f$logL, npar = f$npar, conv = f$conv, seconds = round(proc.time()[["elapsed"]] - t0, 1))
    write.table(row, out, sep = ",", row.names = FALSE, col.names = !file.exists(out), append = file.exists(out))
    cat(sprintf("%-9s %-6s %-9s mu = %8.5f (SE %.5f)  logL = %9.3f  %6.1fs\n", s, model, k, row$mu, row$se_mu, row$logL, row$seconds))
  }
}
