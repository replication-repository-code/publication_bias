# Single-level RoBMA-PSMA that treats every estimate as independent (no clustering by trial), with the package's
# default ensemble and sampler settings, as a reader running RoBMA() with its defaults would (response to Reviewer 3 in the
# revision notes; not reported in the paper). This reproduces the estimate of about 0.043 for Kaiser et al.'s data; the
# paper's multilevel and trial-averaged models account for the dependence. With the default sampler (15,000 draws), these
# fits mix less well than the paper's (see mu_Rhat and mu_ESS in the output files).
# data/kaiser_raw.csv: the 458 behavior estimates (outcome codes B-G) as in Kaiser et al.'s file, before duplicates are removed.
# Run from this folder: Rscript robma_single_level.R   (about 10-20 minutes per data set)
source("robma_functions.R")
for (f in c("kaiser_raw", "K_only", "C_full")) {
  d <- read.csv(paste0("data/", f, ".csv"))
  fit <- RoBMA(yi = d$yi, sei = d$sei, measure = "SMD", parallel = TRUE, seed = 1, silent = TRUE)
  out <- robma_summary(fit, paste(f, "unclustered"))
  sc <- summary(fit, conditional = TRUE)$estimates_conditional
  out$mu_conditional <- if (!is.null(sc) && "mu" %in% rownames(sc)) sc["mu", "Mean"] else NA
  write.csv(out, paste0("output/unclustered_", f, ".csv"), row.names = FALSE); print(out); rm(fit); gc()
}
