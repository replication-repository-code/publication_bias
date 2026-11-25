# Table WA13: RoBMA-PSMA under alternative priors (trial averages) and the default multilevel model, behavior and knowledge.
# Run from this folder: Rscript table_WA13.R   (several hours; requires JAGS). Each fit is saved in output/ and reused if present.
source("robma_functions.R")
runs <- list(list("C_full", "avg", "default", NA), list("C_full", "avg", "nullnormal", NA), list("C_full", "avg", "effprior", 0.1), list("C_full", "ml", "default", NA),
             list("C_know", "avg", "default", NA), list("C_know", "avg", "nullnormal", NA), list("C_know", "avg", "effprior", 0.2), list("C_know", "ml", "default", NA))
out <- list()
for (r in runs) {
  lab <- paste(r[[1]], r[[2]], r[[3]], sep = "_"); f <- paste0("output/robma_", lab, ".rds")
  fit <- if (file.exists(f)) readRDS(f) else fit_robma(read.csv(paste0("data/", r[[1]], ".csv")), r[[2]], r[[3]], r[[4]])
  if (!file.exists(f)) saveRDS(fit, f); out[[lab]] <- robma_summary(fit, lab); print(out[[lab]][, c("model", "mu", "mu_lo", "mu_hi", "tau", "I2", "BF_effect", "BF_bias")]); rm(fit); gc()
}
write.csv(do.call(rbind, out), "output/table_WA13.csv", row.names = FALSE)
