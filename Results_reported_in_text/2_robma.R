# Multilevel and trial-averaged RoBMA-PSMA for the data sets discussed in the text. Requires JAGS; many hours in total.
# Fits are saved in output/ and reused if present. Run from this folder: Rscript 2_robma.R
source("robma_functions.R"); out <- list()
for (s in c("C_full", "K_only", "C_know", "K_know", "C_R3", "C_knowobj", "C_noR6", "C_noR5R6")) for (type in c("ml", "avg")) {
  if (s %in% c("C_knowobj", "C_noR6", "C_noR5R6") && type == "avg") next
  f <- paste0("output/robma_", s, "_", type, ".rds")
  fit <- if (file.exists(f)) readRDS(f) else fit_robma(read.csv(paste0("data/", s, ".csv")), type); if (!file.exists(f)) saveRDS(fit, f)
  out[[paste(s, type)]] <- robma_summary(fit, paste(s, type)); m <- robma_models(fit, type)
  write.csv(m, paste0("output/robma_", s, "_", type, "_models.csv"), row.names = FALSE); print(out[[paste(s, type)]]); rm(fit); gc()
}
write.csv(do.call(rbind, out), "output/robma_summaries.csv", row.names = FALSE)
