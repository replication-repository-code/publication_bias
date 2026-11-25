# Table WA12: most probable RoBMA-PSMA models, multilevel (panel A) and trial-averaged (panel B), combined behavior data.
# Run from this folder: Rscript table_WA12.R   (about 1-3 hours; requires JAGS)
source("robma_functions.R")
d <- read.csv("data/C_full.csv")
for (type in c("ml", "avg")) {
  f <- paste0("output/robma_C_full_", type, ".rds")
  fit <- if (file.exists(f)) readRDS(f) else fit_robma(d, type); if (!file.exists(f)) saveRDS(fit, f)
  m <- robma_models(fit, type); write.csv(m, paste0("output/table_WA12_", type, "_models.csv"), row.names = FALSE)
  write.csv(robma_summary(fit, paste("C_full", type)), paste0("output/robma_C_full_", type, "_summary.csv"), row.names = FALSE)
  print(head(m, 5)); rm(fit); gc()
}
