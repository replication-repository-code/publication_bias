# Figure WA6, step 2: multilevel RoBMA-PSMA (estimates nested in trials). Requires JAGS.
# Run from this folder: Rscript 2_robma.R   (about 1-3 hours). The fitted object is saved as output/robma_C_know_ml.rds;
# output/robma_C_know_ml_summary.csv is included, so step 3 can be run without refitting.
source("robma_functions.R")
d <- read.csv("data/C_know.csv")
f <- "output/robma_C_know_ml.rds"; fit <- if (file.exists(f)) readRDS(f) else fit_robma(d, type = "ml"); if (!file.exists(f)) saveRDS(fit, f)
write.csv(robma_summary(fit, "C_know multilevel"), "output/robma_C_know_ml_summary.csv", row.names = FALSE)
sink("output/robma_C_know_ml_printout.txt"); print(summary(fit)); sink()
