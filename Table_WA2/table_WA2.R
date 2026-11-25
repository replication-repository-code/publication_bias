# Table WA2: uncorrected random-effects models (estimates nested in trials; CR2 standard errors).  Run: Rscript table_WA2.R
source("re_functions.R")
out <- do.call(rbind, lapply(c("C_full", "C_R3", "K_only", "C_know", "K_know"), function(s) re_model(read.csv(paste0("data/", s, ".csv")), s)))
write.csv(out, "output/table_WA2.csv", row.names = FALSE); print(out, digits = 3)
# Robustness (text of web appendix B): a separate random effect for papers, crossed with trials or with trials nested in papers
pl <- do.call(rbind, lapply(c("C_full", "C_know"), function(s) {
  d <- read.csv(paste0("data/", s, ".csv")); d$es <- factor(seq_len(nrow(d)))
  f2 <- rma.mv(yi, sei^2, random = ~ 1 | cluster/es, data = d, method = "REML", sparse = TRUE)
  f3a <- rma.mv(yi, sei^2, random = list(~ 1 | paper, ~ 1 | cluster/es), data = d, method = "REML", sparse = TRUE)
  f3b <- rma.mv(yi, sei^2, random = ~ 1 | paper/cluster/es, data = d, method = "REML", sparse = TRUE)
  data.frame(sample = s, model = c("Trials and estimates (body of the paper)", "Papers crossed with trials", "Trials nested in papers"),
             g = c(f2$b, f3a$b, f3b$b), sigma2 = c(paste(signif(f2$sigma2, 3), collapse = "; "), paste(signif(f3a$sigma2, 3), collapse = "; "),
             paste(signif(f3b$sigma2, 3), collapse = "; ")), logLik = c(logLik(f2), logLik(f3a), logLik(f3b))) }))
write.csv(pl, "output/table_WA2_paper_level.csv", row.names = FALSE); print(pl, digits = 4)
