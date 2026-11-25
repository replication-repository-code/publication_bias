# Andrews and Kasy t models and uncorrected random effects for the data sets discussed in the text
# (web appendix A: subsets of outcomes; web appendix G: Kaiser et al.'s data alone; web appendix R: knowledge,
# objective knowledge alone, and Kaiser et al.'s knowledge data alone; web appendix F: financial behaviors only).
# Run from this folder: Rscript 1_andrews_kasy_and_random_effects.R   (1-2 hours)
source("ak_functions.R"); source("re_functions.R")
ak <- list(); re <- list()
for (s in c("C_noR6", "C_noR5R6", "K_only", "C_knowobj", "K_know", "C_R3")) {
  d <- read.csv(paste0("data/", s, ".csv")); ak[[s]] <- ak_table(d, models = "t", sample_name = s); re[[s]] <- re_model(d, s)
}
write.csv(do.call(rbind, ak), "output/andrews_kasy_t_other_samples.csv", row.names = FALSE)
write.csv(do.call(rbind, re), "output/random_effects_other_samples.csv", row.names = FALSE)
