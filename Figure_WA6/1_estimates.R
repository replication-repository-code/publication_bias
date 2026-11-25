# Figure WA6, step 1: uncorrected random-effects estimate and Andrews and Kasy (2019) t models.
# Run from this folder: Rscript 1_estimates.R   (about 10-20 minutes)
source("ak_functions.R"); source("re_functions.R")
d <- read.csv("data/C_know.csv")
re <- re_model(d, "C_know"); write.csv(re, "output/random_effects.csv", row.names = FALSE)
tab <- ak_table(d, models = "t", spec_names = setdiff(names(specs), "none"), sample_name = "C_know")
write.csv(tab, "output/andrews_kasy_t.csv", row.names = FALSE)
