# Table WA17: relative publication probabilities from equation 2 (t), financial knowledge, combined data.  Run: Rscript table_WA17.R
source("ak_functions.R")
d <- read.csv("data/C_know.csv")
tab <- ak_table(d, models = "t", spec_names = c("asym196", "asym1645", "asymboth"), sample_name = "C_know")
out <- tab[, c("spec", "betas", "se_betas")]; write.csv(out, "output/table_WA17.csv", row.names = FALSE); print(out)
