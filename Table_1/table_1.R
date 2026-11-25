# Table 1: relative publication probabilities from asymmetric equation 2 (t-distributed true effects), combined data.
# Run from this folder: Rscript table_1.R   (about 5-10 minutes)
source("ak_functions.R")
d <- read.csv("data/C_full.csv")
tab <- ak_table(d, models = "t", spec_names = c("asym196", "asym1645", "asymboth"), sample_name = "C_full")
regions <- list(asym196 = c("(-inf, -1.96)", "[-1.96, 0)", "[0, 1.96)"), asym1645 = c("(-inf, -1.645)", "[-1.645, 0)", "[0, 1.645)"),
                asymboth = c("(-inf, -1.96)", "[-1.96, -1.645)", "[-1.645, 0)", "[0, 1.645)", "[1.645, 1.96)"))
out <- do.call(rbind, lapply(names(regions), function(s) { r <- tab[tab$spec == s, ]
  data.frame(panel = c(asym196 = "A: 5% level", asym1645 = "B: 10% level", asymboth = "C: both")[s], region = regions[[s]],
             estimate = round(as.numeric(strsplit(r$betas, "; ")[[1]]), 3), se = round(as.numeric(strsplit(r$se_betas, "; ")[[1]]), 3)) }))
write.csv(out, "output/table_1.csv", row.names = FALSE); write.csv(tab, "output/table_1_model_details.csv", row.names = FALSE); print(out)
