# Table WA5 (column Uncorrected g [95% CI]) and the difference between journal articles and working papers.
# Run after 1_counts.py, from this folder: Rscript 2_random_effects.R
source("re_functions.R"); out <- list(); dif <- list()
for (dom in c("behavior", "knowledge")) {
  d <- read.csv(paste0("data/pubstatus_", dom, ".csv"))
  for (p in 1:0) { r <- re_model(d[d$published == p, ], paste(dom, ifelse(p == 1, "journal", "working paper")))
    out[[paste(dom, p)]] <- data.frame(domain = dom, published = p, k = r$k, trials = r$trials, g = r$mu, lo = r$lo, hi = r$hi) }
  d$es <- factor(seq_len(nrow(d))); f <- rma.mv(yi, sei^2, mods = ~ published, random = ~ 1 | cluster/es, data = d, method = "REML", sparse = TRUE)
  ct <- coef_test(f, vcov = "CR2", cluster = d$cluster)
  dif[[dom]] <- data.frame(domain = dom, difference_journal_minus_wp = ct$beta[2], se_CR2 = ct$SE[2], p_Satterthwaite = ct$p_Satt[2])
}
r <- do.call(rbind, out); write.csv(r, "output/table_WA5_random_effects.csv", row.names = FALSE); print(r, digits = 3)
dd <- do.call(rbind, dif); write.csv(dd, "output/table_WA5_difference.csv", row.names = FALSE); print(dd, digits = 3)
