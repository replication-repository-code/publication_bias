# Figure WA2 (funnel plot of trial-averaged effect sizes, behavior) and the regression tests for funnel-plot asymmetry
# quoted in web appendix K (combined data and Kaiser et al.'s data alone; knowledge too, not reported).
# Run from this folder: Rscript funnel.R
source("re_functions.R")
out <- list()
for (s in c("C_full", "K_only", "C_know", "K_know")) {
  ft <- funnel_test(read.csv(paste0("data/", s, ".csv")), s); out[[s]] <- ft$row
  if (s == "C_full") {
    png("output/figure_WA2.png", width = 6, height = 5.2, units = "in", res = 300)
    par(mar = c(4.5, 4.5, 1, 1), family = "serif")
    funnel(ft$fit, level = c(90, 95, 99), shade = c("white", "gray60", "gray80"), refline = 0, legend = list(show = "cis"),
           xlab = "Effect size (Hedges' g, trial average)", ylab = "Standard error", pch = 19,
           col = ifelse(grepl("^K", ft$agg$cluster), "gray25", "#1F5AA6"))
    legend("bottomright", legend = c("Kaiser et al. (2022)", "Newer trials"), pch = 19, col = c("gray25", "#1F5AA6"), bty = "n", cex = 0.85)
    dev.off()
  }
}
out <- do.call(rbind, out); write.csv(out, "output/funnel_tests.csv", row.names = FALSE); print(out, digits = 3)
