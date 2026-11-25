# Table WA15 (behavior) and the knowledge results quoted in web appendix R: subgroup estimates with ONE selection function
# per dimension (Andrews and Kasy 2019, "conditioning on covariates"):
# subgroup-specific locations of the t distribution of true effects, common scale, degrees of freedom, and publication probabilities.
# Behavior: equation 2, cutoff 1.645. Knowledge: three-category version (few significantly negative estimates).
# The columns "Uncorrected g" and "Corrected g, separate selection" of Table WA15 come from the folder
# Figure_WA4_and_Table_WA14 (behavior; for knowledge, Web_appendix_R_knowledge_subgroups).
# Run from this folder: Rscript table_WA15.R   (20-40 minutes)
source("ak_likelihood_KH.R"); source("ak_functions.R"); source("ak_general.R")
pooled <- function(D, var, levels_, lab, domain) {
  if (!is.null(levels_)) D <- D[D[[var]] %in% levels_, ]
  g <- factor(D[[var]], levels = unique(D[[var]])); M <- model.matrix(~ 0 + g)
  if (domain == "knowledge") { cu <- c(0, 1.645); mp <- c(1, 2, 0) } else { cu <- c(-1.645, 0, 1.645); mp <- c(1, 2, 3, 0) }
  f <- fit_gen(D$yi, D$sei, as.character(D$cluster), cu, mp, "t", M)
  z <- D$yi / D$sei
  data.frame(domain = domain, dimension = lab, group = levels(g), k = as.vector(table(g)),
             trials = as.vector(tapply(D$cluster, g, function(x) length(unique(x)))),
             sig_negative_estimates = as.vector(tapply(z < -1.645, g, sum)),
             sig_negative_trials = as.vector(tapply(seq_along(z), g, function(i) length(unique(D$cluster[i][z[i] < -1.645])))),
             g = f$P[1:ncol(M)], se = f$se[1:ncol(M)], lo = f$P[1:ncol(M)] - 1.96 * f$se[1:ncol(M)], hi = f$P[1:ncol(M)] + 1.96 * f$se[1:ncol(M)],
             common_probabilities_SE = fmt_b(f), logL = f$logL, AIC = -2 * f$logL + 2 * f$npar, conv = f$conv)
}
B <- read.csv("data/C_full.csv"); K <- read.csv("data/C_know.csv")
B$outcome4 <- ifelse(B$otype %in% c("B", "C", "D"), c(B = "Borrowing", C = "Budgeting and planning", D = "Saving and investing")[B$otype],
                     "Other behaviors and downstream outcomes")
out <- list(pooled(B, "outcome4", NULL, "Outcome", "behavior"),
            pooled(B, "interv", c("Classroom", "Counseling", "Nudge", "Online"), "Intervention", "behavior"),
            pooled(B, "lowinc", NULL, "Sample", "behavior"),
            pooled(B, "topj", NULL, "Journal", "behavior"),
            pooled(K, "interv", c("Classroom", "Online"), "Intervention", "knowledge"),
            pooled(K, "lowinc", NULL, "Sample", "knowledge"),
            pooled(K, "topj", NULL, "Journal", "knowledge"))
r <- do.call(rbind, out); write.csv(r, "output/table_WA15.csv", row.names = FALSE); print(r, digits = 3)
