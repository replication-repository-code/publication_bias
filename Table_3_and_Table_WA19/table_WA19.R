# Table WA19: publication probabilities and likelihood-ratio tests of the symmetry restrictions, combined knowledge data
# (C_know.csv; the knowledge version of table 3, web appendix R), t-distributed true effects, cutoffs at 1.645 and 1.96.
# Run from this folder: Rscript table_WA19.R   (1-2 minutes)
source("ak_likelihood_KH.R"); source("ak_functions.R"); source("ak_general.R")
d <- read.csv("data/C_know.csv"); X <- d$yi; s <- d$sei; cl <- as.character(d$cluster)
out <- list()
for (c0 in c(1.645, 1.96)) {
  cu <- c(-c0, 0, c0)    # regions: Z < -c (p1), [-c, 0) (p2), [0, c) (p3), Z >= c (reference, fixed at 1)
  m <- list("Eq. 2: p1, p2, p3 free" = c(1, 2, 3, 0),
            "Significant results alike: p1 = 1" = c(0, 1, 2, 0),
            "Nonsignificant results alike: p2 = p3" = c(1, 2, 2, 0),
            "Eq. 1: p1 = 1 and p2 = p3" = c(0, 1, 1, 0))
  fits <- lapply(m, function(mp) fit_gen(X, s, cl, cu, mp))
  for (k in names(fits)) {
    f <- fits[[k]]; LR <- 2 * (fits[[1]]$logL - f$logL); ddf <- fits[[1]]$npar - f$npar
    out[[paste(c0, k)]] <- data.frame(cutoff = c0, model = k, probabilities_SE = fmt_b(f), g = f$P[1], se = f$se[1], logL = f$logL,
      LR_vs_eq2 = if (ddf > 0) LR else NA, df = if (ddf > 0) ddf else NA, p = if (ddf > 0) pchisq(LR, ddf, lower.tail = FALSE) else NA)
  }
}
r <- do.call(rbind, out); rownames(r) <- NULL; write.csv(r, "output/table_WA19.csv", row.names = FALSE); print(r, digits = 4)
