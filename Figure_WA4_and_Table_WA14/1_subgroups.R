# Figure WA4 and Table WA14, step 1: subgroup estimates for behavior: uncorrected random effects and Andrews and Kasy
# equation 2 (t, cutoff 1.645); the three-category version when a subgroup has fewer than five significant negative estimates.
# Run from this folder: Rscript 1_subgroups.R   (about 15-30 minutes)
source("ak_functions.R"); source("re_functions.R")
D <- read.csv("data/C_full.csv")
domain <- "behavior"
groups <- if (domain == "behavior") list(
  "Outcome: borrowing" = quote(otype == "B"), "Outcome: budgeting and planning" = quote(otype == "C"), "Outcome: saving and investing" = quote(otype == "D"),
  "Outcome: other behaviors and downstream outcomes" = quote(!otype %in% c("B", "C", "D")),
  "Intervention: classroom" = quote(interv == "Classroom"), "Intervention: counseling" = quote(interv == "Counseling"),
  "Intervention: educational nudge" = quote(interv == "Nudge"), "Intervention: online" = quote(interv == "Online"),
  "Sample: low-income" = quote(lowinc == "Low-income"), "Sample: general population" = quote(lowinc == "General population"),
  "Journal: top economics" = quote(topj == "Top economics journal"), "Journal: other" = quote(topj == "Other outlet"), "All" = quote(TRUE)) else list(
  "Intervention: classroom" = quote(interv == "Classroom"), "Intervention: online" = quote(interv == "Online"),
  "Sample: low-income" = quote(lowinc == "Low-income"), "Sample: general population" = quote(lowinc == "General population"),
  "Journal: top economics" = quote(topj == "Top economics journal"), "Journal: other" = quote(topj == "Other outlet"), "All" = quote(TRUE))
out <- list()
for (gn in names(groups)) {
  d <- D[eval(groups[[gn]], D), ]; z <- d$yi / d$sei; nneg <- sum(z < -1.645)
  spec <- if (nneg >= 5) "asym1645" else "asym1645_3"; sp <- if (spec == "asym1645") specs$asym1645 else specs3$asym1645_3
  f <- try(fit_one(d$yi, d$sei, as.character(d$cluster), sp, "t"), silent = TRUE); r <- try(re_model(d), silent = TRUE)
  ok <- !inherits(f, "try-error") && f$conv == 0
  out[[gn]] <- data.frame(domain = domain, group = gn, k = nrow(d), trials = length(unique(d$cluster)), n_signeg = nneg,
    re_mu = if (inherits(r, "try-error")) NA else r$mu, re_se = if (inherits(r, "try-error")) NA else r$se_CR2,
    re_lo = if (inherits(r, "try-error")) NA else r$lo, re_hi = if (inherits(r, "try-error")) NA else r$hi,
    ak_spec = spec, ak_mu = if (ok) f$P[1] else NA, ak_se = if (ok) f$se[1] else NA,
    ak_betas = if (ok) paste(sprintf("%.3f", f$betas), collapse = "; ") else NA,
    ak_se_betas = if (ok) paste(sprintf("%.3f", f$se_betas), collapse = "; ") else NA, converged = ok)
  print(out[[gn]][, c("group", "k", "n_signeg", "re_mu", "ak_spec", "ak_mu", "ak_se")], digits = 3)
}
write.csv(do.call(rbind, out), "output/table_WA14.csv", row.names = FALSE)
