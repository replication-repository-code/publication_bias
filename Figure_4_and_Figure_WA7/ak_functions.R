# Andrews and Kasy (2019) step-function selection models with t-distributed or normally distributed true effects.
# t model (used in the body of the paper): Andrews and Kasy's t specification. True effects are t-distributed,
#   theta ~ mu + tau * t(df), and each estimate is normal around its true effect, X | theta ~ N(theta, SE^2)
#   (ak_t.R). It nests the normal model (df -> infinity).
# Normal model: true effects ~ N(mu, tau^2) and estimate | true effect ~ N(true effect, SE^2), so the
#   latent (pre-selection) distribution of each estimate is N(mu, SE^2 + tau^2) (likelihood code in ak_likelihood_KH.R).
# Uses ak_likelihood_KH.R and ak_t.R (must be in the same folder).
# Each model is estimated from several starting values; the maximum of the likelihood is kept.
# Standard errors are clustered by trial. Fits are cached in output/cache/ (keyed by the data and the specification) and
# reused when a script is run again; delete that folder to refit everything.

source("ak_likelihood_KH.R"); source("ak_t.R")

# Selection specifications. cu = cutoffs on Z; sym = symmetric (equation 1) or asymmetric (equation 2).
specs <- list(
  none      = list(cu = c(1.96), sym = TRUE, sel = FALSE),
  sym1645   = list(cu = c(1.645), sym = TRUE, sel = TRUE),
  asym1645  = list(cu = c(-1.645, 0, 1.645), sym = FALSE, sel = TRUE),
  sym196    = list(cu = c(1.96), sym = TRUE, sel = TRUE),
  asym196   = list(cu = c(-1.96, 0, 1.96), sym = FALSE, sel = TRUE),
  symboth   = list(cu = c(1.645, 1.96), sym = TRUE, sel = TRUE),
  asymboth  = list(cu = c(-1.96, -1.645, 0, 1.645, 1.96), sym = FALSE, sel = TRUE))
# Three-category asymmetric versions (all negative estimates pooled into one region), used when
# very few estimates are significantly negative (subgroups; knowledge in web appendix I)
specs3 <- list(
  asym1645_3 = list(cu = c(0, 1.645), sym = FALSE, sel = TRUE),
  asym196_3  = list(cu = c(0, 1.96), sym = FALSE, sel = TRUE),
  asymboth_3 = list(cu = c(0, 1.645, 1.96), sym = FALSE, sel = TRUE))
spec_label <- c(none = "No selection", sym1645 = "Eq. 1 symmetric, 1.645", asym1645 = "Eq. 2 asymmetric, 1.645",
                sym196 = "Eq. 1 symmetric, 1.96", asym196 = "Eq. 2 asymmetric, 1.96",
                symboth = "Eq. 1 symmetric, both", asymboth = "Eq. 2 asymmetric, both")

# Cache: the value of expr is stored in output/cache/<tag>_<md5 of args>.rds and read from there if present.
ak_cached <- function(tag, args, expr) {
  dir.create("output/cache", showWarnings = FALSE, recursive = TRUE)
  tf <- tempfile(); saveRDS(args, tf, version = 2, compress = FALSE); key <- unname(tools::md5sum(tf)); unlink(tf)
  f <- file.path("output/cache", paste0(tag, "_", key, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  v <- expr; tmp <- paste0(f, ".tmp", Sys.getpid()); saveRDS(v, tmp); file.rename(tmp, f); v
}

# Normal model (Andrews and Kasy's likelihood with normal true effects, ak_likelihood_KH.R).
fit_normal <- function(X, s, cl, sp) {
  cu <- sp$cu; sym <- sp$sym; Tp <- Tpowers_fun(X / s, cu, sym); nb <- if (sp$sel) length(cu) else 0
  full <- function(P) {           # P = (mu, tau, publication probabilities)
    beta <- if (sp$sel) c(P[(length(P) - nb + 1):length(P)], 1) else rep(1, length(cu) + 1)
    VariationVarianceLogLikelihood(P[1], P[2], beta, cu, sym, X, s, Tp, df = Inf)
  }
  obj <- function(P) { v <- full(P)$LLH; if (!is.finite(v)) 1e10 else v }
  starts <- lapply(list(c(0, 1), c(0, 0.1), c(0.1, 0.3), c(mean(X), sd(X))), function(st) c(st, rep(1, nb)))
  best <- NULL
  for (st in starts) {
    o <- try(nlminb(st, obj, lower = c(-Inf, rep(0, length(st) - 1)), upper = rep(Inf, length(st)),
                    control = list(eval.max = MaxEval, iter.max = MaxIter, abs.tol = Tol)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$objective < best$objective)) best <- o
  }
  P <- best$par
  V <- try(RobustVariance(stepsize, length(X), P, full, cluster_ID = cl), silent = TRUE)
  se <- if (inherits(V, "try-error")) rep(NA, length(P)) else sqrt(pmax(diag(V), 0))
  list(P = P, se = se, logL = -best$objective, conv = best$convergence, pw = as.vector(full(P)$logL), npar = length(P),
       betas = if (nb) P[(length(P) - nb + 1):length(P)] else numeric(0), se_betas = if (nb) se[(length(P) - nb + 1):length(P)] else numeric(0))
}

# Fit one model. X = estimates, s = standard errors, cl = trial identifiers, sp = element of specs, model = "t" or "normal".
# The t model is also started from the solution of the normal model.
fit_one <- function(X, s, cl, sp, model) {
  ak_cached(paste0("ak_", model), list(X = X, s = s, cl = cl, cu = sp$cu, sym = sp$sym, sel = sp$sel, model = model), {
    if (model == "normal") fit_normal(X, s, cl, sp) else {
      nf <- fit_normal(X, s, cl, sp); nb <- if (sp$sel) length(sp$cu) else 0
      map <- if (sp$sel) c(seq_along(sp$cu), 0) else rep(0, length(sp$cu) + 1)
      ex <- list(c(nf$P[1:2], 3, nf$betas), c(nf$P[1], nf$P[2] / 2, 1.5, nf$betas))
      fit_t(X, s, cl, sp$cu, sp$sym, map, starts_extra = ex)
    }
  })
}

# Nesting of selection specifications (a nested in b)
nest_spec <- function(a, b) a == b || a == "none" || (b == "asymboth") ||
  (a == "sym1645" && b %in% c("symboth", "asym1645")) || (a == "sym196" && b %in% c("symboth", "asym196")) ||
  (a == "symboth" && b == "asymboth")
nested <- function(a, b) nest_spec(a$spec, b$spec) && (a$model == b$model || (a$model == "normal" && b$model == "t")) && !(a$spec == b$spec && a$model == b$model)
lr_p <- function(LR, ddf, boundary) {
  if (ddf <= 0) return(NA)
  if (!boundary) return(pchisq(max(LR, 0), ddf, lower.tail = FALSE))
  0.5 * (if (ddf - 1 > 0) pchisq(max(LR, 0), ddf - 1, lower.tail = FALSE) else as.numeric(LR <= 0)) + 0.5 * pchisq(max(LR, 0), ddf, lower.tail = FALSE)
}

# Fit a set of models to one data set and return one row per model, with likelihood-ratio tests,
# a clustered Vuong test against the minimum-AIC model, and LR tests of t versus normal. The t-versus-normal p-value uses a
# 50:50 mixture of chi2(0) and chi2(1), because df = infinity (the normal model) lies on the boundary of the parameter space.
ak_table <- function(d, models = c("normal", "t"), spec_names = names(specs), sample_name = "") {
  X <- d$yi; s <- d$sei; cl <- as.character(d$cluster); n <- length(X); fits <- list()
  for (m in models) for (k in spec_names) {
    sp <- if (k %in% names(specs)) specs[[k]] else specs3[[k]]
    f <- fit_one(X, s, cl, sp, m); f$spec <- k; f$model <- m; fits[[paste(m, k)]] <- f
    cat(sprintf("%-8s %-11s mu = %7.4f (SE %.4f)  logL = %.2f\n", m, k, f$P[1], f$se[1], f$logL))
  }
  tab <- do.call(rbind, lapply(fits, function(f) data.frame(
    sample = sample_name, k = n, clusters = length(unique(cl)), model = f$model, spec = f$spec,
    mu = f$P[1], se_mu = f$se[1], tau = f$P[2], se_tau = f$se[2], df_t = if (f$model == "t") f$P[3] else NA,
    betas = paste(sprintf("%.3f", f$betas), collapse = "; "), se_betas = paste(sprintf("%.3f", f$se_betas), collapse = "; "),
    logL = f$logL, npar = f$npar, AIC = -2 * f$logL + 2 * f$npar, conv = f$conv)))
  rownames(tab) <- NULL
  b <- which.min(tab$AIC); fb <- fits[[b]]; tab$best <- seq_len(nrow(tab)) == b; tab$dAIC <- tab$AIC - tab$AIC[b]
  tab$vuong_z <- NA; tab$vuong_p <- NA
  for (i in seq_along(fits)) if (i != b) {
    dd <- tapply(fb$pw - fits[[i]]$pw, cl, sum); z <- mean(dd) / (sd(dd) / sqrt(length(dd)))
    tab$vuong_z[i] <- z; tab$vuong_p[i] <- pnorm(z, lower.tail = FALSE)
  }
  tab$LR_t_vs_normal <- NA; tab$p_t_vs_normal <- NA
  if (all(c("t", "normal") %in% models)) for (i in seq_along(fits)) {
    f <- fits[[i]]; if (f$model == "t") {
      tab$LR_t_vs_normal[i] <- 2 * (f$logL - fits[[paste("normal", f$spec)]]$logL)
      tab$p_t_vs_normal[i] <- lr_p(tab$LR_t_vs_normal[i], 1, boundary = TRUE) }
  }
  tab
}

# Likelihood-ratio tests of equation 2 versus equation 1 at the same cutoffs (chi-square, df = difference in the number of
# publication probabilities; equation 1 imposes p1 = 1 and p2 = p3 on equation 2, interior restrictions)
lr_eq2_vs_eq1 <- function(tab, model = "t") {
  g <- function(s, v) tab[[v]][tab$model == model & tab$spec == s]
  do.call(rbind, lapply(list(c("1.645", "asym1645", "sym1645"), c("1.96", "asym196", "sym196"), c("both", "asymboth", "symboth")), function(p) {
    if (!all(c(p[2], p[3]) %in% tab$spec[tab$model == model])) return(NULL)
    LR <- 2 * (g(p[2], "logL") - g(p[3], "logL")); ddf <- g(p[2], "npar") - g(p[3], "npar")
    data.frame(model = model, cutoffs = p[1], LR_eq2_vs_eq1 = LR, df = ddf, p = lr_p(LR, ddf, boundary = FALSE)) }))
}
