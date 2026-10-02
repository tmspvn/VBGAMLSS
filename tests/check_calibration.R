# Rscript tests/check_calibration.R  (from the package root)
suppressMessages({ library(gamlss2); library(gamlss.dist) })
source("R/Utilities.R"); source("R/Lobo_cross_validation.R")
set.seed(1)

# formula parsing + single-level guard
f <- "Y ~ sex + pb(age) + s(protocol, bs='re') + s(scanner_hcp,bs = \"re\") | 1 + s(protocol,bs='re') | 1 | 1"
re <- re_vars_by_param(f)
stopifnot(identical(re$mu, c("protocol", "scanner_hcp")), identical(re$sigma, "protocol"), length(re$nu) == 0)
g <- capture.output(f2 <- drop_degenerate_re(f, data.frame(protocol = c("a", "b"), scanner_hcp = "pooled")))
stopifnot(!grepl("scanner_hcp", f2), grepl("s\\(protocol", f2))
check_formula_LHS(f2)

# held-out batch: 2 new scanners with known link-scale offsets on top of the population model.
# Families differ in links / mean availability to check nothing is distribution-specific.
recover <- function(family, rdist, par0, dmu, dsig, lims) {
  fam <- gamlss2:::complete_family(family)
  n <- 400; grp <- factor(rep(c("s1", "s2"), each = n / 2))
  lf <- function(p) gamlss2:::make.link2(fam$links[[p]])$linkfun
  li <- function(p, eta) fam$map2par(setNames(list(eta), p))[[p]]
  par <- lapply(par0, function(v) rep(v, n))
  truth <- par
  truth$mu    <- li("mu", lf("mu")(par$mu) + dmu[grp])
  truth$sigma <- li("sigma", lf("sigma")(par$sigma) + dsig[grp])
  y <- do.call(rdist, c(list(n), truth))
  y[c(5, 300)] <- NaN # masked values must not break it
  nfit <- c(par, list(y = y, yhat = rep(0, n)))
  folds <- sample(rep_len(1:5, n))
  cal <- calibrate_heldout(nfit, folds, list(mu = grp, sigma = grp), fam)
  est_dmu  <- tapply(lf("mu")(cal$mu) - lf("mu")(par$mu), grp, mean)
  est_dsig <- tapply(lf("sigma")(cal$sigma) - lf("sigma")(par$sigma), grp, mean)
  ok <- is.finite(y)
  nll <- function(m) -sum(fam$pdf(y[ok], lapply(m[fam$names], function(v) v[ok]), log = TRUE))
  cat(sprintf("%-7s links %-28s dmu %s -> %s | dsig %s -> %s | NLL %.0f -> %.0f\n", fam$family,
              paste(fam$links, collapse = "/"), paste(dmu, collapse = ","), paste(round(est_dmu, 2), collapse = ","),
              paste(round(dsig, 2), collapse = ","), paste(round(est_dsig, 2), collapse = ","), nll(nfit), nll(cal)))
  stopifnot(abs(est_dmu - dmu) < lims[1], abs(est_dsig - dsig) < lims[2], nll(cal) < nll(nfit), all(is.finite(cal$yhat)))
  cal_mu <- calibrate_heldout(nfit, folds, list(mu = grp), fam) # mu-only RE: sigma untouched
  stopifnot(identical(cal_mu$sigma, par$sigma))
}
recover(SHASHo, rSHASHo, list(mu = 1, sigma = 0.2, nu = 0.3, tau = 1.2), c(s1 = .5, s2 = -.3), log(c(s1 = 1.5, s2 = .7)), c(.06, .15))
recover(BCPEo, rBCPEo, list(mu = 5, sigma = 0.15, nu = 1, tau = 2), c(s1 = .2, s2 = -.1), log(c(s1 = 1.4, s2 = .8)), c(.04, .15)) # no $mean
recover(BE, rBE, list(mu = 0.4, sigma = 0.3), c(s1 = .6, s2 = -.4), c(s1 = .4, s2 = -.3), c(.12, .15))                 # logit links
cat("calibration checks passed\n")
