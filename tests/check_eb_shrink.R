# Rscript tests/check_eb_shrink.R  (from the package root; needs ebnm; ~3 min)
suppressMessages({ library(gamlss2); library(gamlss.dist); library(future); library(doFuture) })
source("R/Utilities.R"); source("R/Core.R"); source("R/Support.R")
set.seed(4)

# mask 12x12x6 with a few holes; voxel order = which(mask) as in images2matrix()
mask <- array(TRUE, c(12, 12, 6)); mask[3, 3, ] <- FALSE; mask[10, 8, 2:4] <- FALSE
xyz <- which(mask, arr.ind = TRUE); nvox <- nrow(xyz)
truth <- 1 + 0.03 * xyz[, 2] + 0.6 * (xyz[, 1] > 6)       # smooth gradient + sharp step at x = 6|7
n <- 80
covs <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)))
img <- sapply(seq_len(nvox), function(v) rNO(n, truth[v] + 0.1 * (covs$sex == "M"), 0.3))
noisy <- sample(nvox, round(nvox / 3)); img[-(1:10), noisy] <- NA  # a third of voxels: 10 subjects only
img[, 5] <- NA                                                       # one voxel that cannot be fit
wd <- file.path(tempdir(), "eb"); dir.create(wd)
m <- vbgamlss(img, "Y ~ sex | 1", covs, g.family = NO, num_cores = 2, cachedir = wd, show_progress = FALSE,
              future_plan_strategy = "multisession", force_constraints = NULL)
stopifnot(m[[5]]$status == "too_few_obs")
m_bytes <- unclass(m)

b_raw <- vapply(seq_len(nvox), function(i) if (i == 5) NA else m[[i]]$coefficients$mu[["(Intercept)"]], 0)
err  <- function(b, idx) mean(abs(b[idx] - truth[idx]), na.rm = TRUE)
edge <- which(xyz[, 1] %in% 6:7)
shrunk <- list()
for (pf in c("normal", "point_laplace")) {
  s <- eb_shrink(m, mask, prior_family = pf, kernel_size = 5, num_cores = 2)
  b_eb <- vapply(seq_len(nvox), function(i) if (i == 5) NA else s[[i]]$coefficients$mu[["(Intercept)"]], 0)
  md   <- vapply(seq_len(nvox), function(i) if (is.null(s[[i]]$eb_shrink)) NA else s[[i]]$eb_shrink$prior_mode[["mu.(Intercept)"]], 0)
  lam  <- 1 - (b_eb - md) / (b_raw - md) # fraction shrunk towards the local mode
  cat(sprintf("%-14s |error| all %.3f -> %.3f | noisy %.3f -> %.3f | precise %.3f -> %.3f | step edge %.3f -> %.3f | median shrink noisy %.2f, precise %.2f\n",
              pf, err(b_raw, -5), err(b_eb, -5), err(b_raw, noisy), err(b_eb, noisy), err(b_raw, -c(noisy, 5)), err(b_eb, -c(noisy, 5)),
              err(b_raw, edge), err(b_eb, edge), median(lam[noisy], na.rm = TRUE), median(lam[-noisy], na.rm = TRUE)))
  stopifnot(err(b_eb, -5) < err(b_raw, -5), err(b_eb, noisy) < 0.7 * err(b_raw, noisy),
            median(lam[noisy], na.rm = TRUE) > median(lam[-noisy], na.rm = TRUE), # noisy voxels move more
            err(b_eb, edge) <= 1.1 * err(b_raw, edge))                           # step not blurred away
  shrunk[[pf]] <- list(s = s, b_eb = b_eb)
}
s <- shrunk$point_laplace$s; b_eb <- shrunk$point_laplace$b_eb

# original untouched; shrunk copy predicts with the shrunk coefficients; failed voxel left alone
stopifnot(identical(unclass(m), m_bytes), inherits(s, "vbgamlss"), isTRUE(s[[5]]$error))
i <- noisy[1]; nd <- data.frame(sex = factor(c("F", "M"), levels = c("F", "M")))
pm <- predict(s[[i]], newdata = nd, type = "parameter")$mu
stopifnot(isTRUE(all.equal(pm[1], b_eb[i])), isTRUE(all.equal(s[[i]]$eb_shrink$original[["mu.(Intercept)"]], b_raw[i])),
          s[[i]]$eb_shrink$posterior_sd[["mu.(Intercept)"]] > 0)

# sigma only: mu untouched, sigma changed
ss <- eb_shrink(m, mask, params = "sigma", num_cores = 2)
stopifnot(identical(ss[[i]]$coefficients$mu, m[[i]]$coefficients$mu),
          !identical(ss[[i]]$coefficients$sigma, m[[i]]$coefficients$sigma))

# wrong mask is caught
e <- tryCatch(eb_shrink(m, mask[, , 1:5]), error = function(e) conditionMessage(e))
stopifnot(grepl("voxels, the model has", e))
cat("eb_shrink checks passed\n")
