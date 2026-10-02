# Rscript tests/check_ants.R  (from the package root; needs ANTsR + ebnm; ~2 min)
suppressMessages({ library(ANTsR); library(gamlss2); library(gamlss.dist); library(future); library(doFuture) })
for (f in c("R/Utilities.R", "R/Core.R", "R/Support.R", "R/Mapping.R")) source(f)
set.seed(5)
td <- file.path(tempdir(), "ants"); dir.create(td, showWarnings = FALSE)

# synthetic 10x10x6 grid, 2 mm, non-zero origin; mask with holes; known pattern per voxel
dims <- c(10, 10, 6)
mask_arr <- array(1, dims); mask_arr[c(1, 10), , ] <- 0; mask_arr[5, 5, 2:4] <- 0
ref <- function(a) { img <- as.antsImage(a); antsSetSpacing(img, c(2, 2, 2)); antsSetOrigin(img, c(-10, 5, 3)); img }
mask_path <- file.path(td, "mask.nii.gz"); antsImageWrite(ref(mask_arr), mask_path)
xyz <- which(mask_arr > 0, arr.ind = TRUE); nvox <- nrow(xyz)
pattern <- array(0, dims); pattern[mask_arr > 0] <- 1 + 0.05 * xyz[, 1] + 0.02 * xyz[, 2] + 0.01 * xyz[, 3]
n <- 40; covs <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)))
subj_imgs <- lapply(seq_len(n), function(s) {
  a <- pattern + 0.1 * (covs$sex[s] == "M") + array(rnorm(prod(dims), 0, 0.05), dims)
  a[xyz[1, , drop = FALSE]] <- 0 # first mask voxel never valid (constraints exclude 0) -> failed voxel 1; NIfTI can't keep NaN
  a * mask_arr
})
paths <- vapply(seq_len(n), function(s) { p <- file.path(td, sprintf("s%02d.nii.gz", s)); antsImageWrite(ref(subj_imgs[[s]]), p); p }, "")
img4d <- as.antsImage(array(unlist(subj_imgs), c(dims, n)))
antsSetSpacing(img4d, c(2, 2, 2, 1)); antsSetOrigin(img4d, c(-10, 5, 3, 0))
path4d <- file.path(td, "all4d.nii.gz"); antsImageWrite(img4d, path4d)

# images2matrix: both input paths agree, voxel order = which(mask) (column-major)
X_list <- as.matrix(images2matrix(as.list(paths), mask_path))
X_4d   <- as.matrix(images2matrix(path4d, mask_path))
stopifnot(dim(X_list) == c(n, nvox), isTRUE(all.equal(X_list, X_4d, check.attributes = FALSE)))
tol <- 1e-6 # NIfTI stores float32
expect <- t(vapply(subj_imgs, function(a) a[mask_arr > 0], numeric(nvox)))
stopifnot(isTRUE(all.equal(X_list, expect, check.attributes = FALSE, tolerance = tol)))
cat("images2matrix: list and 4D paths agree, voxel order matches the mask\n")

# fit (voxel 1 fails: all NaN)
m <- vbgamlss(X_list, "Y ~ sex | 1", covs, g.family = NO, num_cores = 2, cachedir = td, show_progress = FALSE,
              future_plan_strategy = "multisession", force_constraints = c(1e-8, Inf))
stopifnot(m[[1]]$status == "too_few_obs")

# coefficient maps land on the right voxels
fn <- map_model_coefficients(m, mask_path, file.path(td, "coef"), return_files = TRUE)
mu_map <- as.array(antsImageRead(grep("par-MU_coef-Intercept", fn, value = TRUE)))
b <- vapply(2:nvox, function(i) m[[i]]$coefficients$mu[["(Intercept)"]], 0)
stopifnot(isTRUE(all.equal(mu_map[mask_arr > 0][-1], b, tolerance = tol)), is.na(mu_map[mask_arr > 0][1]) || mu_map[mask_arr > 0][1] == 0)
cat("map_model_coefficients:", length(fn), "maps, values on the right voxels, failed voxel empty\n")

# prediction and z-score maps
p <- predict(m, newdata = covs, ptype = "parameter", num_cores = 2)
fp <- map_model_predictions(p, mask_path, file.path(td, "pred"), return_files = TRUE)
fp2 <- map_model_predictions(p, mask_path, file.path(td, "predsub"), index = c(2, 5), return_files = TRUE)
mu3 <- as.array(antsImageRead(grep("subj-3_.*par-MU", fp, value = TRUE)))
stopifnot(isTRUE(all.equal(mu3[mask_arr > 0][-1], vapply(2:nvox, function(i) p[[i]]$mu[3], 0), tolerance = tol)))
stopifnot(length(fp2) == 2 * 2, all(grepl("subj-(2|5)_", fp2)), length(fp) == n * 2) # 2 params (mu, sigma)
cat("map_model_predictions: all subjects and index subset ok\n")

z <- zscore.vbgamlss(p, X_list, num_cores = 2)
fz <- map_zscores(z, mask_path, file.path(td, "z"), return_files = TRUE)
z4 <- as.array(antsImageRead(fz, 4))
stopifnot(dim(z4)[4] == n, isTRUE(all.equal(z4[, , , 7][mask_arr > 0][-1], vapply(2:nvox, function(i) z[[i]][7], 0), tolerance = tol)))
fz3 <- map_zscores(z, mask_path, file.path(td, "z3"), index = c(1, 4), output_4D = FALSE, return_files = TRUE)
stopifnot(length(fz3) == 2)
cat("map_zscores: 4D and 3D subset ok\n")

# eb_shrink: mask as a path gives the same result as the array
if (requireNamespace("ebnm", quietly = TRUE)) {
  s_path <- eb_shrink(m, mask_path, num_cores = 2, min_neighbors = 5)
  s_arr  <- eb_shrink(m, mask_arr,  num_cores = 2, min_neighbors = 5)
  stopifnot(identical(s_path[[10]]$coefficients, s_arr[[10]]$coefficients), !is.null(s_path[[10]]$eb_shrink))
  cat("eb_shrink: mask path == mask array\n")
}
cat("ANTsR checks passed\n")
