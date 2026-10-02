# Rscript tests/check_core.R  (from the package root; a few minutes, 2 local workers)
suppressMessages({ library(gamlss2); library(gamlss.dist); library(future); library(doFuture) })
source("R/Utilities.R"); source("R/Core.R")
set.seed(1)

# TRY: a fit that warns then errors runs ONCE and the error is captured
runs <- 0
r <- TRY({ runs <- runs + 1; warning("w1"); stop("boom") })
stopifnot(runs == 1, identical(r$error, "boom"), identical(r$warnings, "w1"), identical(r$value, NA))
runs <- 0
r <- TRY({ runs <- runs + 1; warning("w2"); 42 })
stopifnot(runs == 1, r$value == 42, is.null(r$error))

# synthetic data: 12 voxels, voxel 10 has only 3 valid observations, voxel 11 makes gamlss2 error
n <- 150
covs <- data.frame(age = runif(n, 10, 80), site = factor(sample(c("a", "b", "c"), n, TRUE)),
                   junk = matrix(rnorm(n * 40), n))  # unused columns must not reach workers
img <- sapply(1:12, function(v) rSHASHo(n, 1 + 0.01 * covs$age + c(a = 0, b = .3, c = -.2)[covs$site], 0.3, 0.1, 1))
img[sample(n, 10), 3] <- NA
img[-(1:3), 10] <- NA
img[1, 11] <- Inf
cache_root <- file.path(tempdir(), "vbg_check"); dir.create(cache_root)
fo <- "Y ~ age + site | 1 | 1 | 1"  # kept simple: a pb() + re SHASH fit takes ~40 s here
run <- function(cache_id = "outer1", ...) vbgamlss(img, fo, covs, g.family = SHASHo, num_cores = 2, cachedir = cache_root,
                                                cache_id = cache_id, debug = TRUE, show_progress = FALSE,
                                                future_plan_strategy = "multisession", force_constraints = c(-Inf, Inf), ...)

m1 <- run()
cdir <- attr(m1, "cachedir")
reg_path <- file.path(cdir, ".vbgamlss.registry")
reg <- qs2::qs_read(reg_path)
print(table(reg$status))
stopifnot(all(reg$fitted), reg$status[10] == "too_few_obs", reg$status[11] == "error",
          grepl("Inf", reg$message[11]), all(reg$status[-(10:11)] %in% c("ok", "warned", "nonconverged")))
stopifnot(isTRUE(m1[[10]]$error), m1[[10]]$vxl == 10)

# shards: one file per worker, no per-voxel files
shards <- list.files(file.path(cdir, ".voxfits"))
cat("shard files:", length(shards), "\n")
stopifnot(length(shards) >= 1, length(shards) <= 2, all(grepl("^shard\\..*\\.bin$", shards)))

# failed voxel keeps everything its fit had and replays alone; clean fits keep nothing extra
dc <- m1[[11]]$debug_call
stopifnot(nrow(dc$data) == n, identical(dc$maxit, c(100, 33)), dc$control$eps == 1e-5, is.null(m1[[1]]$debug_call))
err <- tryCatch(refit_voxel(m1, 11), error = function(e) conditionMessage(e))
stopifnot(identical(err, "NA/NaN/Inf in 'y'"))
fixed <- refit_voxel(m1, 11, data = transform(dc$data, Y = ifelse(is.finite(Y), Y, NA)), control = list(eps = 1e-4))
stopifnot(inherits(fixed, "gamlss2"))

# saved voxel models are small and still predict
cat(sprintf("voxel models: max %.1f KB\n", max(reg$nbytes) / 2^10))
stopifnot(max(reg$nbytes) < 2^20)
p <- predict(m1[[1]], newdata = covs[1:3, ], type = "parameter")$mu
stopifnot(length(p) == 3, all(is.finite(p)))

# simulated kill: registry missed voxels 5-12, and a shard ends in a half-written record
reg$fitted[5:12] <- FALSE; reg$shard[5:12] <- NA; qs2::qs_save(reg, reg_path)
con <- file(file.path(cdir, ".voxfits", shards[1]), "ab")
writeBin(c(SHARD_MAGIC, 99L, 5000L), con, size = 4); writeBin(as.raw(1:10), con); close(con)

t0 <- Sys.time(); out <- capture.output(m2 <- run()); dt <- as.numeric(Sys.time() - t0, units = "secs")
stopifnot(any(grepl("Resuming from", out)), any(grepl("Recovered 8 voxels", out)), any(grepl("already fitted", out)),
          identical(attr(m2, "cachedir"), cdir), length(list.dirs(cache_root, recursive = FALSE)) == 1,
          identical(unclass(m1)[1:12], unclass(m2)[1:12]))
cat(sprintf("resume took %.1f s: recovered from shards, no refit, identical models\n", dt))

# retry_failed: only the error voxel is refit (and fails the same way); too_few_obs is not retried
out <- capture.output(m3 <- run(retry_failed = TRUE))
stopifnot(any(grepl("Processing  1  voxels", out)), m3[[11]]$status == "error")
cat("retry_failed: refit only the failed voxel\n")

# cache removed once the aggregated model is on disk
drop_vbgamlss_cache(m3)
stopifnot(!dir.exists(cdir))

# several chunks: workers keep appending to their own shard across chunks, same fits
out <- capture.output(mc <- run("chunks", chunk_max_mb = 0.004))
stopifnot(any(grepl("Chunk: 4/4", out)), length(list.files(file.path(attr(mc, "cachedir"), ".voxfits"))) <= 2,
          isTRUE(all.equal(coef(m1[[1]]), coef(mc[[1]]))), isTRUE(all.equal(coef(m1[[12]]), coef(mc[[12]]))))
cat("4 chunks: same fits, shards reused across chunks\n")

# segmentation: voxel 2 has only 3 subjects with the target label
seg <- matrix(1L, n, 12); seg[-(1:3), 2] <- 2L
ms <- run("seg", segmentation = seg, segmentation_target = 1)
stopifnot(ms[[2]]$status == "too_few_obs", ms[[1]]$status == "ok")
cat("segmentation: labels respected\n")

# a cache from the old one-file-per-voxel layout is skipped via cache_id, refused when passed directly
old <- file.path(cache_root, "legacy.vbgamlss.cache.OLDx"); dir.create(old)
qs2::qs_save(data.frame(voxel = 1:12, fitted = TRUE, converged = TRUE, relative_paths = "x", full_paths = "x"),
             file.path(old, ".vbgamlss.registry"))
out <- capture.output(ml <- run("legacy"))
stopifnot(any(grepl("Ignoring old-format cache", out)), attr(ml, "cachedir") != old)
err <- tryCatch(vbgamlss(img, fo, covs, g.family = SHASHo, cachedir = old, num_cores = 2, show_progress = FALSE,
                         future_plan_strategy = "multisession"), error = function(e) conditionMessage(e))
stopifnot(grepl("older vbgamlss", err))
cat("old-format cache: skipped / refused\n")
cat("core checks passed\n")
