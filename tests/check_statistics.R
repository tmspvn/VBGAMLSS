# Rscript tests/check_statistics.R  (from the package root; needs ggplot2, scales, gt; ~10 s)
source("R/Lobo_cross_validation.R"); source("R/Statistics.R"); source("R/Plotting.R")
set.seed(6)

# synthetic NCV result: 2 images x 3 sigma candidates, 4 outer folds x 3 inner folds.
# fold difficulty is shared by all candidates (as in real batches); "s(site, bs='re')" is truly best.
mu <- "sex + pb(age) + s(site, bs='re')"
cands <- c("1", "s(site, bs='re')", "sex + s(site, bs='re')")
effect <- c(0, -0.08, -0.04)
jobs <- expand.grid(k = seq_along(cands), image = c("/x/INPUT/FA_ready/fa.nii.gz", "/x/INPUT/MD_ready/md.nii.gz"), stringsAsFactors = FALSE)
jobs$formula <- paste("Y ~", mu, "|", cands[jobs$k], "| 1 | 1")
difficulty <- list(outer = rnorm(4, -500, 150), inner = matrix(rnorm(12, -900, 50), 3))
fold <- function(gd, miss = 0) list(GD = list(mean = gd), MAE = list(mean = 0.1 + rnorm(1, 0, 1e-4)),
                                    LL = list(mean = gd / 1000), CLL = list(mean = gd / 900), missFitsPerc = miss)
res <- lapply(seq_len(nrow(jobs)), function(i) {
  e <- effect[jobs$k[i]]
  list(outer = lapply(1:4, function(o) fold(difficulty$outer[o] * (1 - e) + rnorm(1, 0, 2), miss = (i == 3) * 0.02)),
       inner = lapply(1:4, function(o) lapply(1:3, function(f) fold(difficulty$inner[f, o] * (1 - e) + rnorm(1, 0, 2)))))
})
res[[6]] <- NA # a failed candidate
attr(res, "jobs") <- data.frame(formula = jobs$formula, image = jobs$image, family = "SHASH")
results <- list(results = res, model_ranking = NULL)

f <- ncv_folds(results)
stopifnot(setequal(unique(f$Y), c("FA", "MD")), setequal(unique(f$label), cands),        # step part detected
          nrow(f) == 5 * (4 + 12) * 4, max(f$missfit, na.rm = TRUE) == 0.02)              # failed candidate dropped

s <- ncv_summary(f)
stopifnot(all(s$label[s$best] == "s(site, bs='re')"), s$missFitPct[s$job == 3] == 2,
          all(abs(tapply(s$LOBO_GD_WinPct, s$Y, sum) - 100) < 1e-9))                      # win % sums to 100
s_inner <- ncv_summary(f, weights = c(LOBO = 0, innerCV = 1))
stopifnot(isTRUE(all.equal(s_inner$Composite, s_inner$mean_rank_innerCV)))

p <- ncv_paired(f, "1")
r <- attr(p, "reference")
best <- p[p$label == "s(site, bs='re')" & p$metric == "GD", ]
stopifnot(all(best$median < -5), all(best$signif), nrow(r) == 2 * 4 * 2,                  # ~ -8 % GD, CI excludes 0
          identical(ncv_paired(f, "Y ~ sex + pb(age) + s(site, bs='re') | 1 | 1 | 1")$median, p$median))  # full formula works
e <- tryCatch(ncv_paired(f, "nope"), error = conditionMessage)
stopifnot(grepl("reference not found", e))

g <- plot_ncv_forest(p)
stopifnot(inherits(g, "ggplot"), !"MAE" %in% levels(g$data$metric))                         # sigma step: no MAE
td <- file.path(tempdir(), "ncvrep")
out <- ncv_report(results, td, name = "SIG")
stopifnot(all(file.exists(file.path(td, c("SIG_folds.csv", "SIG_summary.csv", "SIG_paired.csv", "SIG_forest.png",
                                          "SIG_FA_table.html", "SIG_MD_table.html")))))
cat("statistics checks passed\n")
