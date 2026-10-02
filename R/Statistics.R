################################  Statistics  ##################################
# Summaries of a vbgamlss.model_selection_NCV() run. All metrics are lower-is-better.
# Per-fold value = the fold's describe_stats()$mean (the voxel-wise mean of that metric).


#' Per-fold metrics of every candidate of a nested-CV model selection, as one tidy table
#'
#' @param results A result file of vbgamlss.model_selection_NCV() (path or loaded object).
#' @param registry Only for result files written before the job table was stored with the
#'   results: the run's `.Slurm.Registry` (path or object), which maps candidates to images.
#' @param metrics Metrics to extract.
#' @return data.frame, one row per candidate x CV level x fold x metric: `job`, `Y` (image label),
#'   `family`, `formula`, `label` (only the formula parts that differ between candidates, i.e. the
#'   greedy step being selected), `cv` ("LOBO" / "innerCV"), `outer`, `fold`, `metric`, `value`, and
#'   `missfit` (fraction of voxels that failed, LOBO folds only).
#' @export
ncv_folds <- function(results, registry = NULL, metrics = c("GD", "MAE", "LL", "CLL")) {
  if (is.character(results)) results <- qs2::qs_read(results)
  res  <- if (!is.null(results$results)) results$results else results
  jobs <- attr(res, "jobs")
  if (is.null(jobs)) {
    if (is.null(registry)) stop("this result file has no job table: pass the run's .Slurm.Registry as `registry`")
    if (is.character(registry)) registry <- qs2::qs_read(registry)
    jobs <- data.frame(formula = registry$formula, image = registry$image, family = registry$family)
  }
  jobs$Y <- image_label(jobs$image)
  jobs$label <- varying_parts(jobs$formula)

  fold_mean <- function(fold, m) {
    v <- fold[[m]]
    if (is.null(v)) NA_real_ else as.numeric(v[["mean"]])
  }
  one_cv <- function(folds, cv, outer) {
    d <- expand.grid(fold = seq_along(folds), metric = metrics, stringsAsFactors = FALSE)
    d$value <- mapply(function(f, m) fold_mean(folds[[f]], m), d$fold, d$metric)
    d$missfit <- if (cv == "LOBO") {
      vapply(d$fold, function(f) { v <- folds[[f]]$missFitsPerc; if (is.null(v)) NA_real_ else as.numeric(v) }, 0)
    } else NA_real_
    d$cv <- cv
    d$outer <- if (is.null(outer)) d$fold else outer
    d
  }
  rows <- lapply(seq_along(res), function(i) {
    fit <- res[[i]]
    if (!is.list(fit) || is.null(fit$outer)) return(NULL) # failed candidate
    d <- rbind(one_cv(fit$outer, "LOBO", NULL),
               do.call(rbind, lapply(seq_along(fit$inner), function(o) one_cv(fit$inner[[o]], "innerCV", o))))
    cbind(job = i, jobs[i, c("Y", "family", "formula", "label")], d, row.names = NULL)
  })
  out <- do.call(rbind, rows)
  out[, c("job", "Y", "family", "formula", "label", "cv", "outer", "fold", "metric", "value", "missfit")]
}


#' Per-candidate summary and ranking of a nested-CV model selection
#'
#' Within each image (`Y`): mean over folds of every metric under LOBO and inner CV, their ranks
#' (1 = best), a composite rank, the share of folds each candidate wins (GD), the batch-to-batch SD
#' of the LOBO GD, and the worst-fold missfit.
#' @param folds Output of ncv_folds().
#' @param weights Weights of the mean LOBO rank and the mean inner-CV rank in the composite.
#' @param metrics Metrics that enter the ranks.
#' @return data.frame sorted by `Y` and `Composite`, with `best` flagging the winner per image.
#' @export
ncv_summary <- function(folds, weights = c(LOBO = 0.5, innerCV = 0.5), metrics = c("GD", "MAE", "LL", "CLL")) {
  folds <- folds[folds$metric %in% metrics, ]
  key <- unique(folds[, c("job", "Y", "family", "formula", "label")])
  out <- key
  for (cv in c("LOBO", "innerCV")) for (m in metrics) {
    s <- folds[folds$cv == cv & folds$metric == m, ]
    out[[paste0(cv, "_", m)]] <- tapply(s$value, s$job, mean, na.rm = TRUE)[as.character(out$job)]
  }
  lobo_gd <- folds[folds$cv == "LOBO" & folds$metric == "GD", ]
  out$LOBO_GD_sd <- tapply(lobo_gd$value, lobo_gd$job, stats::sd, na.rm = TRUE)[as.character(out$job)]
  out$missFitPct <- 100 * tapply(lobo_gd$missfit, lobo_gd$job, function(v) suppressWarnings(max(v, na.rm = TRUE)))[as.character(out$job)]
  out$missFitPct[!is.finite(out$missFitPct)] <- NA

  win_pct <- function(s) { # share of folds where the candidate has the lowest value
    w <- tapply(seq_len(nrow(s)), paste(s$cv, s$outer, s$fold), function(r) s$job[r][which.min(s$value[r])])
    w <- unlist(w)
    100 * tabulate(match(w, out$job), nbins = nrow(out)) / max(length(w), 1)
  }
  out$LOBO_GD_WinPct <- out$innerCV_GD_WinPct <- NA_real_
  out$mean_rank_LOBO <- out$mean_rank_innerCV <- out$Composite <- NA_real_
  for (y in unique(out$Y)) {
    ix <- which(out$Y == y)
    for (cv in c("LOBO", "innerCV")) {
      s <- folds[folds$Y == y & folds$cv == cv & folds$metric == "GD", ]
      out[[paste0(cv, "_GD_WinPct")]][ix] <- win_pct(s)[ix]
      rk <- sapply(metrics, function(m) rank(out[[paste0(cv, "_", m)]][ix], ties.method = "min", na.last = "keep"))
      out[[paste0("mean_rank_", cv)]][ix] <- rowMeans(matrix(rk, nrow = length(ix)), na.rm = TRUE)
    }
  }
  out$Composite <- (weights[["LOBO"]] * out$mean_rank_LOBO + weights[["innerCV"]] * out$mean_rank_innerCV) / sum(weights)
  out <- out[order(out$Y, out$Composite), ]
  out$best <- !duplicated(out$Y)
  rownames(out) <- NULL
  out
}


#' Paired per-fold comparison of every candidate against a reference model
#'
#' On every fold, 100 * (metric - metric_ref) / |metric_ref|: batch difficulty is shared by all
#' candidates on a fold, so the paired change isolates the model difference. Negative = better.
#' @param folds Output of ncv_folds().
#' @param reference The reference model: its full formula or its `label` (whitespace-insensitive).
#' @param n_boot Bootstrap replicates for the 95% CI of the median change.
#' @param seed Seed for the bootstrap (the caller's RNG stream is left untouched).
#' @return data.frame, one row per image x metric x CV level x candidate: `median`, `lo`, `hi`
#'   (bootstrap 95% CI), `signif` (CI excludes 0), `n_folds`. Attribute "reference" holds the
#'   reference's own metric per image x metric x CV level: median, SD over folds, bootstrap SE.
#' @export
ncv_paired <- function(folds, reference, n_boot = 2000, seed = 0) {
  norm <- function(x) gsub("\\s|1\\+", "", x)
  is_ref <- norm(folds$formula) == norm(reference) | norm(folds$label) == norm(reference)
  if (!any(is_ref)) stop("reference not found among the candidates; labels are: ",
                         paste(unique(folds$label), collapse = " ; "))
  k <- c("Y", "cv", "outer", "fold", "metric")
  ref <- folds[is_ref, c(k, "value")]
  names(ref)[names(ref) == "value"] <- "ref"
  d <- merge(folds[!is_ref, ], ref, by = k)
  d$pct <- 100 * (d$value - d$ref) / abs(d$ref)

  boot_median <- function(v) { v <- v[is.finite(v)]; apply(matrix(sample(v, n_boot * length(v), TRUE), n_boot), 1, stats::median) }
  with_seed(seed, {
    g <- split(d$pct, d[, c("Y", "metric", "cv", "job")], drop = TRUE)
    st <- do.call(rbind, lapply(names(g), function(nm) {
      v <- g[[nm]]; b <- boot_median(v)
      data.frame(id = nm, median = stats::median(v[is.finite(v)]), lo = stats::quantile(b, .025, names = FALSE),
                 hi = stats::quantile(b, .975, names = FALSE), n_folds = sum(is.finite(v)))
    }))
    r <- split(ref$ref, ref[, c("Y", "metric", "cv")], drop = TRUE)
    ref_stats <- do.call(rbind, lapply(names(r), function(nm) {
      v <- r[[nm]][is.finite(r[[nm]])]
      data.frame(id = nm, median = stats::median(v), sd = stats::sd(v), se = stats::sd(boot_median(v)))
    }))
  })
  key <- unique(d[, c("Y", "metric", "cv", "job", "family", "formula", "label")])
  key$id <- interaction(key[, c("Y", "metric", "cv", "job")], drop = TRUE)
  out <- merge(key, st, by = "id")
  out$signif <- out$lo > 0 | out$hi < 0
  rk <- unique(folds[is_ref, c("Y", "metric", "cv", "label", "formula")])
  rk$id <- as.character(interaction(rk[, c("Y", "metric", "cv")], drop = TRUE))
  attr_ref <- merge(rk, ref_stats, by = "id")
  out$id <- NULL; attr_ref$id <- NULL
  out <- out[order(out$Y, out$metric, out$cv, out$median), ]
  rownames(out) <- NULL
  attr(out, "reference") <- attr_ref
  out
}


# image label from its path: ".../INPUT/MK_ready/mk_cropped.nii.gz" -> "MK", else the file name
image_label <- function(image) {
  lab <- ifelse(grepl("INPUT/[^/_]+", image), sub("^.*INPUT/([^/_]+).*$", "\\1", image),
                sub("\\.nii(\\.gz)?$", "", basename(image)))
  toupper(lab)
}


# the "|" parts that differ between candidates (= the greedy step being selected), without "1 +"
varying_parts <- function(formulas) {
  parts <- lapply(strsplit(sub("^[^~]*~", "", formulas), "|", fixed = TRUE), function(p) trimws(sub("^\\s*1\\s*\\+", "", p)))
  np <- max(lengths(parts))
  parts <- lapply(parts, function(p) { length(p) <- np; p })
  vary <- vapply(seq_len(np), function(j) length(unique(vapply(parts, `[`, "", j))) > 1, logical(1))
  if (!any(vary)) vary[1] <- TRUE
  vapply(parts, function(p) paste(p[vary], collapse = "  |  "), "")
}
