#
#
#
#
#
#
#
#
#
#


# Define the plotting function
plot_gamlss_shash <- function(mu = 0, sigma = 1, nu = 1, tau = 1, x_range = c(-5, 5)) {
  # Create a sequence of x values
  x_vals <- seq(x_range[1], x_range[2], length.out = 500)

  # Create a grid of all parameter combinations passed to the function
  params <- expand.grid(mu = mu, sigma = sigma, nu = nu, tau = tau)

  # Calculate the density for each combination of parameters using gamlss.dist
  df_list <- lapply(1:nrow(params), function(i) {
    m <- params$mu[i]
    s <- params$sigma[i]
    n <- params$nu[i]
    t <- params$tau[i]

    # Using dSHASH from the gamlss.dist package
    y_vals <- dSHASH(x = x_vals, mu = m, sigma = s, nu = n, tau = t)

    # Create a descriptive label for the legend
    label <- sprintf("μ=%g, σ=%g, ν=%g, τ=%g", m, s, n, t)

    data.frame(x = x_vals, y = y_vals, label = as.factor(label))
  })

  # Combine the lists into a single data frame for ggplot
  df <- do.call(rbind, df_list)

  # Generate the overlaid plots
  ggplot(df, aes(x = x, y = y, color = label)) +
    geom_line(linewidth = 1) +
    labs(
      title = "SHASH Distribution",
      subtitle = "Comparing density curves across varying parameters",
      x = "x",
      y = "Density",
      color = "Parameters"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "vertical")
}

# --- Examples of how to use the function ---

# 1. Varying the Scale (sigma)
# plot_gamlss_shash(sigma = c(0.5, 1, 2))

# 2. Varying the Skewness (nu)
# plot_gamlss_shash(nu = c(-1, 0, 1, 2))

# 3. Varying the Tailweight/Kurtosis (tau)
# plot_gamlss_shash(tau = c(0.5, 0.8, 1, 1.5), x_range = c(-8, 8))

# 4. Combining multiple variations at once
# plot_gamlss_shash(nu = c(-1, 1), tau = c(0.5, 1.5), x_range = c(-6, 6))



# ==============================================================================
# Nested-CV model selection: forest plot, summary table, one-call report
# ==============================================================================

#' Forest plot of a nested-CV step: every candidate against a reference model
#'
#' Rows = candidates (simplest at the top); x = paired per-fold % change vs the reference
#' (median and bootstrap 95% CI, pseudo-log axis, left of 0 = better); LOBO and inner CV side by
#' side; one row of panels per image. The reference row (shaded) shows its own metric:
#' median (SD over folds, bootstrap SE). MAE is shown only when the mu formula is among the parts
#' being compared (it measures the location fit).
#' @param paired Output of ncv_paired().
#' @param metrics Metrics to show, in column order.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @export
plot_ncv_forest <- function(paired, metrics = c("MAE", "GD", "CLL"), title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("scales", quietly = TRUE))
    stop("plot_ncv_forest needs ggplot2 and scales")
  ref <- attr(paired, "reference")
  mu_part <- function(fo) trimws(sub("^\\s*1\\s*\\+", "", strsplit(sub("^[^~]*~", "", fo), "|", fixed = TRUE)[[1]][1]))
  if (length(unique(vapply(c(paired$formula, ref$formula), mu_part, ""))) < 2) metrics <- setdiff(metrics, "MAE")
  metrics <- intersect(metrics, unique(paired$metric))
  d <- paired[paired$metric %in% metrics, ]
  ref <- ref[ref$metric %in% metrics, ]

  ref_label <- ref$label[1]
  ord <- complexity_order(unique(c(d$label, ref_label)))
  ylab <- gsub(",\\s*bs\\s*=\\s*['\"]re['\"]", "", ord)
  ylab[ord == ref_label] <- paste0(ylab[ord == ref_label], "  (reference)")
  off <- c(LOBO = -0.15, innerCV = 0.15)
  d$y   <- match(d$label, ord) + off[d$cv]
  ref$y <- match(ref$label, ord) + off[ref$cv]
  ref$txt <- sprintf("%s: %.3g  (SD %.2g, SE %.2g)", ref$cv, ref$median, ref$sd, ref$se)
  facet_lab <- c(MAE = "MAE", GD = "GD / LL", LL = "LL", CLL = "CLL")
  d$metric   <- factor(d$metric, metrics, facet_lab[metrics])
  ref$metric <- factor(ref$metric, metrics, facet_lab[metrics])
  cols <- c(LOBO = "#1F9BCF", innerCV = "#E8862E")
  k <- match(ref_label, ord)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$median, y = .data$y, colour = .data$cv)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = k - 0.45, ymax = k + 0.45, alpha = 0.08) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey45") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$signif), size = 2.2) +
    ggplot2::geom_point(data = ref, ggplot2::aes(x = 0), shape = 18, size = 3.2) +
    ggplot2::geom_text(data = ref, ggplot2::aes(x = 0, label = .data$txt), hjust = -0.08, size = 2.6, show.legend = FALSE) +
    ggplot2::facet_grid(Y ~ metric, scales = "free_x") +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5, base = 10),
                                breaks = c(-100, -10, -1, 0, 1, 10, 100)) +
    ggplot2::scale_y_reverse(breaks = seq_along(ord), labels = ylab, expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_colour_manual(values = cols, labels = c(LOBO = "LOBO (held-out batch)", innerCV = "inner CV"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "CI includes 0", `TRUE` = "CI excludes 0"), name = NULL) +
    ggplot2::labs(title = title, y = NULL,
                  x = "change vs the reference (shaded row) on the same fold, %   (median, bootstrap 95% CI;  < 0 = better)") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0, face = "bold"))
}


#' Formatted table of a nested-CV summary for one image
#'
#' @param summary Output of ncv_summary().
#' @param Y The image to show.
#' @param title Optional title.
#' @param file Optional .html path to save the table to.
#' @return A gt table (needs the gt package). Values of a diverged fit (more than 10x the column's
#'   median magnitude) are displayed as "Inf"; the underlying numbers are untouched.
#' @export
table_ncv_summary <- function(summary, Y, title = NULL, file = NULL) {
  if (!requireNamespace("gt", quietly = TRUE)) stop("table_ncv_summary needs the gt package")
  d <- summary[summary$Y == Y, ]
  cols <- c("family", "label", "Composite", "mean_rank_LOBO", "mean_rank_innerCV",
            "LOBO_GD", "LOBO_GD_sd", "LOBO_GD_WinPct", "LOBO_MAE", "LOBO_LL", "LOBO_CLL",
            "innerCV_GD", "innerCV_GD_WinPct", "innerCV_MAE", "innerCV_LL", "innerCV_CLL", "missFitPct")
  cols <- intersect(cols, names(d))
  metric_cols <- intersect(c("LOBO_GD", "LOBO_GD_sd", "LOBO_MAE", "LOBO_LL", "LOBO_CLL",
                             "innerCV_GD", "innerCV_MAE", "innerCV_LL", "innerCV_CLL"), cols)
  g <- gt::gt(d[, cols])
  g <- gt::tab_header(g, title = gt::md(paste0("**", Y, "**", if (!is.null(title)) paste(" \u2014", title))),
                      subtitle = "nested CV  \u00b7  metrics averaged over folds, lower = better  \u00b7  sorted by composite rank")
  g <- gt::cols_label(g, family = "Family", label = "Candidate", Composite = "Composite",
                      mean_rank_LOBO = "LOBO rank", mean_rank_innerCV = "inner-CV rank",
                      LOBO_GD = "GD", LOBO_GD_sd = "GD SD (batch)", LOBO_GD_WinPct = "GD win %",
                      LOBO_MAE = "MAE", LOBO_LL = "LL", LOBO_CLL = "CLL",
                      innerCV_GD = "GD", innerCV_GD_WinPct = "GD win %", innerCV_MAE = "MAE",
                      innerCV_LL = "LL", innerCV_CLL = "CLL", missFitPct = "missfit % (worst fold)")
  g <- gt::tab_spanner(g, "Rank (1 = best)", c("Composite", "mean_rank_LOBO", "mean_rank_innerCV"))
  g <- gt::tab_spanner(g, "LOBO \u2014 held-out batch", gt::starts_with("LOBO_"))
  g <- gt::tab_spanner(g, "Inner CV \u2014 within-distribution", gt::starts_with("innerCV_"))
  g <- gt::fmt_number(g, c("Composite", "mean_rank_LOBO", "mean_rank_innerCV"), decimals = 2)
  g <- gt::fmt_number(g, intersect(c("LOBO_GD", "LOBO_GD_sd", "innerCV_GD"), cols), decimals = 1)
  g <- gt::fmt_number(g, intersect(c("LOBO_MAE", "LOBO_LL", "LOBO_CLL", "innerCV_MAE", "innerCV_LL", "innerCV_CLL"), cols), decimals = 4)
  g <- gt::fmt_number(g, intersect(c("LOBO_GD_WinPct", "innerCV_GD_WinPct", "missFitPct"), cols), decimals = 1)
  g <- gt::sub_values(g, columns = metric_cols, replacement = "Inf",
                      fn = function(x) is.finite(x) & abs(x) > 10 * stats::median(abs(x), na.rm = TRUE))
  if (requireNamespace("scales", quietly = TRUE))
    g <- gt::data_color(g, columns = "Composite", palette = c("#2E7D32", "#F5F5F5"), domain = range(d$Composite, na.rm = TRUE))
  if (!is.null(file)) gt::gtsave(g, file)
  g
}


#' One-call report of a nested-CV model-selection step
#'
#' Writes to `outdir`: `<name>_folds.csv` (every per-fold metric), `<name>_summary.csv` (ranks and
#' aggregates), `<name>_paired.csv` and `<name>_forest.png` (every candidate vs the reference), and one
#' `<name>_<Y>_table.html` per image. Prints the best candidate per image.
#' @param results A result file of vbgamlss.model_selection_NCV() (path or object).
#' @param outdir Output directory (created if needed).
#' @param reference Reference model for the paired comparison (full formula or label); default
#'   the simplest candidate.
#' @param registry Only for result files without a stored job table, see ncv_folds().
#' @param name Prefix of the output files.
#' @param weights Composite weights, see ncv_summary().
#' @return Invisibly, list(folds, summary, paired, forest).
#' @export
ncv_report <- function(results, outdir, reference = NULL, registry = NULL, name = "NCV",
                       weights = c(LOBO = 0.5, innerCV = 0.5)) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  out <- function(x) file.path(outdir, paste0(name, "_", x))
  folds <- ncv_folds(results, registry)
  summ  <- ncv_summary(folds, weights = weights)
  utils::write.csv(folds, out("folds.csv"), row.names = FALSE)
  utils::write.csv(summ,  out("summary.csv"), row.names = FALSE)
  if (is.null(reference)) reference <- complexity_order(unique(folds$label))[1]
  paired <- ncv_paired(folds, reference)
  utils::write.csv(paired, out("paired.csv"), row.names = FALSE)

  forest <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    forest <- plot_ncv_forest(paired, title = paste(name, "\u00b7 reference:", reference))
    n_cand <- length(unique(folds$label)); n_y <- length(unique(folds$Y))
    n_met <- length(unique(forest$data$metric))
    ggplot2::ggsave(out("forest.png"), forest, width = 4 + 5 * n_met, height = 1.5 + 0.45 * n_cand * n_y, dpi = 200, limitsize = FALSE)
  }
  if (requireNamespace("gt", quietly = TRUE)) {
    for (y in unique(summ$Y)) table_ncv_summary(summ, y, title = name, file = out(paste0(y, "_table.html")))
  }
  cat("Best candidate per image (composite rank):\n")
  b <- summ[summ$best, ]
  for (i in seq_len(nrow(b))) cat(sprintf("  %-4s %s  (composite %.2f)\n", b$Y[i], b$label[i], b$Composite[i]))
  cat("Written to", outdir, "\n")
  invisible(list(folds = folds, summary = summ, paired = paired, forest = forest))
}


# simplest first: intercept-only, then by number of terms, by-variable terms after their plain version
complexity_order <- function(labels) {
  n_terms <- lengths(regmatches(labels, gregexpr("\\+", labels)))
  labels[order(labels != "1", n_terms, grepl("by\\s*=", labels), nchar(labels))]
}
