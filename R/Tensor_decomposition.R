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
#' Tucker Decomposition with Background Imputation
#'
#' @param tnsr Tensor with K modes
#' @param mask Optional Tensor object of the same dimensions as tnsr (1 = valid, 0 = background). Default is NULL.
#' @param ranks a vector of the modes of the output core Tensor
#' @param max_iter maximum number of iterations if error stays above \code{tol}, requires ~300 iteration for convergence on brain images
#' @param tol relative Frobenius norm error tolerance
#' @param decay decay factor (damping) for the background imputation [1, 0]. This forces the background imputation to smoothly decay toward zero, anchoring mask edges (eg. 0.95)
#' @export
bi_tucker <- function(tnsr, mask=NULL, ranks=NULL, max_iter=500, tol=1e-5, decay=1){
  stopifnot(is(tnsr,"Tensor"))
  if(is.null(ranks)) stop("ranks must be specified")
  if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks<=0)!=0) stop("ranks must be positive")

  # Use ::: because .is_zero_tensor is an unexported internal function in rTensor
  if (rTensor:::.is_zero_tensor(tnsr)) stop("Zero tensor detected")

  # --- 1. Tensor setup---
  if (!is.null(mask)) {
    stopifnot(is(mask, "Tensor"))
    if (!all(mask@modes == tnsr@modes)) stop("Mask dimensions must strictly match tnsr dimensions")

    # Extract the array from the mask tensor for fast logical indexing
    logical_mask <- as.logical(mask@data)
    tnsr@data[!logical_mask] <- 0
  }

  # Save a robust, untouched copy of the original masked data
  # to compute the true final residuals later.
  orig_data <- tnsr@data
  tnsr_norm <- rTensor::fnorm(tnsr)

  # initialization via truncated hosvd
  num_modes <- tnsr@num_modes
  U_list <- vector("list",num_modes)
  for(m in 1:num_modes){
    temp_mat <- rTensor::rs_unfold(tnsr,m=m)@data
    U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
  }

  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  modes <- tnsr@modes
  modes_seq <- 1:num_modes

  # progress bar
  pb <- txtProgressBar(min=0, max=max_iter, style=3)

  # --- 2. Main ALS Loop ---
  while((curr_iter <= max_iter) && (!converged)){
    setTxtProgressBar(pb,curr_iter)

    # Update the background of our working tensor to the model's prediction.
    if (!is.null(mask) && curr_iter > 1) {
      tnsr@data[!logical_mask] <- est@data[!logical_mask] * decay
    }

    for(m in modes_seq){
      # core Z minus mode m
      X <- rTensor::ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
      # truncated SVD of X
      U_list[[m]] <- svd(rTensor::rs_unfold(X,m=m)@data,nu=ranks[m])$u
    }

    # compute core tensor Z
    Z <- rTensor::ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)

    # calculate the current estimate
    est <- rTensor::ttl(Z, U_list, ms=1:num_modes)

    # iteration residual (used only for convergence tracking)
    if (!is.null(mask)) {
      # Calculate norm strictly inside mask
      curr_resid <- sqrt(sum((orig_data[logical_mask] - est@data[logical_mask])^2))
    } else {
      curr_resid <- rTensor::fnorm(tnsr - est)
    }

    fnorm_resid[curr_iter] <- curr_resid

    if (curr_iter > 1 && (abs(curr_resid - fnorm_resid[curr_iter-1]) / tnsr_norm < tol)) {
      converged <- TRUE
      setTxtProgressBar(pb,max_iter)
    } else {
      curr_iter <- curr_iter + 1
    }
  }
  close(pb)

  # --- Convergence Warning ---
  if (!converged) {
    warning(paste("Failed to converge within", max_iter, "iterations for ranks", paste(ranks, collapse = "x")))
  }

  # --- 3. Final Mask Cleanup ---
  # Zero out the hallucinated background so the returned estimate is clean
  if (!is.null(mask)) {
    est@data[!logical_mask] <- 0
  }

  # --- 4. Variance & Residual Calculation (Uncentered & Mask-Safe) ---
  # Compute exact full-dimensional residuals for the final output array
  true_residual_data <- orig_data - est@data

  if (!is.null(mask)) {
    # Isolate only the valid voxels for strict metric calculations
    valid_orig  <- orig_data[logical_mask]
    valid_est   <- est@data[logical_mask]
    valid_resid <- valid_orig - valid_est

    n_elements <- length(valid_orig)
    tnsr_ms    <- sum(valid_orig^2) / n_elements

    true_sse   <- sum(valid_resid^2)
    true_sst   <- sum(valid_orig^2)

  } else {
    # Full tensor calculations if no mask is provided
    n_elements <- prod(tnsr@modes)
    tnsr_ms    <- sum(as.vector(orig_data)^2) / n_elements

    true_sse   <- sum(true_residual_data^2)
    true_sst   <- sum(orig_data^2)
    cum_energy <- sum(est@data^2) / true_sst
  }

  # --- UNCENTERED METRICS ---
  mse        <- true_sse / n_elements
  rmse       <- sqrt(mse)
  smse       <- mse / tnsr_ms

  # Calculate robust final metrics (Fraction of total energy captured)
  explained_energy <- 1 - (true_sse / true_sst)
  final_fnorm_resid  <- sqrt(true_sse)
  norm_percent       <- (1 - (final_fnorm_resid / sqrt(true_sst))) * 100

  # Filter out empty pre-allocated zeros from the convergence tracker
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]

  invisible(list(
    Z = Z,
    U = U_list,
    conv = converged,
    est = est,
    residuals = true_residual_data,
    norm_percent = norm_percent,
    explained_energy = explained_energy,
    MSE = mse,
    RMSE = rmse,
    SMSE = smse,
    fnorm_resid = final_fnorm_resid,
    all_resids = fnorm_resid
  ))
}






#' Save Tucker Decomposition Results to NIfTI with Spatial Metadata
#'
#' Saves the core tensor, factor matrices, residuals, and reconstructed image to NIfTI files.
#' Uses a reference image to ensure the spatial maps overlay correctly in visualization tools.
#'
#' @param tucker_res The output list from the `tucker_masked` function.
#' @param out_dir A string indicating the directory where files should be saved.
#' @param prefix A string prefix to append to the filenames (e.g., "patient1_"). Default is "".
#' @param reference_image_path Optional path to a NIfTI file to copy coordinate system from.
#' @export
save_tucker_results <- function(tucker_res, out_dir = ".", prefix = "", reference_image_path = NULL) {

  library(ANTsR)
  library(rTensor)

  # Ensure the output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Helper function to cleanly construct file paths: /out_dir/prefix_filename.nii.gz
  make_path <- function(suffix) {
    file.path(out_dir, paste0(prefix, suffix))
  }

  # Load reference image if provided
  ref_img <- NULL
  if (!is.null(reference_image_path) && file.exists(reference_image_path)) {
    message(sprintf("Loading spatial reference from: %s", basename(reference_image_path)))
    ref_img <- antsImageRead(reference_image_path)
  } else if (!is.null(reference_image_path)) {
    warning("Reference image path provided but file not found. Saving without spatial metadata.")
  }

  # --- 1. Save the Core Tensor (Z) ---
  core_array <- tucker_res$Z@data
  core_img   <- as.antsImage(core_array)

  out_core <- make_path("tucker_core.nii.gz")
  antsImageWrite(core_img, out_core)
  cat(" - Saved Core:", out_core, "\n")

  # --- 2. Save the 3 Spatial U Basis Matrices ---
  U_list <- tucker_res$U

  for (i in seq_along(U_list)) {
    out_name <- make_path(sprintf("tucker_U%d_basis.nii.gz", i))
    antsImageWrite(as.antsImage(U_list[[i]]), out_name)
    cat(sprintf(" - Saved Factor U%d: %s\n", i, out_name))
  }

  # --- 3. Save Residuals ---
  if (!is.null(tucker_res$residuals)) {
    resid_array <- tucker_res$residuals

    if (!is.null(ref_img)) {
      # Use reference to copy Origin, Spacing, and Direction
      resid_img <- as.antsImage(resid_array, reference = ref_img)
    } else {
      resid_img <- as.antsImage(resid_array)
    }

    out_resid <- make_path("tucker_residuals.nii.gz")
    antsImageWrite(resid_img, out_resid)
    cat(" - Saved Residuals:", out_resid, "\n")
  }

  # --- 4. Save Reconstructed Image ---
  # Rebuild the full spatial tensor mathematically from the compressed core and bases
  recon_tnsr <- ttl(tucker_res$Z, U_list, ms = seq_along(U_list))
  recon_array <- recon_tnsr@data

  if (!is.null(ref_img)) {
    recon_img <- as.antsImage(recon_array, reference = ref_img)
  } else {
    recon_img <- as.antsImage(recon_array)
  }

  out_recon <- make_path("tucker_reconstructed.nii.gz")
  antsImageWrite(recon_img, out_recon)
  cat(" - Saved Reconstructed Image:", out_recon, "\n")
}



#' Estimate Optimal Tucker Rank using Multi-Level Elbow Detection
#'
#' Performs a two-pass search (coarse then fine) to find the optimal rank scale
#' for a Tucker decomposition of a 3D brain image.
#'
#' @param tnsr input rTensor.
#' @param tnsr_mask input rTensor mask.
#' @param max_iter Numeric. Maximum number of iterations for the bi_tucker algorithm. Default is 500.
#' @param tol Numeric. Relative Frobenius norm error tolerance for convergence. Default is 1e-5.
#'
#' @return A list containing:
#' \itemize{
#'   \item{optimal}{A data frame row with the parameters of the detected elbow rank.}
#'   \item{optimal_RMSE}{Numeric. The masked Root Mean Square Error at the optimal scale.}
#'   \item{optimal_SMSE}{Numeric. The masked Standardized Mean Square Error at the optimal scale.}
#'   \item{full_stats}{A data frame containing metrics for all tested scales.}
#'   \item{tnsr_dim}{A numeric vector of the original tensor dimensions.}
#' }
#' @export
estimate_optimal_tucker_rank <- function(tnsr,
                                         tnsr_mask,
                                         max_iter = 500,
                                         tol = 1e-5) {

  # ============================================================================
  # HELPER: FIND ELBOW -> should i normalize axes??
  # ============================================================================
  find_elbow_rank <- function(results_df, x_col = "total_params", y_col = "SMSE") {
    results_df <- results_df[order(results_df[[x_col]]), ]
    x <- results_df[[x_col]]; y <- results_df[[y_col]]
    x1 <- x[1]; y1 <- y[1]; xn <- x[length(x)]; yn <- y[length(y)]

    numerator <- abs((xn - x1) * (y1 - y) - (x1 - x) * (yn - y1))
    denominator <- sqrt((xn - x1)^2 + (yn - y1)^2)

    results_df$elbow_distance <- numerator / denominator
    optimal_index <- which.max(results_df$elbow_distance)

    return(results_df[optimal_index, ])
  }

  # ============================================================================
  # PART 1: LOAD AND PREPARE MASKED DATA
  # ============================================================================
  stopifnot(is(tnsr,"Tensor"))
  stopifnot(is(tnsr_mask,"Tensor"))

  tnsr_dim <- dim(tnsr)

  # Create a 1D logical vector from the mask
  mask_vec <- as.vector(tnsr_mask@data) > 0
  n_mask   <- sum(mask_vec)

  # Calculate baseline variance ONLY within the mask for accurate SMSE
  tnsr_var <- var(as.vector(tnsr@data)[mask_vec])

  # ============================================================================
  # HELPER: EVALUATE A SPECIFIC GRID OF SCALES
  # ============================================================================
  evaluate_grid <- function(scales_to_test) {
    results_list <- list()
    for (i in seq_along(scales_to_test)) {
      scale <- scales_to_test[i]
      current_ranks <- pmax(1, floor(tnsr_dim * scale))
      message(sprintf("  Scale %.2f | Ranks: %s", scale, paste(current_ranks, collapse="x")))

      tucker_res <- bi_tucker(tnsr = tnsr, mask = tnsr_mask, ranks = current_ranks,
                                  max_iter = max_iter, tol = tol)

      expl_energy <- tucker_res$explained_energy * 100

      # Safely extract residuals and compute RSS strictly within the mask
      res_data <- if(is(tucker_res$residuals, "Tensor")) tucker_res$residuals@data else tucker_res$residuals
      res_vec  <- as.vector(res_data)
      rss      <- sum(res_vec[mask_vec]^2)

      k_core    <- prod(current_ranks)
      k_factors <- sum(tnsr_dim * current_ranks)
      k         <- k_core + k_factors

      results_list[[i]] <- data.frame(
        scale = scale,
        converged = tucker_res$conv,
        rank_core = k_core,
        total_params = k,
        rank_str = paste(current_ranks, collapse="x"),
        explained_energy = expl_energy,
        MSE = tucker_res$MSE,
        RMSE = tucker_res$RMSE,
        SMSE = tucker_res$SMSE,
        AIC = 2 * k + n_mask * log(rss/n_mask),
        BIC = k * log(n_mask) + n_mask * log(rss/n_mask)
      )
    }
    return(do.call(rbind, results_list))
  }

  # ============================================================================
  # PART 2: COARSE SEARCH
  # ============================================================================
  message("--- Coarse Search---")
  coarse_scales <- seq(0.05, 0.95, by = 0.1)
  coarse_stats  <- evaluate_grid(coarse_scales)
  coarse_optimal <- find_elbow_rank(coarse_stats)

  message(sprintf(">> Coarse Elbow Scale: %.2f", coarse_optimal$scale))

  # ============================================================================
  # PART 3: FINE SEARCH
  # ============================================================================
  message("--- Fine Search ---")
  fine_scales_raw <- seq(max(0.05, coarse_optimal$scale - 0.09),
                         min(0.95, coarse_optimal$scale + 0.09),
                         by = 0.01)

  fine_scales <- setdiff(round(fine_scales_raw, 3), round(coarse_scales, 3))

  if (length(fine_scales) > 0) {
    fine_stats  <- evaluate_grid(fine_scales)
    final_stats <- rbind(coarse_stats, fine_stats)
  } else {
    final_stats <- coarse_stats
  }

  final_stats <- final_stats[order(final_stats$scale), ]
  true_optimal <- find_elbow_rank(final_stats)

  message(sprintf("\n>> FINAL Optimal Scale: %.2f | Explained Variance: %.2f%% | Masked RMSE: %.4f | Masked SMSE: %.4f",
                  true_optimal$scale, true_optimal$explained_energy, true_optimal$RMSE, true_optimal$SMSE))

  # Return specific optimized values as separate list items for easy downstream access
  return(list(
    optimal = true_optimal,
    optimal_RMSE = true_optimal$RMSE,
    optimal_SMSE = true_optimal$SMSE,
    full_stats = final_stats,
    tnsr_dim = tnsr_dim
  ))
}




#' Plot Tucker Decomposition Elbow Results
#'
#' Generates a dual-axis ggplot2 visualization from the output of
#' estimate_optimal_tucker_rank().
#'
#' @param elbow_results List. The exact object returned by estimate_optimal_tucker_rank().
#' @return A ggplot2 object.
#' @export
plot_tucker_elbow <- function(elbow_results) {

  # Extract data from the results list
  final_stats  <- elbow_results$full_stats
  true_optimal <- elbow_results$optimal

  # 1. Format a clean information label (Total Params removed)
  optimal_info <- sprintf(
    "Elbow Selected:\nScale: %.2f\nRanks: %s\nCore Params: %s\nExpl. Variance: %.2f%%",
    true_optimal$scale,
    true_optimal$rank_str,
    format(true_optimal$rank_core, big.mark = ","),
    true_optimal$explained_energy
  )

  # 2. Robust helper function for the top axis
  # This looks up the exact parameter count from your dataframe
  # rather than trying to recalculate it.
  get_core_params <- function(breaks) {
    sapply(breaks, function(b) {
      # Find the row in final_stats that matches the break scale
      closest_idx <- which.min(abs(final_stats$scale - b))
      format(final_stats$rank_core[closest_idx], big.mark = ",", scientific = FALSE)
    })
  }

  # Where the ticks should appear on the top axis
  top_axis_breaks <- seq(0.1, 0.9, by = 0.2)

  # 3. Generate the plot
  p1 <- ggplot(final_stats, aes(x = scale, y = explained_energy)) +
    geom_line(color = "#2c3e50", linewidth = 1) +
    geom_point(color = "#2c3e50", size = 2) +

    # Highlight optimal point & elbow vertical line
    geom_vline(xintercept = true_optimal$scale, linetype = "dotted", color = "#e74c3c", linewidth = 1) +
    geom_point(data = true_optimal, aes(x = scale, y = explained_energy), color = "#e74c3c", size = 4) +

    # Bottom-Right Info Box
    annotate("label", x = Inf, y = -Inf, label = optimal_info,
             hjust = 1.05, vjust = -0.1, fill = "white", color = "#2c3e50", fontface = "bold", size = 4) +

    # Threshold lines
    geom_hline(yintercept = c(95, 99, 99.5), linetype = "dashed", color = c("#e74c3c", "#f39c12", "#27ae60"), alpha = 0.6) +
    annotate("text", x = min(final_stats$scale), y = 95.2, label = "95%", color = "#e74c3c", hjust = 0, vjust = 0, fontface = "bold") +
    annotate("text", x = min(final_stats$scale), y = 99.2, label = "99%", color = "#f39c12", hjust = 0, vjust = 0, fontface = "bold") +
    annotate("text", x = min(final_stats$scale), y = 99.7, label = "99.5%", color = "#27ae60", hjust = 0, vjust = 0, fontface = "bold") +

    # Dual X-Axis configuration
    scale_x_continuous(
      breaks = seq(0, 1, 0.1),
      name = "Rank Scale",
      sec.axis = sec_axis(
        trans = ~ .,
        breaks = top_axis_breaks,
        labels = function(x) get_core_params(x),
        name = "Number of Core Parameters"
      )
    ) +
    scale_y_continuous(name = "Explained Variance (%)") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title.x.top = element_text(color = "#2980b9", face = "bold", margin = margin(b = 10)),
      axis.text.x.top = element_text(color = "#2980b9")
    ) +
    labs(title = "Tucker Decomposition: Explained Variance vs. Scale & Complexity",
         subtitle = "Top axis displays total core parameter count. Red dot denotes the mathematical elbow.")

  return(p1)
}



# ============================================================================ #
# EXPERIMENTAL
# ============================================================================ #


#' Project a 4D Image onto a 3D Tucker Basis with Advanced Diagnostics
#'
#' Projects a 4D tensor onto a 3D spatial Tucker basis using uncentered data.
#' Calculates global metrics. Optionally saves NIfTI maps for Uncentered R2, MaxAE,
#' mean/median errors, standard deviation/MAD, and an average reshuffled AR1 map.
#'
#' @param img_tnsr An `rTensor` object representing the 4D image.
#' @param basis_list A list of 3 spatial factor matrices (U1, U2, U3).
#' @param mask_tnsr Optional `rTensor` 3D mask (1 = valid, 0 = background).
#' @param ref_image_path Optional path to a reference NIfTI file to copy coordinates.
#' @param compute_diagnostics Logical. If TRUE (default), computes and saves 3D summary maps and 1D metrics.
#' @param n_reshuffles Integer. Number of 4th-dimension reshuffles to average for the baseline AR1 map.
#' @param out_dir String indicating the output directory. Default is ".".
#' @param prefix String prefix for saved filenames. Default is "".
#' @export
project_to_tucker_basis <- function(img_tnsr, basis_list, mask_tnsr = NULL,
                                    ref_image_path = NULL, compute_diagnostics = TRUE,
                                    n_reshuffles = 25, out_dir = ".", prefix = "") {

  library(rTensor)
  library(ANTsR)

  # --- 1. Input Validation & Setup ---
  stopifnot(is(img_tnsr, "Tensor"))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  make_path <- function(suffix) file.path(out_dir, paste0(prefix, suffix))

  ref_img <- NULL
  if (!is.null(ref_image_path) && file.exists(ref_image_path)) {
    ref_img <- antsImageRead(ref_image_path)
  }

  # --- 2. Masking ---
  logical_mask <- NULL
  if (!is.null(mask_tnsr)) {
    mask_arr <- mask_tnsr@data != 0
    if (img_tnsr@num_modes == 4 && length(dim(mask_arr)) == 3) {
      logical_mask <- replicate(img_tnsr@modes[4], mask_arr)
    } else {
      logical_mask <- mask_arr
    }
    img_tnsr@data[!logical_mask] <- 0
  }

  # --- 3. Forward Projection & Reconstruction ---
  target_modes <- seq_along(basis_list)
  transposed_basis <- lapply(basis_list, t)

  core_tnsr <- rTensor::ttl(img_tnsr, transposed_basis, ms = target_modes)
  recon_tnsr <- rTensor::ttl(core_tnsr, basis_list, ms = target_modes)
  resid_array <- img_tnsr@data - recon_tnsr@data

  if (!is.null(logical_mask)) {
    recon_tnsr@data[!logical_mask] <- 0
    resid_array[!logical_mask] <- 0
  }

  # --- 4. Global Metrics: 4D Uncentered Explained Energy ---
  if (!is.null(logical_mask)) {
    valid_orig  <- img_tnsr@data[logical_mask]
    valid_resid <- resid_array[logical_mask]
  } else {
    valid_orig  <- img_tnsr@data
    valid_resid <- resid_array
  }

  true_sse <- sum(valid_resid^2)
  true_sst <- sum(valid_orig^2)
  explained_energy <- 1 - (true_sse / true_sst)
  message(sprintf(">> 4D Cumulative Explained Energy (Uncentered): %.2f%%", explained_energy * 100))

  # --- Helper: Safe Spatial Mapping ---
  create_safe_ants_img <- function(array_data, reference_image) {
    if (is.null(reference_image)) return(as.antsImage(array_data))
    dim_data <- length(dim(array_data))
    dim_ref <- reference_image@dimension

    if (dim_data == dim_ref) {
      return(as.antsImage(array_data, reference = reference_image))
    }
    if (dim_data == 4 && dim_ref == 3) {
      spc <- c(antsGetSpacing(reference_image), 1)
      orig <- c(antsGetOrigin(reference_image), 0)
      dir_3d <- antsGetDirection(reference_image)
      dir_4d <- diag(4)
      dir_4d[1:3, 1:3] <- dir_3d
      return(as.antsImage(array_data, spacing = spc, origin = orig, direction = dir_4d))
    }
    return(as.antsImage(array_data))
  }

  # --- 5. Save Core & 4D Data ---
  antsImageWrite(as.antsImage(core_tnsr@data), make_path("projected_core.nii.gz"))
  antsImageWrite(create_safe_ants_img(resid_array, ref_img), make_path("projected_residuals_4D.nii.gz"))
  antsImageWrite(create_safe_ants_img(recon_tnsr@data, ref_img), make_path("projected_reconstruction_4D.nii.gz"))

  # --- 6. Advanced Diagnostics (Optional) ---
  if (compute_diagnostics && length(dim(resid_array)) == 4) {
    message("--- Computing 3D Summary Maps and 1D Subject Metrics ---")
    T_dim <- dim(resid_array)[4]
    vol_dims <- dim(resid_array)[1:3]

    # A. Standard Spatial Metrics
    mean_resid <- apply(abs(resid_array), c(1, 2, 3), mean, na.rm = TRUE)
    median_ae  <- apply(abs(resid_array), c(1, 2, 3), median, na.rm = TRUE)

    sd_resid   <- apply(resid_array, c(1, 2, 3), sd, na.rm = TRUE)
    mad_resid  <- apply(resid_array, c(1, 2, 3), mad, na.rm = TRUE) # Median Absolute Deviation (robust std)
    median_raw <- apply(resid_array, c(1, 2, 3), median, na.rm = TRUE)

    max_ae     <- apply(abs(resid_array), c(1, 2, 3), max, na.rm = TRUE)

    sd_resid[is.na(sd_resid)] <- 0
    mad_resid[is.na(mad_resid)] <- 0

    # UNCENTERED Voxel-wise R2 Map
    sum_resid2 <- apply(resid_array^2, c(1, 2, 3), sum, na.rm = TRUE)
    sum_y2     <- apply(img_tnsr@data^2, c(1, 2, 3), sum, na.rm = TRUE)
    r2_map     <- 1 - (sum_resid2 / sum_y2)
    r2_map[is.na(r2_map) | sum_y2 == 0] <- 0

    # B. Fast AR1 Computation & Reshuffling Estimation
    message(sprintf("--- Computing Mean AR1 of %d Reshuffles ---", n_reshuffles))

    mask_idx <- if (!is.null(mask_tnsr)) which(mask_tnsr@data != 0) else 1:prod(vol_dims)
    V_dim <- length(mask_idx)
    resid_mat <- matrix(0, nrow = V_dim, ncol = T_dim)

    for (t in seq_len(T_dim)) {
      resid_mat[, t] <- resid_array[, , , t][mask_idx]
    }

    fast_ar1 <- function(mat) {
      N <- ncol(mat) - 1
      mat_t1 <- mat[, 1:N, drop = FALSE]
      mat_t2 <- mat[, 2:(N+1), drop = FALSE]

      sum_t1 <- rowSums(mat_t1)
      sum_t2 <- rowSums(mat_t2)
      sum_t1_sq <- rowSums(mat_t1^2)
      sum_t2_sq <- rowSums(mat_t2^2)
      sum_t1_t2 <- rowSums(mat_t1 * mat_t2)

      num <- N * sum_t1_t2 - (sum_t1 * sum_t2)
      den <- sqrt(pmax(N * sum_t1_sq - sum_t1^2, 0) * pmax(N * sum_t2_sq - sum_t2^2, 0))
      ar1 <- num / den
      ar1[is.na(ar1) | den == 0] <- 0
      return(ar1)
    }

    # Average of Reshuffled AR1 Maps
    reshuffled_ar1_sum <- numeric(V_dim)
    for (i in seq_len(n_reshuffles)) {
      shuffled_mat <- resid_mat[, sample(T_dim)]
      reshuffled_ar1_sum <- reshuffled_ar1_sum + fast_ar1(shuffled_mat)
    }
    reshuffled_ar1_mean_vec <- reshuffled_ar1_sum / n_reshuffles

    # Reconstruct 3D NIfTI map for reshuffled AR1
    ar1_map_reshuffled <- array(0, dim = vol_dims)
    ar1_map_reshuffled[mask_idx] <- reshuffled_ar1_mean_vec

    # Save 3D Maps
    antsImageWrite(create_safe_ants_img(mean_resid, ref_img), make_path("projected_residuals_mean_3D.nii.gz"))
    antsImageWrite(create_safe_ants_img(median_ae, ref_img),  make_path("projected_residuals_medianAE_3D.nii.gz"))
    antsImageWrite(create_safe_ants_img(median_raw, ref_img), make_path("projected_residuals_medianRaw_3D.nii.gz"))

    antsImageWrite(create_safe_ants_img(sd_resid, ref_img),   make_path("projected_residuals_sd_3D.nii.gz"))
    antsImageWrite(create_safe_ants_img(mad_resid, ref_img),  make_path("projected_residuals_mad_3D.nii.gz"))

    antsImageWrite(create_safe_ants_img(max_ae, ref_img),     make_path("projected_residuals_maxAE_3D.nii.gz"))
    antsImageWrite(create_safe_ants_img(r2_map, ref_img),     make_path("projected_R2_map_uncentered_3D.nii.gz"))

    antsImageWrite(create_safe_ants_img(ar1_map_reshuffled, ref_img), make_path("projected_AR1_reshuffled_mean_3D.nii.gz"))

    # C. Subject-wise Metrics (1D CSV)
    rmse_per_subj       <- numeric(T_dim)
    median_se_per_subj  <- numeric(T_dim)
    mean_orig_per_subj  <- numeric(T_dim)
    mean_recon_per_subj <- numeric(T_dim)

    for (t in seq_len(T_dim)) {
      v_orig  <- img_tnsr@data[, , , t]
      v_recon <- recon_tnsr@data[, , , t]
      v_resid <- resid_array[, , , t]

      if (!is.null(logical_mask)) {
        m_t <- logical_mask[, , , t]
        v_orig  <- v_orig[m_t]
        v_recon <- v_recon[m_t]
        v_resid <- v_resid[m_t]
      }

      rmse_per_subj[t]       <- sqrt(mean(v_resid^2, na.rm = TRUE))
      median_se_per_subj[t]  <- median(v_resid, na.rm = TRUE)
      mean_orig_per_subj[t]  <- mean(v_orig, na.rm = TRUE)
      mean_recon_per_subj[t] <- mean(v_recon, na.rm = TRUE)
    }

    diff_signal <- mean_recon_per_subj - mean_orig_per_subj
    percent_signal_change <- (diff_signal / mean_orig_per_subj) * 100

    metrics_df <- data.frame(
      Subject_Index = seq_len(T_dim),
      RMSE = rmse_per_subj,
      Median = median_se_per_subj,
      Mean_Original_Signal = mean_orig_per_subj,
      Mean_Reconstructed_Signal = mean_recon_per_subj,
      Signal_Difference = diff_signal,
      Percent_Signal_Change = percent_signal_change
    )

    csv_out <- make_path("projected_subject_metrics.csv")
    write.csv(metrics_df, csv_out, row.names = FALSE)
    cat(" - Saved 1D subject metrics to:", csv_out, "\n")
  }

  invisible(list(
    core = core_tnsr,
    reconstruction = recon_tnsr,
    residuals = resid_array,
    explained_energy = explained_energy
  ))
}


#' Project Tucker Basis to Reconstruct a Tensor
#'
#' Reconstructs a full 3D or 4D tensor from a Tucker decomposition core
#' and factor matrices. Optionally applies a spatial mask and saves the
#' result as a NIfTI file with matching headers.
#'
#' @param tucker_res List. The output from a Tucker decomposition. Must contain
#'                   the core tensor (named 'Z' or 'core') and factor matrices
#'                   (named 'U' or 'factors').
#' @param ref_image_path Character (Optional). Path to a reference NIfTI image
#'                       to copy spatial metadata (spacing, origin, direction).
#' @param mask_path Character (Optional). Path to a NIfTI mask to apply to the
#'                  reconstructed tensor (sets background to 0).
#' @param out_path Character (Optional). Path to save the output NIfTI file
#'                 (e.g., "/path/to/reconstructed.nii.gz").
#' @return An `rTensor` object containing the reconstructed tensor.
#' @export
reconstruct_tucker_tensor <- function(tucker_res,
                                      ref_image_path = NULL,
                                      mask_path = NULL,
                                      out_path = NULL) {

  library(rTensor)
  library(ANTsR)

  # 1. Extract Core and Factors (Handles both rTensor and custom naming conventions)
  core_tnsr <- if (!is.null(tucker_res$Z)) tucker_res$Z else tucker_res$core
  factors   <- if (!is.null(tucker_res$U)) tucker_res$U else tucker_res$factors

  if (is.null(core_tnsr) || is.null(factors)) {
    stop("Error: tucker_res must contain the core tensor ('Z' or 'core') and factor matrices ('U' or 'factors').")
  }

  # 2. Reconstruct the Tensor via Tensor-Times-List (ttl)
  message("--- Step 1: Reconstructing tensor from Tucker basis ---")
  num_modes <- length(factors)

  # Mathematically: G x_1 A^(1) x_2 A^(2) x_3 A^(3) ...
  recon_tnsr <- ttl(core_tnsr, factors, ms = 1:num_modes)
  recon_array <- recon_tnsr@data

  # 3. Apply Spatial Mask
  if (!is.null(mask_path)) {
    message(sprintf("--- Step 2: Applying mask (%s) ---", basename(mask_path)))
    mask_img <- antsImageRead(mask_path)
    mask_arr <- as.array(mask_img)

    # Broadcast 3D mask to 4D if reconstructing a 4D timeseries/multimodal tensor
    if (length(dim(recon_array)) == 4 && length(dim(mask_arr)) == 3) {
      message("  -> 4D Tensor detected: Broadcasting 3D mask across the 4th dimension.")
      mask_arr <- replicate(dim(recon_array)[4], mask_arr)
    }

    # Verify dimensions before multiplication
    if (!all(dim(recon_array) == dim(mask_arr))) {
      warning("Mask dimensions do not match reconstructed tensor dimensions. Skipping mask application.")
    } else {
      recon_array <- recon_array * mask_arr
    }
  }

  # 4. Save to NIfTI with ANTsR
  if (!is.null(out_path)) {
    message("--- Step 3: Saving to NIfTI ---")

    if (!is.null(ref_image_path)) {
      # Steal the spatial header from the original image
      ref_img <- antsImageRead(ref_image_path)
      recon_img <- as.antsImage(recon_array, reference = ref_img)
    } else {
      warning("No reference image provided. Saving with default ANTs spatial header (1x1x1 spacing).")
      recon_img <- as.antsImage(recon_array)
    }

    antsImageWrite(recon_img, out_path)
    message(sprintf(">> Saved successfully to: %s", out_path))
  } else {
    message(">> Reconstruction complete. (NIfTI save skipped: no out_path provided).")
  }

  # Update the rTensor object with the masked data before returning
  recon_tnsr@data <- recon_array
  return(recon_tnsr)
}












