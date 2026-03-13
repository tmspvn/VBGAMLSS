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
#
#
# == Low-dimensional parameter manifold / low-rank initialization ==
# Obtain good voxel-wise initialization by learning a low-dimensional
# relationship between voxel signal and fitted parameters using a small
# set of fully fitted voxels.
#' @export
manifold_initialization <- function(imageframe,
                                    train_data,
                                    g_formula,
                                    g_family = SHASH,
                                    n_components = 5,
                                    n_clusters = 6,
                                    n_anchors = 500,       # Defaults to sqrt(voxels) if left NULL
                                    save_images = FALSE,   # Set to TRUE to save the NIfTI files
                                    donor_mask = NULL,     # path to donor mask nifti for saving coefs maps
                                    output_dir = getwd(),
                                    num_threads = NULL,
                                    ...) {

  # ---------------------------------------------------------
  # 1. Data Loading
  # ---------------------------------------------------------
  # Ensure imageframe is a matrix
  if (!is.matrix(imageframe)) {
    imageframe <- as.matrix(imageframe)
  }

  # ---------------------------------------------------------
  # 2. PCA & Stratified Sampling (Optimized)
  # ---------------------------------------------------------
  message(sprintf("Running PCA (%d components) and K-Means (%d clusters)", n_components, n_clusters))

  # Transpose so rows = voxels, cols = signal features
  signal_matrix_t <- t(imageframe)

  # Run PCA ONCE
  pca_result <- prcomp(signal_matrix_t, center = TRUE, scale. = FALSE, rank. = n_components)
  latent_signals <- pca_result$x

  # Fast Clustering
  set.seed(42)
  kmeans_result <- kmeans(latent_signals, centers = n_clusters, iter.max = 2000, nstart = 20)
  cluster_labels <- kmeans_result$cluster

  # Stratified Selection
  anchors_per_cluster <- floor(n_anchors / n_clusters)
  anchor_indices <- integer(0)

  for (c in 1:n_clusters) {
    indices_in_cluster <- which(cluster_labels == c)
    n_to_pick <- min(length(indices_in_cluster), anchors_per_cluster)

    if (length(indices_in_cluster) > 1) {
      picked <- sample(indices_in_cluster, size = n_to_pick, replace = FALSE)
    } else {
      picked <- indices_in_cluster
    }
    anchor_indices <- c(anchor_indices, picked)
  }

  message(sprintf("Selected %d anchor voxels.", length(anchor_indices)))

  # ---------------------------------------------------------
  # 3. Optional: Save Anchor Mask
  # ---------------------------------------------------------
  if (save_images) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    anchor_out <- file.path(output_dir, "mask_anchor_voxels.nii.gz")

    donor_mask <- antsImageRead(donor_mask)
    anchor_image <- donor_mask * 0
    is_tissue <- as.array(donor_mask) > 0
    tissue_vals <- rep(0, sum(is_tissue))
    tissue_vals[anchor_indices] <- 1
    anchor_image[is_tissue] <- tissue_vals

    antsImageWrite(anchor_image, anchor_out)
    message("Saved anchor mask to: ", anchor_out)
  }

  # ---------------------------------------------------------
  # 4. Fit Anchor Voxels
  # ---------------------------------------------------------
  message("Running slow nonlinear fits on anchor voxels...")
  imageframe_anchors <- imageframe[, anchor_indices, drop = FALSE]

  D <- vbgamlss(imageframe = imageframe_anchors,
                g.formula = g_formula,
                g.family = g_family,
                num_cores = num_threads,
                train.data = train_data,
                ...
  )

  # Extract Coefficients Dynamically
  unlisted_coeffs <- unlist(D[[1]]$coefficients)
  P_coefs <- length(unlisted_coeffs)
  coef_names <- names(unlisted_coeffs)

  anchor_params <- matrix(NA_real_,
                          nrow = length(D),
                          ncol = P_coefs,
                          dimnames = list(NULL, coef_names))

  for(i in seq_along(D)){
    anchor_params[i, ] <- unlist(D[[i]]$coefficients)
  }

  # ---------------------------------------------------------
  # 5. Train Random Forests & Predict
  # ---------------------------------------------------------
  message("Training Random Forests to map the manifold...")

  train_x <- as.data.frame(latent_signals[anchor_indices, , drop = FALSE])
  full_x <- as.data.frame(latent_signals)
  train_y <- anchor_params

  N_voxels <- nrow(latent_signals)
  predicted_parameters <- matrix(0, nrow = N_voxels, ncol = P_coefs)
  colnames(predicted_parameters) <- coef_names

  for (i in 1:P_coefs) {
    df_train <- cbind(target = train_y[, i], train_x)

    rf_model <- ranger(
      formula = target ~ .,
      data = df_train,
      num.trees = 200,
      min.node.size = 5,
      num.threads = num_threads
    )

    predicted_parameters[, i] <- predict(rf_model, data = full_x)$predictions
    message(sprintf("Model %d/%d ('%s') trained and predicted.", i, P_coefs, coef_names[i]))
  }

  # ---------------------------------------------------------
  # 6. Optional: Save Predicted Maps
  # ---------------------------------------------------------
  if (save_images) {
    message("Saving predicted parameter maps...")

    for (i in 1:P_coefs) {
      clean_name <- gsub("[^A-Za-z0-9]", "_", coef_names[i])
      clean_name <- gsub("_+", "_", clean_name)
      clean_name <- gsub("_$", "", clean_name)

      out_file <- file.path(output_dir, paste0("predicted_", clean_name, ".nii.gz"))

      param_image <- donor_mask * 0
      param_image[as.array(donor_mask) > 0] <- predicted_parameters[, i]
      antsImageWrite(param_image, out_file)
    }
    message("All predicted maps saved.")
  }

  return(predicted_parameters)
}





#' @export
# == average initialization ==
global_initialization <- function(imageframe,
                                  train_data,
                                  g_formula,
                                  g_family = NO,
                                  force_ypositivity = TRUE,
                                  eps = 1e-5,
                                  ...) {
  # Calculate subject-wise mean and attach to covariates
  train_data$Y <- rowMeans(as.matrix(imageframe), na.rm = TRUE)

  if (force_ypositivity) {
    train_data$Y[train_data$Y <= 0] <- eps * 0.001
  }

  # Fit the global model
  glob_fit <- gamlss2::gamlss2(formula = g_formula,
                               data = train_data,
                               family = g_family,
                               light = FALSE,
                               maxit = c(300, 100),
                               control = gamlss2::gamlss2_control(trace = TRUE, eps = eps))

  unlisted_coeffs <- unlist(glob_fit$coefficients)
  P_coefs <- length(unlisted_coeffs)
  coef_names <- names(unlisted_coeffs)

  predicted_parameters <- matrix(
    unlisted_coeffs,
    nrow = dim(imageframe)[2],
    ncol = P_coefs,
    byrow = TRUE
  )

  colnames(predicted_parameters) <- coef_names

  return(predicted_parameters)
}




# == Correlation-length-based multiresolution initialization ==
# Estimate a spatial downsampling factor for voxel fitting based on the
# spatial correlation length of the parameter field.

multiresolution_initialization <- function(imageframe,
                   covs,
                   g_formula,
                   donor_mask,
                   sample_fraction = 0.01, # Default to 1% of voxels
                   max_sample = 1000,      # Memory safeguard for distance matrices
                   save_images = FALSE,
                   output_dir = getwd(),
                   num_threads = 10,
                   ...) {

  # ---------------------------------------------------------
  # 1. Setup & Sample a small subset of voxels
  # ---------------------------------------------------------
  if (!is.matrix(imageframe)) {
    imageframe <- as.matrix(imageframe)
  }

  if (is.character(donor_mask)) {
    donor_mask <- antsImageRead(donor_mask)
  }

  is_tissue <- as.array(donor_mask) > 0
  tissue_idx <- which(is_tissue)
  N_voxels <- length(tissue_idx)

  n_sample <- min(max_sample, max(100, floor(N_voxels * sample_fraction)))
  message(sprintf("Step 1: Sampling %d voxels (~%.2f%%) for initial fit...",
                  n_sample, (n_sample/N_voxels)*100))

  set.seed(42)
  sample_indices <- sample(seq_len(N_voxels), n_sample, replace = FALSE)
  imageframe_sample <- imageframe[, sample_indices, drop = FALSE]

  # Run full nonlinear voxel fit on the sample
  D_sample <- vbgamlss(imageframe = imageframe_sample,
                       g.formula = g_formula,
                       g.family = SHASH, # Assuming identical family to manifold
                       num_cores = num_threads,
                       train.data = covs,
                       ...)

  unlisted_coeffs <- unlist(D_sample[[1]]$coefficients)
  P_coefs <- length(unlisted_coeffs)
  coef_names <- names(unlisted_coeffs)

  sample_params <- matrix(NA_real_, nrow = length(D_sample), ncol = P_coefs)
  for(i in seq_along(D_sample)) {
    sample_params[i, ] <- unlist(D_sample[[i]]$coefficients)
  }

  # ---------------------------------------------------------
  # 2 & 3. Estimate spatial correlation length & variogram range
  # ---------------------------------------------------------
  message("Steps 2 & 3: Estimating spatial correlation length (Variogram)...")

  # Extract physical coordinates
  coords_all <- which(is_tissue, arr.ind = TRUE)
  spacing <- antsGetSpacing(donor_mask)
  coords_phys <- scale(coords_all, center = FALSE, scale = 1/spacing)
  coords_sample <- coords_phys[sample_indices, ]

  # Pairwise spatial distances (h) and parameter differences (d)
  dist_mat <- as.matrix(dist(coords_sample))
  h_vals <- dist_mat[lower.tri(dist_mat)]

  # Standardize parameters so features contribute equally to the difference
  scale_params <- scale(sample_params)
  param_dist_mat <- as.matrix(dist(scale_params))^2
  d_vals <- param_dist_mat[lower.tri(param_dist_mat)]

  # Empirical Semivariogram (30 bins)
  n_bins <- 30
  breaks <- seq(0, max(h_vals), length.out = n_bins + 1)
  bin_centers <- breaks[-1] - diff(breaks)/2
  gamma_h <- numeric(n_bins)

  for(b in 1:n_bins){
    in_bin <- h_vals >= breaks[b] & h_vals < breaks[b+1]
    if(any(in_bin)) gamma_h[b] <- 0.5 * mean(d_vals[in_bin])
  }

  valid_bins <- !is.na(gamma_h) & gamma_h > 0
  bin_centers <- bin_centers[valid_bins]
  gamma_h <- gamma_h[valid_bins]

  # Determine plateau L (95% of the 90th percentile of variance)
  target_variance <- 0.95 * quantile(gamma_h, 0.90, na.rm = TRUE)
  plateau_idx <- which(gamma_h >= target_variance)[1]
  if(is.na(plateau_idx)) plateau_idx <- length(gamma_h)

  L <- bin_centers[plateau_idx]
  message(sprintf("Estimated spatial correlation length L = %.2f mm", L))

  # ---------------------------------------------------------
  # 4 & 5. Choose coarse grid & compute downsampling factor
  # ---------------------------------------------------------
  coarse_spacing <- pmax(L / 2, spacing) # Ensures we don't upsample
  downsample_factor <- coarse_spacing / spacing

  # Convert to integer strides for grid sampling
  strides <- pmax(1, round(downsample_factor))
  message(sprintf("Steps 4 & 5: Downsampling strides set to X:%d, Y:%d, Z:%d",
                  strides[1], strides[2], strides[3]))

  # ---------------------------------------------------------
  # 6. Build coarse grid & Fit
  # ---------------------------------------------------------
  message("Step 6: Fitting model on the coarse grid...")

  # Identify coarse grid voxels (1-based index mapping)
  keep_x <- (coords_all[, 1] - 1) %% strides[1] == 0
  keep_y <- (coords_all[, 2] - 1) %% strides[2] == 0
  keep_z <- (coords_all[, 3] - 1) %% strides[3] == 0

  coarse_indices <- which(keep_x & keep_y & keep_z)
  message(sprintf("Coarse grid contains %d voxels.", length(coarse_indices)))

  imageframe_coarse <- imageframe[, coarse_indices, drop = FALSE]

  D_coarse <- vbgamlss(imageframe = imageframe_coarse,
                       g.formula = g_formula,
                       g.family = SHASH,
                       num_cores = num_threads,
                       train.data = covs,
                       ...)

  coarse_params <- matrix(NA_real_, nrow = length(D_coarse), ncol = P_coefs)
  for(i in seq_along(D_coarse)) {
    coarse_params[i, ] <- unlist(D_coarse[[i]]$coefficients)
  }

  # ---------------------------------------------------------
  # 7 & 8. Interpolate parameters to full resolution
  # ---------------------------------------------------------
  message("Step 7 & 8: Trilinear interpolation to full resolution...")

  predicted_parameters <- matrix(0, nrow = N_voxels, ncol = P_coefs)
  colnames(predicted_parameters) <- coef_names

  dim_fine <- dim(donor_mask)
  dim_coarse <- ceiling(dim_fine / strides)

  # Calculate mapped array coordinates for the coarse fits
  coarse_array_coords <- 1 + (coords_all[coarse_indices, ] - 1) %/%
    matrix(rep(strides, each = length(coarse_indices)), ncol = 3)

  for (i in 1:P_coefs) {
    # Initialize an empty array for the specific parameter
    coarse_array <- array(0, dim = dim_coarse)
    coarse_array[coarse_array_coords] <- coarse_params[, i]

    # Convert to a temporary ANTsImage (scaled spatial spacing)
    coarse_ants <- as.antsImage(coarse_array, spacing = spacing * strides)

    # Interpolate back to original full resolution space
    # interpType = 0 corresponds to linear (trilinear in 3D)
    interp_ants <- resampleImageToTarget(coarse_ants, donor_mask, interpType = 0)

    # Extract back into matrix mapped to tissue indices
    predicted_parameters[, i] <- as.array(interp_ants)[tissue_idx]
  }

  # ---------------------------------------------------------
  # Optional: Save Predicted Maps
  # ---------------------------------------------------------
  if (save_images) {
    message("Saving CLBMRI predicted parameter maps...")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    for (i in 1:P_coefs) {
      clean_name <- gsub("[^A-Za-z0-9]", "_", coef_names[i])
      clean_name <- gsub("_+", "_", clean_name)
      clean_name <- gsub("_$", "", clean_name)

      out_file <- file.path(output_dir, paste0("clbmri_init_", clean_name, ".nii.gz"))

      param_image <- donor_mask * 0
      param_image[is_tissue] <- predicted_parameters[, i]
      antsImageWrite(param_image, out_file)
    }
    message("All initialized maps saved.")
  }

  return(predicted_parameters)
}








