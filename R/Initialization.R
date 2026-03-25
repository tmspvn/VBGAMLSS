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

  # Return the fitted values on the response scale
  return(fitted(glob_fit))
}












