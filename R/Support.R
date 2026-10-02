################################  Support  #####################################








# ------------------------------------------------------------------------------
#' Save VBGAMLSS models from files with the specified prefix.
#'
#' @param model_list vbgamlss fitted model.
#' @param filename Prefix of the files to save the VBGAMLSS models.
#' @param voxel Index/s of the Voxel/s to save model for. If NULL, save the entire model.
#' @export
save_model <- function(model_list, filename, voxel=NULL) {
  # save
  if (is.null(voxel)) {
    qs2::qs_save(model_list, file = glue(filename, ".vbgamlss"))
  } else {
    warning("Untested, possibly not working as intended!")
    qs2::qs_save(model_list[[voxel]], file = glue(filename, ".voxel{voxel}.vbgamlss"))
  }
  cat('Saved VBGAMLSS model to:', filename)
}





# ------------------------------------------------------------------------------
#' Load VBGAMLSS models from files with the specified prefix.
#'
#' @param filepath path to serialized model.
#' @return A structure containing the loaded VBGAMLSS models.
#' @export
load_model <- function(filepath) {
  # load
  return(structure(qs2::qs_read(filepath), class = "vbgamlss"))
}





#' Predict method for vbgamlss objects
#'
#' @param object vbgamlss object to make predictions from.
#' @param newdata New data to use for predictions. Must be data.frame.
#' @param num_cores Number of CPU cores to use for parallel processing.
#'   Defaults to available cores if not provided.
#' @param ptype Type of prediction to make: "parameter", "link", "response", "terms".
#'   Defaults to "parameter".
#' @param segmentation Optional data.frame or matrix with region labels (e.g., tissue types)
#'   mapping to the voxels in the object.
#' @param segmentation_target Optional integer. If provided, the prediction is only
#'   performed for voxels where the segmentation matches this value.
#' @param terms Optional character vector. What terms to include in prediction. See gamlss2 doc.
#' @param what Optional character. Specifies which distribution parameter to predict
#'   (e.g., "mu", "sigma").
#' @param future_plan_strategy Character. The parallel back-end strategy to use.
#'   Defaults to "future.mirai::mirai_cluster".
#' @param chunk_max_mb Integer. Maximum size in MB for each data chunk sent to workers.
#'   Defaults to 1024.
#' @param ... Additional arguments passed to the underlying predict.gamlss2 function.
#'
#' @return A structure of class "vbgamlss.predictions" containing a list of
#'   prediction results for each voxel.
#'
#' @import future.apply
#' @import future
#' @importFrom qs2 qs_save qs_read qs_deserialize
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads omp_get_max_threads omp_set_num_threads
#' @importFrom itertools isplitIndices
#' @export
predict.vbgamlss <- function(object,
                             newdata=NULL,
                             num_cores=NULL,
                             ptype='parameter',
                             segmentation=NULL,
                             segmentation_target=NULL,
                             terms = NULL,
                             what = NULL,
                             future_plan_strategy = "future.mirai::mirai_cluster",
                             chunk_max_mb = 1024,
                             se.fit=FALSE,
                             ...) {

  # Checks
  if (missing(object)) stop("vbgamlss object is missing")
  if (is.null(num_cores)) num_cores <- future::availableCores()

  if (!is.null(segmentation)) {
    if (!is.data.frame(segmentation) && !is.matrix(segmentation)) {
      stop("Error: segmentation must be a data.frame or matrix")
    }
    segmentation <- as.matrix(segmentation)
  }

  # Extract family object for workers
  familyobj <- restore_family(first_fitted(object))$family
  fname <- familyobj$family

  # Prepare for parallelization
  if (any(grepl("mirai", future_plan_strategy))) {
    mirai::daemons(num_cores)
    future::plan(strategy=future_plan_strategy)
  } else {
    future::plan(strategy=future_plan_strategy, workers=num_cores)
  }
  options(future.globals.maxSize=10*1024^3) # max 10GB otherwise it crashes

  # Chunk
  Nchunks <- estimate_nchunks(object, chunk_max_Mb=chunk_max_mb)
  chunked_indices <- as.list(itertools::isplitIndices( length(object), chunks=Nchunks))

  # Cache chunk for parallel processing
  prepared_chunks <- lapply(seq_along(chunked_indices), function(i) {
    idx_chunk <- chunked_indices[[i]]
    raw_bytes_list <- lapply(idx_chunk, function(idx) .subset2(object, idx))

    # Generate unique temp file
    chunk_path <- tempfile(pattern = sprintf("vbgamlss_%s_predchunk",
                                             rand_names(1)),
                           tmpdir = tempdir(),
                           fileext = ".qs")
    qs2::qs_save(raw_bytes_list, file = chunk_path)

    list(
      file_path = chunk_path,
      seg_data  = if (!is.null(segmentation)) segmentation[, idx_chunk, drop = FALSE] else NULL
    )
  })

  # Cleanup \tmp on exit
  on.exit({
    unlink(sapply(prepared_chunks, `[[`, "file_path"))
  }, add = TRUE)

  # remove large objects that are no longer needed
  rm(object)
  if (exists("segmentation")) rm(segmentation)
  gc(verbose = FALSE)


  # Actual parallel call
  chunk_predictions <- future.apply::future_lapply(prepared_chunks, function(chunk) {

    # avoid competing workers
    if (RhpcBLASctl::blas_get_num_procs() > 1L) RhpcBLASctl::blas_set_num_threads(1L)
    if (RhpcBLASctl::omp_get_num_procs() > 1L)  RhpcBLASctl::omp_set_num_threads(1L)

    # Load one temporary chunk
    chunk_raw_bytes <- qs2::qs_read(chunk$file_path)

    # process the chunk sequentially
    chunk_res <- lapply(seq_along(chunk_raw_bytes), function(k) {

      # get single model and deserialize
      raw_data <- chunk_raw_bytes[[k]]
      if (is.null(raw_data)) return(NA)
      vxlgamlss <- qs2::qs_deserialize(raw_data)

      if (!isTRUE(vxlgamlss$error) && !is.null(vxlgamlss$vxl)) {

        # segmentation?
        vxlgamlss$family <- familyobj
        vxl_newdata <- newdata
        vxl_degfre <- vxlgamlss$df

        if (!is.null(chunk$seg_data)) {
          vxl_newdata$tissue <- chunk$seg_data[, k]
          if (!is.null(segmentation_target)) {
            vxl_newdata <- vxl_newdata[vxl_newdata$tissue == segmentation_target, , drop = FALSE]
          }
        }

        # Actual predict call
        pred_res <- tryCatch({
          predict(vxlgamlss, newdata = vxl_newdata, type = ptype,
                  terms = terms, what = what, se.fit = se.fit, ...)
        }, error = function(e) structure(e$message, class = "try-error"))

        if (inherits(pred_res, "try-error")) {
          message("Prediction failed for voxel: ", pred_res)
          return(NA)
        }

        if (!is.null(what) && length(what) == 1) {
          l <- list()
          l[[what]] <- as.numeric(pred_res)
        } else {
          l <- as.list(pred_res)
        }

        l$family <- fname
        l$vxl <- vxlgamlss$vxl
        l$df <- vxl_degfre
        return(l)

      } else {
        return(NA)
      }
    })

    return(chunk_res)

  }, future.seed = TRUE)

  # Return
  predictions <- unlist(chunk_predictions, recursive = FALSE)
  gc(verbose = FALSE)
  return(structure(predictions, class = "vbgamlss.predictions"))
}






# ------------------------------------
#' Compute Z-scores for vbgamlss predictions given Y voxel data and image mask.
#' @param predictions A vbgamlss.predictions object.
#' @param yimageframe A data.frame or matrix of observed response values.
#' @param num_cores Number of CPU cores to use for parallel processing.
#' @return A structure containing z-scores.
#' @import future.apply
#' @import future
#' @export
zscore.vbgamlss <- function(predictions, yimageframe, num_cores=NULL){

  if (missing(predictions)) { stop("vbgamlss.predictions is missing")}
  if (missing(yimageframe)) { stop("yimageframe is missing")}
  if (is.null(num_cores)) {num_cores <- future::availableCores()}

  # Matrix conversion for much faster column subsetting
  if (!is.data.frame(yimageframe) && !is.matrix(yimageframe)) {
    stop("Error: yimageframe must be a data.frame or matrix")
  }
  yimageframe <- as.matrix(yimageframe)

  # parallel function
  do.zscore <- function(obj) {
    pred <- obj$pred
    yval <- obj$yvxldat

    # Gracefully handle failed predictions (NA) from the fitting phase
    if (length(pred) == 1 && is.na(pred[1])) {
      return(rep(NA, length(yval)))
    }

    # get number of params
    lpar <- sum(names(pred) %in% c("mu", "sigma", "nu", "tau"))
    qfun <- paste("p", pred$family, sep="")

    # explicit 'q' argument prevents positional mismatches
    if (lpar == 1) {
      newcall <- call(qfun, q = yval, mu = pred$mu)
    } else if (lpar == 2) {
      newcall <- call(qfun, q = yval, mu = pred$mu, sigma = pred$sigma)
    } else if (lpar == 3) {
      newcall <- call(qfun, q = yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu)
    } else {
      newcall <- call(qfun, q = yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu, tau = pred$tau)
    }

    cdf <- eval(newcall)

    # Bound CDF to prevent Inf/-Inf z-scores from perfect 0 or 1 probabilities
    cdf[cdf < 1e-12] <- 1e-12
    cdf[cdf > (1 - 1e-12)] <- 1 - 1e-12

    rqres <- qnorm(cdf)
    return(rqres)
  }

  # compute chunk size
  Nchunks <- estimate_nchunks(yimageframe)

  # predict setup
  future::plan(strategy="future::cluster", workers=num_cores)

  zscores <- list()
  chunked_indices <- as.list(itertools::isplitIndices(ncol(yimageframe), chunks=Nchunks))

  for (i in seq_along(chunked_indices)){
    ichunk <- chunked_indices[[i]]
    cat(paste0("Chunk: ", i, "/", Nchunks, " (Computing Z-scores)\n"))

    # Extract chunk data
    y_chunk <- yimageframe[, ichunk, drop=FALSE]
    pred_chunk <- predictions[ichunk]

    # Pack iterable chunk for workers
    iterable_chunk <- lapply(seq_along(ichunk), function(j) {
      list(pred = pred_chunk[[j]], yvxldat = as.numeric(y_chunk[, j]))
    })

    # compute z-scores in parallel using future framework
    subzs <- future.apply::future_lapply(iterable_chunk,
                                         do.zscore,
                                         future.seed = TRUE,
                                         future.packages = c("gamlss.dist",
                                                             "gamlss",
                                                             "gamlss2"))
    zscores <- c(zscores, subzs)
  }

  gc()
  return(structure(zscores, class = "vbgamlss.zscores"))
}














# ------------------------------------
#' Empirical-Bayes spatial shrinkage of voxel-wise linear coefficients (post hoc)
#'
#' Takes a fitted vbgamlss model and returns a shrunk copy; the input is left untouched.
#' For every voxel v and linear coefficient (e.g. mu intercept, sexM) with estimate b_v and
#' standard error s_v (stored at fit time by vbgamlss, `$se_linear`), a local prior g_v is fitted
#' with ebnm (marginal maximum likelihood, mode estimated) on the OTHER mask voxels in a
#' kernel_size^3 cube around v, and b_v is replaced by its posterior mean under g_v. Noisy voxels
#' are pulled towards their neighbourhood; with an adaptive prior (default "point_laplace") large,
#' well-measured departures (e.g. at tract boundaries) are kept. Smooth (e.g. pb) and random-effect
#' terms are left as fitted.
#'
#' @param object A fitted vbgamlss model.
#' @param mask The mask used to build the imageframe (path, antsImage or 3D array); it defines the voxel order.
#' @param params Distribution parameters whose linear coefficients are shrunk. Default "mu"; shrinking
#'   sigma/nu/tau changes the predicted spread, so re-check z-score calibration.
#' @param prior_family Any ebnm prior family: "point_laplace" (default, adaptive, ~16 ms per voxel and
#'   coefficient), "normal" (classic normal-normal EB, ~7 ms), "normal_scale_mixture" (ash, ~130 ms), ...
#' @param kernel_size Odd side of the cubic kernel in voxels (5 = 5x5x5).
#' @param min_neighbors Voxels with fewer valid neighbours in the kernel are left unshrunk.
#' @param num_cores Cores (forked) for the per-voxel prior fits and for re-writing the voxel models.
#' @param save_model Optional path; saved as <save_model>.vbgamlss.
#' @return The shrunk vbgamlss copy. Each changed voxel model gains `$eb_shrink` (original, posterior_sd,
#'   prior_mode per coefficient); attribute "eb_shrink" summarises the shift towards the prior mode.
#' @export
eb_shrink <- function(object, mask, params = "mu", prior_family = "point_laplace", kernel_size = 5,
                      min_neighbors = 10, num_cores = 1L, save_model = NULL) {

  if (!requireNamespace("ebnm", quietly = TRUE)) stop("eb_shrink needs the ebnm package: install.packages('ebnm')")
  if (kernel_size %% 2 != 1) stop("kernel_size must be odd")
  if (is.character(mask)) mask <- ANTsR::antsImageRead(mask)
  mask_arr <- as.array(mask) > 0
  coords <- which(mask_arr, arr.ind = TRUE) # same voxel order as images2matrix()
  nvox <- length(object)
  if (nrow(coords) != nvox) stop("mask has ", nrow(coords), " voxels, the model has ", nvox)

  # 1. estimates and standard errors (one read per voxel)
  cat("Reading coefficients and standard errors\n")
  est <- pbmcapply::pbmclapply(seq_len(nvox), function(i) {
    g <- object[[i]]
    if (is.null(g) || isTRUE(g$error)) return(NULL)
    if (is.null(g$se_linear)) return(NULL)
    b  <- unlist(g$coefficients[params]) # "mu.(Intercept)", ...
    se <- g$se_linear
    names(se) <- sub(".p.", ".", names(se), fixed = TRUE) # "mu.p.(Intercept)" -> "mu.(Intercept)"
    list(b = b, s = se[names(b)])
  }, mc.cores = num_cores)
  if (all(vapply(est, is.null, logical(1)))) {
    stop("No voxel has stored standard errors ($se_linear): refit the model with this VBGAMLSS version")
  }
  cn <- unique(unlist(lapply(est, function(e) names(e$b))))
  B <- S <- matrix(NA_real_, nvox, length(cn), dimnames = list(NULL, cn))
  for (i in which(!vapply(est, is.null, logical(1)))) {
    B[i, names(est[[i]]$b)] <- est[[i]]$b
    S[i, names(est[[i]]$s)] <- est[[i]]$s
  }
  bad <- !is.finite(S) | S <= 0
  B[bad] <- NA; S[bad] <- NA
  rm(est)

  # 2. neighbours: the other mask voxels in the cube (NA outside the mask)
  h <- (kernel_size - 1) / 2
  offs <- as.matrix(expand.grid(-h:h, -h:h, -h:h))
  offs <- offs[rowSums(abs(offs)) > 0, , drop = FALSE]
  lookup <- array(0L, dim(mask_arr)); lookup[mask_arr] <- seq_len(nvox)
  d <- dim(mask_arr)
  nb <- matrix(vapply(seq_len(nrow(offs)), function(k) {
    cc <- sweep(coords, 2, offs[k, ], "+")
    inside <- cc[, 1] >= 1 & cc[, 1] <= d[1] & cc[, 2] >= 1 & cc[, 2] <= d[2] & cc[, 3] >= 1 & cc[, 3] <= d[3]
    out <- rep(NA_integer_, nvox)
    out[inside] <- lookup[cc[inside, , drop = FALSE]]
    out[out == 0L] <- NA_integer_
    out
  }, integer(nvox)), nrow = nvox)

  # 3. local empirical-Bayes shrinkage with ebnm: prior fitted on the neighbours, posterior for the voxel
  cat("Fitting local", prior_family, "priors\n")
  post <- pbmcapply::pbmclapply(seq_len(nvox), function(v) {
    res <- matrix(NA_real_, 3, length(cn), dimnames = list(c("mean", "sd", "mode"), cn))
    for (j in cn) {
      if (is.na(B[v, j])) next
      idx <- nb[v, ]
      idx <- idx[!is.na(idx) & !is.na(B[idx, j])]
      if (length(idx) < min_neighbors) next
      res[, j] <- tryCatch({
        g <- ebnm::ebnm(B[idx, j], S[idx, j], prior_family = prior_family, mode = "estimate",
                        output = "fitted_g")$fitted_g
        p <- ebnm::ebnm(B[v, j], S[v, j], prior_family = prior_family, g_init = g, fix_g = TRUE,
                        output = c("posterior_mean", "posterior_sd"))$posterior
        c(p$mean, p$sd, g$mean[1])
      }, error = function(e) rep(NA_real_, 3)) # failed prior fit: voxel left unshrunk
    }
    res
  }, mc.cores = num_cores)
  newB <- B
  PSD <- MODE <- B * NA
  for (v in seq_len(nvox)) {
    ok <- !is.na(post[[v]]["mean", ])
    newB[v, ok] <- post[[v]]["mean", ok]; PSD[v, ok] <- post[[v]]["sd", ok]; MODE[v, ok] <- post[[v]]["mode", ok]
  }
  rm(post)
  L <- 1 - (newB - MODE) / (B - MODE) # fraction of the way to the prior mode (0 = kept, 1 = fully shrunk)

  # 4. write the shrunk coefficients into a copy of the model
  changed <- which(rowSums(!is.na(PSD)) > 0)
  cat("Writing", length(changed), "of", nvox, "voxel models\n")
  raws <- pbmcapply::pbmclapply(changed, function(i) {
    g <- qs2::qs_deserialize(.subset2(object, i)) # stored form, no restore_family
    for (p in params) {
      nm <- names(g$coefficients[[p]])
      v  <- newB[i, paste0(p, ".", nm)]
      g$coefficients[[p]][nm[!is.na(v)]] <- v[!is.na(v)]
    }
    g$eb_shrink <- list(original = B[i, ], posterior_sd = PSD[i, ], prior_mode = MODE[i, ], prior_family = prior_family)
    qs2::qs_serialize(g)
  }, mc.cores = num_cores)
  out <- object
  for (k in seq_along(changed)) out[[changed[k]]] <- raws[[k]]
  attr(out, "eb_shrink") <- list(params = params, prior_family = prior_family, kernel_size = kernel_size,
                                 min_neighbors = min_neighbors, n_shrunk = length(changed),
                                 shrink_fraction = apply(L, 2, stats::quantile, c(.05, .25, .5, .75, .95), na.rm = TRUE))
  cat("Fraction of the way to the local prior mode (0 = kept, 1 = fully shrunk):\n")
  print(round(attr(out, "eb_shrink")$shrink_fraction, 3))

  if (!is.null(save_model)) {
    qs2::qs_save(out, file = paste0(save_model, ".vbgamlss"), compress_level = 0L)
    cat("Model saved: ", paste0(save_model, ".vbgamlss"), "\n")
  }
  out
}






























############################## ========== ######################################
                             # DEPRECATED #
############################## ========== ######################################

#' # ------------------------------------------------------------------------------
#' #' Predict method for vbgamlss objects
#' #'
#' #' @param object vbgamlss object to make predictions from.
#' #' @param newdata New data to use for predictions. Must be data.frame
#' #' @param num_cores Number of CPU cores to use for parallel processing.
#' #'   Defaults to one less than the total available cores if not provided.
#' #' @param ptype Type of prediction to make: "parameter", "link", "response", "terms". Defaults to "parameter".
#' #' @param segmentation Optional image/path with region labels.
#' #' @param segmentation_target Optional. Integer to evaluate (eg 1).
#' #' @param afold Optional. Integer, fold index when predicting CV folds (for internal use).
#' #' @param ... Additional arguments passed to the predict.gamlss2 function.
#' #' @return A structure containing predictions.
#' #' @import future.apply
#' #' @import future
#' #' @export
#' predict.vbgamlss <- function(object,
#'                              newdata=NULL,
#'                              num_cores=NULL,
#'                              ptype='parameter',
#'                              segmentation=NULL,
#'                              segmentation_target=NULL,
#'                              afold=NULL,
#'                              terms = NULL,
#'                              what = NULL,
#'                              future_plan_strategy = "future.mirai::mirai_cluster",
#'                              chunk_max_mb = 1024,
#'                              ...){
#'
#'   if (missing(object)) { stop("vbgamlss object is missing")}
#'   if (is.null(num_cores)) {num_cores <- future::availableCores()}
#'
#'   # check segmentation
#'   if (!is.null(segmentation)){
#'     if (!is.data.frame(segmentation) && !is.matrix(segmentation)) {
#'       stop("Error: segmentation must be a data.frame or matrix")
#'     }
#'     segmentation <- as.matrix(segmentation)
#'   }
#'
#'   # record fam obj if missing
#'   familyobj <- restore_family(first_fitted(object))$family
#'   fname <- familyobj$family
#'
#'   # compute chunk size
#'   Nchunks <- estimate_nchunks(object, chunk_max_Mb=chunk_max_mb)
#'
#'   # predict setup
#'   if (any(grepl("mirai", future_plan_strategy))) {
#'     mirai::daemons(num_cores)
#'     future::plan(strategy=future_plan_strategy)
#'   } else {
#'     future::plan(strategy=future_plan_strategy, workers=num_cores)
#'   }
#'   options(future.globals.maxSize=10*1024^3) # 10 GB max per prediction
#'
#'   # get blas omp values
#'   master_blas <- RhpcBLASctl::blas_get_num_procs()
#'   master_omp  <- RhpcBLASctl::omp_get_max_threads()
#'
#'   # split indices to match Core.R chunking
#'   chunked_indices <- as.list(itertools::isplitIndices(length(object), chunks=Nchunks))
#'
#'   # Pre-allocate the predictions list
#'   predictions <- vector("list", length(object))
#'
#'   # chunks loop
#'   for (i in seq_along(chunked_indices)){
#'     idx_chunk <- chunked_indices[[i]]
#'     cat(paste0("Chunk: ", i, "/", Nchunks, " (Predicting)\n"))
#'
#'     # 1. PREPROCESSING OUTSIDE WORKERS (Matching vbgamlss dispatch logic)
#'     if (!is.null(segmentation)) {
#'       voxelseg_chunked <- segmentation[, idx_chunk, drop = FALSE]
#'     }
#'
#'     prep_func <- function(k) {
#'       true_idx <- idx_chunk[k]
#'       vxl_newdata <- newdata
#'
#'       if (!is.null(segmentation)) {
#'         vxl_newdata$tissue <- voxelseg_chunked[, k]
#'         if (!is.null(segmentation_target)) {
#'           vxl_newdata <- vxl_newdata[vxl_newdata$tissue == segmentation_target, , drop = FALSE]
#'         }
#'       }
#'
#'       list(
#'         raw_data = .subset2(object, true_idx),
#'         newdata  = vxl_newdata
#'       )
#'     }
#'
#'     prepared_voxel_data <- lapply(seq_along(idx_chunk), prep_func)
#'
#'     # Clean up to free memory before parallel execution
#'     if (!is.null(segmentation)) rm(voxelseg_chunked)
#'     gc(verbose = FALSE)
#'
#'
#'     # 2. PARALLEL EXECUTION
#'     subpr <- future.apply::future_lapply(prepared_voxel_data, function(vxl_item) {
#'
#'       # WORKER THREAD CONTROL
#'       if (RhpcBLASctl::blas_get_num_procs() > 1L)
#'       {RhpcBLASctl::blas_set_num_threads(1L)}
#'       if (RhpcBLASctl::omp_get_num_procs() > 1L)
#'       {RhpcBLASctl::omp_set_num_threads(1L)}
#'
#'       # Explicitly deserialize inside the worker
#'       if (is.null(vxl_item$raw_data)) return(NA)
#'       vxlgamlss <- qs2::qs_deserialize(vxl_item$raw_data)
#'
#'       # process only if properly fitted (no error flag)
#'       if (!isTRUE(vxlgamlss$error) && !is.null(vxlgamlss$vxl)) {
#'
#'         # recon family object
#'         vxlgamlss$family <- familyobj
#'
#'         # predict using the pre-assembled voxel-specific newdata
#'         pred_res <- tryCatch({
#'           predict(vxlgamlss,
#'                   newdata = vxl_item$newdata,
#'                   type = ptype,
#'                   terms = terms,
#'                   what = what,
#'                   ...)
#'         }, error = function(e) {
#'           structure(e$message, class = "try-error")
#'         })
#'
#'         # If prediction failed, return NA
#'         if (inherits(pred_res, "try-error")) {
#'           return(NA)
#'         }
#'
#'         # Structure the list correctly
#'         if (!is.null(what) && length(what) == 1) {
#'           l <- list()
#'           l[[what]] <- as.numeric(pred_res)
#'         } else {
#'           l <- as.list(pred_res)
#'         }
#'
#'         l$family <- fname
#'         l$vxl <- vxlgamlss$vxl
#'
#'         return(l)
#'
#'       } else {
#'         return(NA)
#'       }
#'     }, future.seed = TRUE)
#'
#'     # Put chunk results into the pre-allocated list
#'     predictions[idx_chunk] <- subpr
#'     gc(verbose = FALSE)
#'
#'     # Reverse blas and openmp threads control
#'     if (RhpcBLASctl::blas_get_num_procs() != master_blas)
#'     {RhpcBLASctl::blas_set_num_threads(master_blas)}
#'     if (RhpcBLASctl::omp_get_num_procs() != master_omp)
#'     {RhpcBLASctl::omp_set_num_threads(master_omp)}
#'
#'   }
#'
#'   gc(verbose = FALSE)
#'   return(structure(predictions, class = "vbgamlss.predictions"))
#' }
#'
#'
#'

