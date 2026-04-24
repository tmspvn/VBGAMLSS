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
#' @param newdata New data to use for predictions. Must be data.frame
#' @param num_cores Number of CPU cores to use for parallel processing.
#'   Defaults to one less than the total available cores if not provided.
#' @param ptype Type of prediction to make: "parameter", "link", "response", "terms". Defaults to "parameter".
#' @param segmentation Optional image/path with region labels.
#' @param segmentation_target Optional. Integer to evaluate (eg 1).
#' @param afold Optional. Integer, fold index when predicting CV folds (for internal use).
#' @param ... Additional arguments passed to the predict.gamlss2 function.
#' @return A structure containing predictions.
#' @import future.apply
#' @import future
#' @export
predict.vbgamlss <- function(object,
                             newdata=NULL,
                             num_cores=NULL,
                             ptype='parameter',
                             segmentation=NULL,
                             segmentation_target=NULL,
                             afold=NULL,
                             terms = NULL,
                             what = NULL,
                             future_plan_strategy = "future.mirai::mirai_cluster",
                             chunk_max_mb = 1024,
                             logdir = NULL,
                             ...){

  if (missing(object)) { stop("vbgamlss object is missing")}
  if (is.null(num_cores)) {num_cores <- future::availableCores()}

  # check segmentation
  if (!is.null(segmentation)){
    if (!is.data.frame(segmentation) && !is.matrix(segmentation)) {
      stop("Error: segmentation must be a data.frame or matrix")
    }
    segmentation <- as.matrix(segmentation)
  }

  # record fam obj if missing
  familyobj <- restore_family(object[[1]])$family
  fname <- familyobj$family

  # compute chunk size
  Nchunks <- estimate_nchunks(object, chunk_max_Mb=chunk_max_mb)

  # predict setup
  if (any(grepl("mirai", future_plan_strategy))) {
    mirai::daemons(num_cores)
    future::plan(strategy=future_plan_strategy)
  } else {
    future::plan(strategy=future_plan_strategy, workers=num_cores)
  }
  options(future.globals.maxSize=50*1024^3) # 10 GB max per prediction

  # get blas omp values
  master_blas <- RhpcBLASctl::blas_get_num_procs()
  master_omp  <- RhpcBLASctl::omp_get_max_threads()

  # split indices to match Core.R chunking
  chunked_indices <- as.list(itertools::isplitIndices(length(object), chunks=Nchunks))

  # PRE-EXTRACT CHUNKS AND BUNDLE DATA
  cat("Preparing data chunks\n")
  prepared_chunks <- lapply(chunked_indices, function(idx_chunk) {
    list(
      raw_bytes = lapply(idx_chunk, function(idx) .subset2(object, idx)),
      seg_data  = if (!is.null(segmentation)) segmentation[, idx_chunk, drop = FALSE] else NULL
    )
  })

  # NUKE THE MASTER OBJECTS TO PREVENT MEMORY LEAKS TO WORKERS
  rm(object)
  if (exists("segmentation")) rm(segmentation)
  gc(verbose = FALSE)

  cat(paste0("Predicting across ", Nchunks, " chunks...\n"))

  # PARALLELIZE OVER THE LIST DIRECTLY
  # By passing `prepared_chunks` as X, future splits it and sends only 1 chunk per worker.
  # PARALLELIZE OVER THE LIST DIRECTLY
  chunk_predictions <- future.apply::future_lapply(prepared_chunks, function(chunk) {

    withCallingHandlers({

      # WORKER THREAD CONTROL
      if (RhpcBLASctl::blas_get_num_procs() > 1L) {RhpcBLASctl::blas_set_num_threads(1L)}
      if (RhpcBLASctl::omp_get_num_procs() > 1L)  {RhpcBLASctl::omp_set_num_threads(1L)}

      # SEQUENTIAL PROCESSING WITHIN THE WORKER
      chunk_res <- lapply(seq_along(chunk$raw_bytes), function(k) {

        raw_data <- chunk$raw_bytes[[k]]
        if (is.null(raw_data)) return(NA)

        vxlgamlss <- qs2::qs_deserialize(raw_data)

        # process only if properly fitted
        if (!isTRUE(vxlgamlss$error) && !is.null(vxlgamlss$vxl)) {

          vxlgamlss$family <- familyobj

          # Subset safely inside the sequential loop
          vxl_newdata <- newdata
          if (!is.null(chunk$seg_data)) {
            vxl_newdata$tissue <- chunk$seg_data[, k]
            if (!is.null(segmentation_target)) {
              vxl_newdata <- vxl_newdata[vxl_newdata$tissue == segmentation_target, , drop = FALSE]
            }
          }

          pred_res <- tryCatch({
            predict(vxlgamlss, newdata = vxl_newdata, type = ptype, terms = terms, what = what, ...)
          }, error = function(e) structure(e$message, class = "try-error"))

          if (inherits(pred_res, "try-error")) {
            message("Prediction failed for voxel ", k, ": ", pred_res)
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
          return(l)

        } else {
          return(NA)
        }
      })

      return(chunk_res)

    }, error = function(e) {
      # HARDCODED PATH FOR TRACEBACK
      cat(sprintf("\n=== ERROR CAUGHT [%s] ===\n%s\n--- TRACEBACK ---\n",
                  Sys.time(), conditionMessage(e)),
          file = log_path, append = TRUE)

      capture.output(print(sys.calls()),
                     file = file.path(logdir,
                                      paste0(rand_names(1), '.vbgamlss_predict.error')),
                     append = TRUE)
    })

  }, future.seed = TRUE)

  # Reverse blas and openmp threads control (Master session)
  if (RhpcBLASctl::blas_get_num_procs() != master_blas)
    {RhpcBLASctl::blas_set_num_threads(master_blas)}
  if (RhpcBLASctl::omp_get_num_procs() != master_omp)
    {RhpcBLASctl::omp_set_num_threads(master_omp)}

  # 4. FLATTEN RESULTS
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



#' Apply kernel james-stein shrinkage to intercept coefficients of a vbgamlss
#'
#' @export
james_stein <- function(object,
                        mask,
                        yimageframe,
                        save_model=NULL,
                        radius_voxels=4,
                        num_cores=NULL){

  # Internal
  JS <- function(theta_hat, sigma2) {
    N <- length(theta_hat)
    theta_bar <- mean(theta_hat)
    # SSD sum of squared deviations
    V <- sum((theta_hat - theta_bar)^2)
    # c+
    c_plus <- pmax(0, 1 - ((N - 2) * sigma2) / V)
    # Calculate the shrunken estimates
    theta_js <- theta_bar + c_plus * (theta_hat - theta_bar)
    return(theta_js)
  }

  # checks
  if (is.character(mask)) {
    if (!file.exists(mask)) stop("File does not exist: ", mask)
    mask <- antsImageRead(mask)
  }

  # Matrix conversion for much faster column subsetting
  if (!is.data.frame(yimageframe) && !is.matrix(yimageframe)) {
    stop("Error: yimageframe must be a data.frame or matrix")
  }

  yimageframe <- as.matrix(yimageframe)
  mask_array <- as.array(mask)
  vox_coords <- which(mask_array > 0, arr.ind = TRUE)

  # parallel setup
  if (is.null(num_cores)) {num_cores <- future::availableCores()}
  future::plan(strategy="future::cluster", workers=num_cores)
  options(future.globals.maxSize=10*1024^3) # 10 GB max

  # Finds 'k' neighbors
  neighbors <- dbscan::kNN(vox_coords, k = round(4.188 *(radius_voxels^3)))$id

  # ----------------------------------------------------------------------------
  # Must loop like this otherwise it copies the whole model a ton of times
  nvox <- dim(neighbors)[1]
  all_voxels_js_coefs <- vector("list", nvox)

  # Per voxel do:
  for (vxl in 1:nvox) {

    kernel_idx <- neighbors[vxl,]
    # Extract ONLY the raw bytes for the subset of voxel for the kernel
    kernel_mods <- lapply(kernel_idx, function(idx) .subset2(object, idx))
    kernel_y <- yimageframe[, kernel_idx]

    # parallelise extraction
    kernel_coefs <- future.apply::future_lapply(seq_along(kernel_mods), function(k) {

      # Grab the raw bytes for this specific worker
      raw_data <- kernel_mods[[k]]
      # Explicitly deserialize inside the worker
      if (is.null(raw_data)) return(NA)
      vxlgamlss <- qs2::qs_deserialize(raw_data)

      # get error
      e <- kernel_y[, k] - predict(vxlgamlss, type = "response")

      # extract coefficients
      vxlcoefs <- unlist(coef(vxlgamlss))
      vxlcoefs['var.err'] <-  var(e)

      # return coefficients
      vxlcoefs

    })

    # bind output
    kernel_coefs <- as.data.frame(do.call(rbind, kernel_coefs))

    # compute James-stein shrinkage
    sigma2_noise <- kernel_coefs[['var.err']]
    kernel_coefs[["var.err"]] <- NULL

    vxl_js_matrix <- apply(kernel_coefs, 2, function(x) JS(x, sigma2_noise))
    vxl_js <- vxl_js_matrix[1, ]

    # store
    all_voxels_js_coefs[[vxl]] <- vxl_js

    cat("Voxel", vxl, 'of', nvox, fill = TRUE)
    flush.console()
  }

  # bind output again
  all_voxels_js_coefs <- as.data.frame(do.call(rbind, all_voxels_js_coefs))
  coefs_names <- names(all_voxels_js_coefs)

  # ----------------------------------------------------------------------------
  # Loop again to assign the new shrieked coefficients
  for (vxl in 1:nvox) {

    # Subset and deserialize
    new_coef <- all_voxels_js_coefs[vxl,]
    raw_data <- .subset2(object, vxl)
    if (is.null(raw_data)) next

    vxlgamlss <- qs2::qs_deserialize(raw_data)

    # Map values back to the internal vxlgamlss structure
    for (par in names(vxlgamlss$coefficients)) {
      # Identify indices belonging to this parameter
      pattern <- paste0("^", par, "\\.")
      par_vals <- new_coef[grep(pattern, names(new_coef))]
      # Remove the "mu." prefix so names match internal (Intercept), sexM, etc.
      names(par_vals) <- gsub(pattern, "", names(par_vals))
      # Assign back to the model object
      vxlgamlss$coefficients[[par]] <- par_vals
    }

    # reserialize and replace
    object[[vxl]] <- qs2::qs_serialize(vxlgamlss)

    cat("Voxel", vxl, 'of', nvox, fill = TRUE)
    flush.console()
  }

  # Done save model back
  if (!is.null(save_model)) {
    qs2::qs_save(object, # FIXED: Save 'object', not undefined 'models'
                 file = paste0(save_model, '.JS.vbgamlss'),
                 compress_level = 0L) # uncompressed
    cat('Model saved: ', paste0(save_model, '.JS.vbgamlss'), '\n')
    return(NULL)
  } else {
    return(object)
  }
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
#'   familyobj <- restore_family(object[[1]])$family
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

