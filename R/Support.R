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





# ------------------------------------------------------------------------------
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
                             ...){

  if (missing(object)) { stop("vbgamlss object is missing")}
  if (is.null(num_cores)) {num_cores <- future::availableCores()}

  # check segmentation
  if (!is.null(segmentation)){
    if (!is.data.frame(segmentation) && !is.matrix(segmentation)) {
      stop("Error: segmentation must be a data.frame or matrix")
    }
  }

  # subset the dataframe if the input is a fold from CV
  #if (!is.null(afold)){
  #  if (!is.logical(afold) && !is.integer(afold)) {
  #    stop("Error: afold must be either logical or integer vector.")
  #  }
  #  if (!is.null(segmentation)) segmentation <- segmentation[afold, , drop=FALSE]
  #}

  # record fam obj if missing
  familyobj <- restore_family(object[[1]])$family
  fname <- familyobj$family

  # compute chunk size
  Nchunks <- estimate_nchunks(object)

  # predict setup
  future::plan(strategy="future::cluster", workers=num_cores)
  options(future.globals.maxSize=10*1024^3) # 10 GB max per prediction

  # split indices to match Core.R chunking
  chunked_indices <- as.list(itertools::isplitIndices(length(object), chunks=Nchunks))

  # Pre-allocate the predictions list
  predictions <- vector("list", length(object))

  # chunks loop
  for (i in seq_along(chunked_indices)){
    idx_chunk <- chunked_indices[[i]]
    cat(paste0("Chunk: ", i, "/", Nchunks, " (Predicting)\n"))

    # Extract ONLY the raw bytes for this chunk in the master session.
    # .subset2 bypasses the S3 method, keeping them as lightweight raw bytes.
    chunk_raw_bytes <- lapply(idx_chunk, function(idx) .subset2(object, idx))

    # parallel call using future_lapply
    subpr <- future.apply::future_lapply(seq_along(idx_chunk), function(k) {

      # Grab the raw bytes for this specific worker
      raw_data <- chunk_raw_bytes[[k]]

      # Explicitly deserialize inside the worker, bypasses the S3 inheritance problem.
      if (is.null(raw_data)) return(NA)
      vxlgamlss <- qs2::qs_deserialize(raw_data)

      # process only if properly fitted (no error flag)
      if (!isTRUE(vxlgamlss$error) && !is.null(vxlgamlss$vxl)) {

        # recon family object
        vxlgamlss$family <- familyobj

        # if multi tissue add
        if (!is.null(segmentation)){
          # Note: k is the index within the chunk. To get the true voxel index for segmentation:
          true_idx <- idx_chunk[k]
          newdata$tissue <- segmentation[, true_idx]
          if (!is.null(segmentation_target)){
            newdata <- newdata[newdata$tissue == segmentation_target, , drop=FALSE]
          }}

        # predict
        pred_res <- tryCatch({
          predict(vxlgamlss,
                  newdata = newdata,
                  type = ptype,
                  terms = terms,
                  what = what,
                  ...)
        }, error = function(e) {
          # Return an object of class "try-error"
          structure(e$message, class = "try-error")
        })

        # If the prediction failed, return NA immediately for this voxel
        if (inherits(pred_res, "try-error")) {
          return(NA)
        }

        # Structure the list correctly
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
    }, future.seed = TRUE)

    # Put the chunk results into the pre-allocated list
    predictions[idx_chunk] <- subpr
    gc(verbose = F)
  }

  gc(verbose = FALSE)
  return(structure(predictions, class = "vbgamlss.predictions"))
}






# ------------------------------------------------------------------------------
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

# predict.vbgamlss <- function(object,
#                              newdata=NULL,
#                              num_cores=NULL,
#                              ptype='parameter',
#                              segmentation=NULL,
#                              segmentation_target=NULL,
#                              afold=NULL,
#                              terms = NULL,
#                              what = NULL,
#                              ...){
#
#   if (missing(object)) { stop("vbgamlss object is missing")}
#   if (is.null(num_cores)) {num_cores <- future::availableCores()}
#
#   # check segmentation
#   if (!is.null(segmentation)){
#     if (!is.data.frame(segmentation) && !is.matrix(segmentation)) {
#       stop("Error: segmentation must be a data.frame or matrix")
#     }
#   }
#
#   # subset the dataframe if the input is a fold from CV
#   if (!is.null(afold)){
#     if (!is.logical(afold) && !is.integer(afold)) {
#       stop("Error: afold must be either logical or integer vector.")
#     }
#     if (!is.null(segmentation)) segmentation <- segmentation[afold, , drop=FALSE]
#   }
#
#   fname <- as.character(object[[1]]$family)
#   familyobj <- gamlss2:::complete_family(get(fname))
#
#   # compute chunk size
#   Nchunks <- estimate_nchunks(object)
#
#   # predict setup
#   future::plan(strategy="future::cluster", workers=num_cores)
#
#   # split indices to match Core.R chunking
#   chunked_indices <- as.list(itertools::isplitIndices(length(object), chunks=Nchunks))
#   predictions <- list()
#
#   # chunks loop
#   for (i in seq_along(chunked_indices)){
#     idx_chunk <- chunked_indices[[i]]
#     cat(paste0("Chunk: ", i, "/", Nchunks, " (Predicting)"), fill=TRUE)
#
#     # parallel call using future.apply to match Core.R framework
#     subpr <- future.apply::future_lapply(idx_chunk, function(idx) {
#
#       vxlgamlss <- object[[idx]]
#
#       # process only if properly fitted (no error flag)
#       if (!isTRUE(vxlgamlss$error) && !is.null(vxlgamlss$vxl)) {
#
#         # recon family object (Core.R strips it to string to save memory)
#         vxlgamlss$family <- familyobj
#
#         # if multi tissue add
#         if (!is.null(segmentation)){
#           newdata$tissue <- segmentation[, idx]
#           if (!is.null(segmentation_target)){
#             newdata <- newdata[newdata$tissue == segmentation_target, , drop=FALSE]
#           }}
#
#         # predict
#         pred_res <- predict(vxlgamlss,
#                             newdata = newdata,
#                             type = ptype,
#                             terms = terms,
#                             what = what,
#                             ...)
#
#         # Structure the list correctly
#         if (!is.null(what) && length(what) == 1) {
#           l <- list()
#           l[[what]] <- as.numeric(pred_res)
#         } else {
#           l <- as.list(pred_res)
#         }
#
#         l$family <- fname
#         l$vxl <- vxlgamlss$vxl
#
#         return(l)
#
#       } else {
#         return(NA)
#       }
#     }, future.seed = TRUE)
#
#     predictions <- c(predictions, subpr)
#   }
#
#   gc()
#   return(structure(predictions, class = "vbgamlss.predictions"))
# }
#
