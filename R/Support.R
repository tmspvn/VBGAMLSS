################################  Support  #####################################








#' Save VBGAMLSS models from files with the specified prefix.
#'
#' @param model_list vbgamlss fitted model.
#' @param filename Prefix of the files to save the VBGAMLSS models.
#' @param voxel Index/s of the Voxel/s to save model for. If NULL, save the entire model.
#' @export
save_model <- function(model_list, filename, voxel=NULL) {
  # save
  if (is.null(voxel)) {
    saveRDS(model_list, file = glue(filename, ".vbgamlss"))
  } else {
    warning("possibly not working as intended!")
    saveRDS(model_list[[voxel]], file = glue(filename, ".voxel{voxel}.vbgamlss"))
  }
  cat('Saved VBGAMLSS model to:', filename)
}


#' Load VBGAMLSS models from files with the specified prefix.
#'
#' @param filepath path to serialized model.
#' @return A structure containing the loaded VBGAMLSS models.
#' @export
load_model <- function(filepath) {
  # load
  return(structure(readRDS(filepath), class = "vbgamlss"))
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
#' @param subsample Optional integer vector of subject indices to predict.
#' @param ... Additional arguments passed to the predict.gamlss2 function.
#' @return A structure containing predictions.
#' @export
predict.vbgamlss <- function(object, newdata=NULL, num_cores=NULL, ptype='parameter',
                             segmentation=NULL, segmentation_target=NULL,
                             afold=NULL, subsample=NULL,
                             ...){
  if (missing(object)) { stop("vbgamlss is missing")}
  if (is.null(num_cores)) {num_cores <- availableCores()}


  # check segmentation
  if (!is.null(segmentation)){
    if (!is.data.frame(segmentation)) {stop("Error: segmentaion must be a data.frame")}
  }


  # subset the dataframe if the input is a fold from CV
  #   a fold must be a boolean vector of length of number of subjects (image 4th dim)
  if (!is.null(afold)){
    if (!is.logical(afold) && !is.integer(afold)) {
      stop("Error: afold must be either logical or integer vector.")
    }
    segmentation <- segmentation[afold,]
  }


  # subset the dataframe if a subsampling scheme is provided
  if (!is.null(subsample)){
    if (!is.numeric(subsample)) {
      stop("Error: subsample must be a numeric vector of indeces of length sum(mask>0).")
    }
    segmentation <- segmentation[,subsample]
  }


  # compute chunk size
  Nchunks <- estimate_nchunks(object)


  # predict
  plan(strategy="future::cluster", workers=num_cores)
  future.opt <- list(packages=c('gamlss2', 'gamlss'))
  fname <- as.character(object[[1]]$family)
  familyobj <- gamlss2:::complete_family(get(fname))
  predictions <- c()
  chunked = as.list(isplitVector(object, chunks=Nchunks))

  # chunks loop
  i = 1
  for (chunk in chunked){
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    i<- i+1

    # parallel call
    subpr <- pbmclapply(seq(length(chunk)),
                        function(i) {
                          vxlgamlss <- chunk[[i]]

                          # process only if properly fitted else NA
                          if ("gamlss2" %in% class(vxlgamlss)) {

                            # recon family object
                            vxlgamlss$family <- familyobj

                            # if multi tissue add
                            newdata_ <- newdata
                            if (!is.null(segmentation)){
                              newdata_$tissue <- segmentation[,i]
                              if (! is.null(segmentation_target)){
                                newdata_ <-  newdata_[newdata_$tissue == segmentation_target,]
                              }}

                            # predict
                            l <- as.list(predict(vxlgamlss,
                                                 newdata = newdata_,
                                                 type=ptype,
                                                 ...))
                            l$family <- fname
                            l$vxl <- vxlgamlss$vxl
                            l
                          } else {
                            NA
                          }
                        },
                        mc.cores=num_cores)
    predictions <- c(predictions, subpr)
  }
  gc()
  return(structure(predictions, class = "vbgamlss.predictions"))
}


#' Compute Z-scores for vbgamlss predictions given Y voxel data and image mask.
#'
#' @param predictions Predictions from vbgamlss model.
#' @param yimageframe data.frame object of the response variables to compute z-score of.
#' @param num_cores Number of CPU cores to use for parallel processing. Defaults to all available cores.
#' @return A structure containing Z-scores.
#' @export
zscore.vbgamlss <- function(predictions, yimageframe, num_cores=NULL){
  if (missing(predictions)) { stop("vbgamlss.predictions is missing")}
  if (missing(yimageframe)) { stop("yimageframe is missing")}
  if (is.null(num_cores)) {num_cores <- availableCores()}

  if (!is.data.frame(yimageframe)){stop("Error: yimageframe must be a data.frame")}

  # parallel function
  do.zscore <- function(obj) {
    pred <- obj$pred
    yval <- obj$yvxldat
    # get number of params
    lpar <- sum(names(pred) %in% c("mu", "sigma", "nu", "tau"))
    qfun <- paste("p",pred$family,sep="")
    if (lpar == 1) {
      newcall <- call(qfun, yval, mu = pred$mu)
    } else if (lpar == 2) {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma)
    } else if (lpar == 3) {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu)
    } else {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu, tau = pred$tau)
    }
    cdf <- eval(newcall)
    rqres <- qnorm(cdf)
    return(rqres)
  }

  # compute chunk size
  Nchunks <- estimate_nchunks(yimageframe)
  # parallel call
  zscores = c()
  chunked = as.list(isplitIndices(ncol(yimageframe), chunks=Nchunks))
  i = 1
  for (ichunk in chunked){
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    i<- i+1
    # chunk yimageframe and predictions
    iterable_chunk <- mapply(list,
                             pred = predictions[ichunk],
                             yvxldat = yimageframe[,ichunk],
                             SIMPLIFY=F)
    # compute z-scores
    subzs <- pbmclapply(iterable_chunk,
                        do.zscore,
                        mc.cores=num_cores)
    zscores <- c(zscores, subzs)
  }

  return(structure(zscores, class = "vbgamlss.zscores"))
}



#' Compute voxel-wise Z-scores from VBGAMLSS coefficients maps
#'
#' Given a subject image and predicted parameter maps for a GAMLSS family, compute zscore at each voxel inside a mask and return a z-score image.
#'
#' @param yimage A NIFTI image (or path readable by \code{antsImageRead}) containing the observed data to score.
#' @param mask An NIFTI image (or path) defining voxels to evaluate (nonzero entries are used).
#' @param family Character scalar naming a \pkg{gamlss.dist} family with an available CDF function \code{p<family>} (for example \code{"NO"}, \code{"SHASH"}). Default is \code{"SHASH"}.
#' @param mu_hat,sigma_hat,nu_hat,tau_hat NIFTI images (or paths) with the fitted parameter maps corresponding to the chosen family. Supply only the parameters required by the family.
#' @param num_cores Integer number of parallel workers. Defaults to \code{parallelly::availableCores()}.
#' @return An z-score image with voxel-wise z-scores.
#' @export
zscore.map.vbgamlss <- function(yimage,
                                mask,
                                family='SHASH',
                                mu_hat,
                                sigma_hat=NULL,
                                nu_hat=NULL,
                                tau_hat=NULL,
                                num_cores=NULL){
  if (is.null(num_cores)) {num_cores <- availableCores()}

  # make df
  frame = do.call(rbind,
                  list(y = images2matrix( yimage, mask),
                       mu = images2matrix( mu_hat, mask),
                       sigma = images2matrix( sigma_hat, mask),
                       nu = images2matrix( nu_hat, mask),
                       tau = images2matrix( tau_hat, mask)
                  ))
  frame <- as.data.frame(t(frame))

  # parallel function
  do.zscore <- function(i) {
    # get number of params
    lpar <- sum(names(frame) %in% c("mu", "sigma", "nu", "tau"))
    yval <- frame[i, 'y']
    pred <- frame[i, ]

    qfun <- paste("p", family, sep="")
    if (lpar == 1) {
      newcall <- call(qfun, yval, mu = pred$mu)
    } else if (lpar == 2) {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma)
    } else if (lpar == 3) {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu)
    } else {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu, tau = pred$tau)
    }
    cdf <- eval(newcall)
    rqres <- qnorm(cdf)
    return(rqres)
  }

  # compute z-scores
  subzs <- pbmclapply(1:dim(frame)[1],
                      do.zscore,
                      mc.cores=num_cores)
  zmap <- matrixToImages(matrix(unlist(subzs), nrow = 1), antsImageRead(mask,3))
  return(zmap[[1]])
}










































############################## ========== ######################################
                             # DEPRECATED #
############################## ========== ######################################

#'
#' #' @export
#' save_model_old <- function(model_list, file_prefix, num_cores = NULL) {
#'   warning("This function is deprecated. Use save_model instead.")
#'   # check and make folder
#'   parts <- unlist(strsplit(file_prefix, "/"))
#'   directory <- paste(parts[-length(parts)], collapse = "/")
#'   if (!dir.exists(directory)) {
#'     dir.create(directory, recursive = TRUE)
#'     cat("Folder created:", directory, "\n")
#'   } else {
#'     stop("Folder already exists:", directory, "\n")
#'   }
#'
#'   if (is.null(num_cores)) {num_cores <- availableCores()}
#'   # save parallel
#'   noout <- pbmclapply(1:length(model_list),
#'                       function(i) {
#'                         filename <- paste0(file_prefix, ".", model_list[[i]]$vxl)
#'                         saveRDS(model_list[[i]], file = filename)
#'                         NULL
#'                       },
#'                       mc.cores=num_cores)
#' }
#'
#'
#'
#' #' Load VBGAMLSS models from files with the specified prefix.
#' #'
#' #' @param file_prefix Prefix of the files containing the VBGAMLSS models.
#' #' @param num_cores Number of CPU cores to use for parallel processing.
#' #'   Defaults to all available cores if not provided.
#' #' @return A structure containing the loaded VBGAMLSS models.
#' #' @export
#' load_model_individually <- function(file_prefix, num_cores = NULL) {
#'   warning("This function is deprecated. Use load_model instead.")
#'   if (is.null(num_cores)) {num_cores <- availableCores()}
#'   # parse path:
#'   parts <- unlist(strsplit(file_prefix, "/"))
#'   directory <- paste(parts[-length(parts)], collapse = "/")
#'   prefix <- gsub("\\'", "", parts[length(parts)])
#'   # List files with the specified prefix
#'   file_paths <- list.files(directory, pattern = paste0("^", prefix))
#'   # Extract integer suffix from file names
#'   integer_suffixes <- sapply(file_paths, function(s) {as.numeric(gsub("^.*\\.(\\d+)$", "\\1", s))})
#'   # check if missing voxel models
#'   expected_numbers <- seq(1, length(integer_suffixes))
#'   missing_numbers <- setdiff(expected_numbers, integer_suffixes)
#'   cat(paste0('VBGAMLSS model directory:', directory))
#'   if (length(missing_numbers) == 0) {
#'     cat(paste0('Found models for ', length(integer_suffixes), ' voxels.'))
#'   } else {
#'     cat(paste0("Some models are missing. can't find model[s] for voxel: ", missing_numbers))
#'   }
#'   integer_suffixes <- as.list(sort(integer_suffixes))
#'   # run load
#'   plan(stategy="future::cluster", workers=num_cores)
#'   handlers(global = TRUE)
#'   handlers("pbmcapply")
#'   p <- progressor(along=integer_suffixes)
#'   # parallel load
#'   models <- foreach(i = integer_suffixes) %dofuture% {
#'     rds <- readRDS(paste0(file_prefix, ".", i))
#'     p()
#'     rds
#'   }
#'   gc()
#'   return(structure(models, class = "vbgamlss"))
#' }
#'
#' #' @export
#' load_model_individually_chunks <- function(file_prefix, num_cores = NULL, chunk_max_mb=256,
#'                                     index=NULL) {
#'   warning("This function is deprecated. Use load_model instead.")
#'   if (is.null(num_cores)) {num_cores <- availableCores()}
#'   # parse path:
#'   parts <- unlist(strsplit(file_prefix, "/"))
#'   directory <- paste(parts[-length(parts)], collapse = "/")
#'   prefix <- gsub("\\'", "", parts[length(parts)])
#'   # List files with the specified prefix
#'   file_paths <- list.files(directory, pattern = paste0("^", prefix))
#'   # Extract integer suffix from file names
#'   integer_suffixes <- sapply(file_paths, function(s) {as.numeric(gsub("^.*\\.(\\d+)$", "\\1", s))})
#'
#'   # check if missing voxel models
#'   expected_numbers <- seq(1, length(integer_suffixes))
#'   missing_numbers <- setdiff(expected_numbers, integer_suffixes)
#'   cat(paste0('VBGAMLSS model directory:', directory), fill=T)
#'
#'   if (length(missing_numbers) == 0) {
#'     cat(paste0('Found models for ', length(integer_suffixes), ' voxels.'), fill=T)
#'     cat(paste0('Expect long post-processing.'), fill=T)
#'   } else {
#'     warning(paste0("Some models are missing. Can't find model[s] for voxel: ", missing_numbers))
#'   }
#'   integer_suffixes <- as.list(sort(integer_suffixes))
#'
#'   # subset the voxel models to load by index
#'   if (! is.null(index)){
#'     integer_suffixes <- integer_suffixes[index]
#'     cat('Index is set, loading only the voxels models number: ',index, '\n', fill=T)
#'   }
#'   # parallel call
#'   plan(stategy="future::cluster", workers=num_cores)
#'   Nchunks = estimate_nchunks(paste0(file_prefix, ".", 1),
#'                              from_files=TRUE,
#'                              chunk_max_Mb=chunk_max_mb)
#'   chunked = as.list(isplitIndices(length(integer_suffixes), chunks=Nchunks))
#'   models = c()
#'   i = 1
#'   for (ichunk in chunked){
#'     cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
#'     i<- i+1
#'     integer_suffixes_chunked <- integer_suffixes[ichunk]
#'     p <- progressor(length(integer_suffixes_chunked))
#'     # parallel loop
#'     loaded <- foreach(i = integer_suffixes_chunked) %dofuture% {
#'       rds <- readRDS(paste0(file_prefix, ".", i))
#'       p()
#'       rds
#'     }
#'     models <- c(models, loaded)
#'   }
#'   gc()
#'   return(structure(models, class = "vbgamlss"))
#' }
