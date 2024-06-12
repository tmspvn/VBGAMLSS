#' Run a GAMLSS2 on all voxels of a NIfTI image within a mask. 
# TODOS:
#       - add batching on the rest of the functions, especially the fitting, 
#         and loading
#       - convert inputting from 4D image to list of 3rd images to be open by 
#         antsR imagesToMatrix(imageList, mask)

library(future)
library(future.apply)
library(gamlss)
library(gamlss2) # devtools::install_github("gamlss-dev/gamlss2")
library(pbmcapply)
library(voxel)
library(doFuture)
library(itertools)
library(ANTsR)
library(progressr)
options(progressr.enable=TRUE)


################################  Core  ########################################

vbgamlss <- function(image, mask, g.formula, train.data, g.family=NO, 
                     segmentation=NULL,
                     num_cores=NULL, ...) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  if (missing(train.data)) { stop("subjData is missing")}
  
  if (class(g.formula) != "character") { stop("g.formula class must be character")}
  if (is.null(num_cores)) {num_cores <- availableCores()}
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  voxeldata <- images2matrix(image, mask)
  if (!is.null(segmentation)){segmentation <- images2matrix(segmentation, mask)}
  gc()

  # Coerce update on left hand g.formula for parallel fitting
  g.fo <- update.formula(g.formula, 'Y ~ .')
  
  # parallel settings
  plan(cluster, workers = num_cores)
  handlers(global = TRUE)
  handlers("pbmcapply")
  p <- with_progress(progressor(ncol(voxeldata)))
  future.opt <- list(packages=c('gamlss2'))
  # foreach call
  models <- foreach(vxlcol = seq_along(voxeldata),
                    .options.future = future.opt,
                    .combine=c
                    ) %dofuture% {
                               # fit specific voxel
                               vxl_train_data <- train.data
                               vxl_train_data$Y <- voxeldata[,vxlcol]
                               # if multi tissue add
                               if (!is.null(segmentation)){
                                 vxl_train_data$tissue <- segmentation[,vxlcol]}
                               # gamlss
                               g <- gamlss2::gamlss2(formula=g.fo,
                                                     data=vxl_train_data,
                                                     family=g.family,
                                                     light=TRUE, 
                                                     trace=FALSE,
                                                     ...)
                               g$control <- NULL
                               g$family <- g$family[1] # reset via gamlss2:::complete_family()
                               g$vxl <- vxlcol
                               p() # update progressbar
                               list(g)
                    }
  gc()
  return(structure(models, class = "vbgamlss"))
}


vbgamlss_chunks <- function(image, mask, g.formula, train.data, g.family=NO,
                            segmentation=NULL, 
                            num_cores=NULL, 
                            chunk_max_mb=64,
                            ...) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  if (missing(train.data)) { stop("subjData is missing")}
  
  if (class(g.formula) != "character") { stop("g.formula class must be character")}
  if (is.null(num_cores)) {num_cores <- availableCores()}
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  voxeldata <- images2matrix(image, mask)
  if (!is.null(segmentation)){segmentation <- images2matrix(segmentation, mask)}
  gc()
  
  # Coerce update on left hand g.formula for parallel fitting
  g.fo <- update.formula(g.formula, 'Y ~ .')
  
  # parallel settings
  plan(cluster, workers = num_cores, rscript_libs = .libPaths())
  handlers(global = TRUE)
  handlers("pbmcapply")
  future.opt <- list(packages=c('gamlss2'))
  # compute chunk size
  Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
  # parallel call
  chunked = as.list(isplitIndices(ncol(voxeldata), chunks=Nchunks))
  i = 1
  models <- list()
  for (ichunk in chunked){
      cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
      i<- i+1
      voxeldata_chunked <- voxeldata[,ichunk]
      p <- progressor(ncol(voxeldata_chunked))
      # fit vbgamlss 
      submodels <- foreach(vxlcol = seq_along(voxeldata_chunked),
                        .options.future = future.opt,
                        .combine=c
      ) %dofuture% {
        # fit specific voxel
        vxl_train_data <- train.data
        vxl_train_data$Y <- voxeldata_chunked[,vxlcol]
        # if multi tissue add
        if (!is.null(segmentation)){
          vxl_train_data$tissue <- segmentation[,vxlcol]}
        # gamlss
        g <- gamlss2::gamlss2(formula=g.fo,
                              data=vxl_train_data,
                              family=g.family,
                              light=TRUE, 
                              trace=FALSE,
                              ...)
        g$control <- NULL
        g$family <- g$family[1] # re-set via gamlss2:::complete_family()
        g$vxl <- vxlcol
        p() # update progressbar
        list(g)
    }
    models <- c(models, submodels)
  }
  gc()
  return(structure(models, class = "vbgamlss"))
}


save_model <- function(model_list, file_prefix, num_cores = NULL) {
  # check and make folder
  parts <- unlist(strsplit(file_prefix, "/"))
  directory <- paste(parts[-length(parts)], collapse = "/")
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    cat("Folder created:", directory, "\n")
  } else {
    stop("Folder already exists:", directory, "\n")
  }

  if (is.null(num_cores)) {num_cores <- availableCores()}
  # save parallel
  noout <- pbmclapply(1:length(model_list), 
           function(i) {
              filename <- paste0(file_prefix, ".", model_list[[i]]$vxl)
              saveRDS(model_list[[i]], file = filename)
              NULL
           },
           mc.cores=num_cores)
  }


#' Load VBGAMLSS models from files with the specified prefix.
#'
#' @param file_prefix Prefix of the files containing the VBGAMLSS models.
#' @param num_cores Number of CPU cores to use for parallel processing.
#'   Defaults to all available cores if not provided.
#' @return A structure containing the loaded VBGAMLSS models.
#' @importFrom future plan handlers progressor
#' @importFrom parallel availableCores cluster
#' @import foreach %dofuture% readRDS
#' @export
load_model <- function(file_prefix, num_cores = NULL) {
  if (is.null(num_cores)) {num_cores <- availableCores()}
  # parse path:
  parts <- unlist(strsplit(file_prefix, "/"))
  directory <- paste(parts[-length(parts)], collapse = "/")
  prefix <- gsub("\\'", "", parts[length(parts)])
  # List files with the specified prefix
  file_paths <- list.files(directory, pattern = paste0("^", prefix))
  # Extract integer suffix from file names
  integer_suffixes <- sapply(file_paths, function(s) {as.numeric(gsub("^.*\\.(\\d+)$", "\\1", s))})
  # check if missing voxel models
  expected_numbers <- seq(1, length(integer_suffixes))
  missing_numbers <- setdiff(expected_numbers, integer_suffixes)
  cat(paste0('VBGAMLSS model directory:', directory))
  if (length(missing_numbers) == 0) {
    cat(paste0('Found models for ', length(integer_suffixes), ' voxels.'))
  } else {
    cat(paste0("Some models are missing. can't find model[s] for voxel: ", missing_numbers))
  }
  integer_suffixes <- as.list(sort(integer_suffixes))
  # run load
  plan(cluster, workers = num_cores)
  handlers(global = TRUE)
  handlers("pbmcapply")
  p <- progressor(along=integer_suffixes)
  # parallel load
  models <- foreach(i = integer_suffixes) %dofuture% {
    rds <- readRDS(paste0(file_prefix, ".", i))
    p()
    rds
  }
  gc()
  return(structure(models, class = "vbgamlss"))
}


load_model_chunks <- function(file_prefix, num_cores = NULL) {
  if (is.null(num_cores)) {num_cores <- availableCores()}
  # parse path:
  parts <- unlist(strsplit(file_prefix, "/"))
  directory <- paste(parts[-length(parts)], collapse = "/")
  prefix <- gsub("\\'", "", parts[length(parts)])
  # List files with the specified prefix
  file_paths <- list.files(directory, pattern = paste0("^", prefix))
  # Extract integer suffix from file names
  integer_suffixes <- sapply(file_paths, function(s) {as.numeric(gsub("^.*\\.(\\d+)$", "\\1", s))})
  
  # check if missing voxel models
  expected_numbers <- seq(1, length(integer_suffixes))
  missing_numbers <- setdiff(expected_numbers, integer_suffixes)
  cat(paste0('VBGAMLSS model directory:', directory), fill=T)
  
  if (length(missing_numbers) == 0) {
    cat(paste0('Found models for ', length(integer_suffixes), ' voxels.'), fill=T)
  } else {
    cat(paste0("Some models are missing. Can't find model[s] for voxel: ", missing_numbers), fill=T)
  }
  integer_suffixes <- as.list(sort(integer_suffixes))
  # parallel call
  Nchunks = estimate_nchunks(paste0(file_prefix, ".", 1), from_files=TRUE)
  chunked = as.list(isplitIndices(length(integer_suffixes), chunks=Nchunks))
  models = c()
  i = 1
  for (ichunk in chunked){
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    i<- i+1
    # compute z-scores
    integer_suffixes_chunked <- integer_suffixes[ichunk]
    p <- progressor(length(integer_suffixes_chunked))
    loaded <- foreach(i = integer_suffixes_chunked) %dofuture% {
      rds <- readRDS(paste0(file_prefix, ".", i))
      p()
      rds
    }
    models <- c(models, loaded)
  }
  gc()
  return(structure(models, class = "vbgamlss"))
}


#' Predict method for vbgamlss objects
#'
#' @param object vbgamlss object to make predictions from.
#' @param newdata New data to use for predictions.
#' @param num_cores Number of CPU cores to use for parallel processing.
#'   Defaults to one less than the total available cores if not provided.
#' @param ... Additional arguments to be passed to the \code{\link{predict}} function.
#' @return A structure containing predictions.
#' @export
predict.vbgamlss <- function(object, newdata=NULL, num_cores=NULL, ...){
  if (missing(object)) { stop("vbgamlss is missing")}
  if (is.null(num_cores)) {num_cores <- availableCores()-1}
  # compute chunk size
  Nchunks <- estimate_nchunks(object)
  # predict
  plan(cluster, workers = num_cores)
  future.opt <- list(packages=c('gamlss2', 'gamlss'))
  fname <- as.character(object[[1]]$family)
  familyobj <- gamlss2:::complete_family(get(fname))
  predictions <- c()
  chunked = as.list(isplitVector(object, chunks=Nchunks))
  i = 1
  for (chunk in chunked){
      cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
      i<- i+1
      subpr <- pbmclapply(chunk, 
                         function(vxlgamlss) {
                          vxlgamlss$family <- familyobj
                          l <- as.list(predict(vxlgamlss, 
                                               newdata = newdata,
                                               ...))
                          l$family <- fname
                          l$vxl <- vxlgamlss$vxl
                          l
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
#' @param yimage R object or filename containing image data.
#' @param ymask R object or filename containing mask data.
#' @param num_cores Number of CPU cores to use for parallel processing. Defaults to all available cores.
#' @return A structure containing Z-scores.
#' @import oro.nifti
#' @importFrom parallel availableCores
#' @importFrom parallel mclapply
#' @export
zscore.vbgamlss <- function(predictions, yimage, ymask, num_cores=NULL){
  if (missing(predictions)) { stop("vbgamlss.predictions is missing")}
  if (missing(yimage)) { stop("vbgamlss.predictions is missing")}
  if (is.null(num_cores)) {num_cores <- availableCores()}
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  yvoxeldata <- images2matrix(yimage, ymask)
  gc()
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
  Nchunks <- estimate_nchunks(yvoxeldata)
  # parallel call
  zscores = c()
  chunked = as.list(isplitIndices(ncol(yvoxeldata), chunks=Nchunks))
  i = 1
  for (ichunk in chunked){
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    i<- i+1
    # chunk yvoxeldata and predictions
    iterable_chunk <- mapply(list, 
                             pred = predictions[ichunk],
                             yvxldat = yvoxeldata[,ichunk], 
                             SIMPLIFY=F)
    # compute z-scores
    subzs <- pbmclapply(iterable_chunk, 
                        do.zscore, 
                        mc.cores=num_cores)
    zscores <- c(zscores, subzs)
  }
  
  return(structure(zscores, class = "vbgamlss.zscores"))
}



################################  Utilities  ###################################

#' Estimate number of chunks based on object size
#'
#' @param object R object for chunk size estimation.
#' @return Estimated number of chunks with max 256MB per job.
#' @examples estimate_nchunks(data)
#' @export
estimate_nchunks <- function(object, from_files=F, chunk_max_Mb=256) {
  # compute chunk size of max 256 mb per job
  if (! from_files){
      memory_size_mb <- object.size(object) / (1024^2)
      Nchunks <- ceiling(memory_size_mb / chunk_max_Mb) 
    } else {
      memory_size_mb <- file.info(object)$size / (1024^2)
      Nchunks <- ceiling(memory_size_mb / chunk_max_Mb) 
    }
  return(Nchunks)
  }


images2matrix <- function(image_list, mask) {
  # imageList is a list containing images.  Mask is a mask image Returns matrix of
  # dimension (numImages, numVoxelsInMask)
  if (missing(image_list)) { stop("image_list is missing")}
  if (missing(mask)) { stop("mask is missing")}
  # load mask
  mask_ <- antsImageRead(mask, 3)
  mask_arr <- as.array(mask_) > 0
  numVoxels <- length(which(mask_arr))
  # check if list of images of image
  if (class(image_list) == "list"){
    cat('Converting from list of images paths', fill=T)
    numImages <- length(image_list)
    dataMatrix <- matrix(nrow = numImages, ncol = numVoxels)
    for (i in 1:length(imageList)) 
      {dataMatrix[i, ] <- as.numeric(imageList[[i]], mask = mask_arr)}
    
  } else {
    if (class(image_list) != "character") { 
      stop("image_list is not a list of nifties or a path to a 4D image")}
    cat('Converting from 4D image path', fill=T)
    large_image <- antsImageRead(image_list, 4)
    numImages <- dim(large_image)[4]
    dataMatrix <- matrix(nrow = numImages, ncol = numVoxels)
    for (i in 1:numImages) 
      {dataMatrix[i, ] <- as.numeric(as.antsImage(drop(large_image[,,,i])), mask = mask_arr)}
    
  }
  return(as.data.frame(dataMatrix))
}


matrix2image <- function(matrix, mask){
  if (missing(matix)) { stop("matrix is missing")}
  if (missing(mask)) { stop("mask is missing")}
  return()
}


map_model_coefficients <- function(fittedobj, mask, filename, return_files=FALSE){
  if (class(fittedobj) != "vbgamlss") { stop("fittedobj must be of class vbgamlss.")}
  message('Warning, specific coefficients from special model terms cannot be map (e.g. pb(), s())')
  nvox <- length(fittedobj)
  first_mod_coefs <- unlist(fittedobj[[1]]$coefficients)
  name_coefs <- names(first_mod_coefs)
  ncoefs <- length(first_mod_coefs)
  coefs_mat <- matrix(nrow = ncoefs, ncol = nvox)
  for (i in 1:nvox) {
    coefs_mat[,i] <- unlist(fittedobj[[i]]$coefficients)
  }
  # convert mat to maps
  coef_maps_images <- matrixToImages(coefs_mat, antsImageRead(mask, 3))
  # save files
  fnames <- c()
  for (i in 1:length(coef_maps_images)) {
    # parse coef name
    parts <- strsplit(name_coefs[i], "\\.")[[1]]
    if (parts[2] != '(Intercept)') {parts[2] = paste0('(', parts[2], ')')}
    fname <- paste0(filename, '_par-', toupper(parts[1]), '_coef-', parts[2],'.nii.gz')
    # save
    antsImageWrite(coef_maps_images[[i]], fname)
    fnames[i] <- fname
  }
  if (return_files) {return(fnames)}
}


map_model_predictions <- function(obj, mask, filename, index=NULL, 
                                  return_files=FALSE){
  if (class(obj) != "vbgamlss.predictions") { stop("obj must be of class vbgamlss.predictions")}
  if (! is.null(index)) {cat(paste0('Mapping predictions index: ', list(index)), fill=T)}
  # prepare usefull info
  nvox <- length(obj)
  first_pred <- obj[[1]]
  family <- first_pred$family
  name_param <- names(first_pred)[! names(first_pred) %in% c('family', 'vxl')]
  nparam <- length(name_param)
  # save a subset?
  if (is.null(index)) {
    nsubj <- length(first_pred[[name_param[1]]])
    subj = 1:nsubj # all
  } else {
    subj = index # subset
    nsubj = length(subj)
  }
  # process
  for (pname in name_param) {
    param_mat <- matrix(nrow=nsubj, ncol=nvox)
    for (ic in 1:nvox) {
      # voxel   #parameter  #subjects
      param_mat[, ic] <- obj[[ic]][[pname]][subj]
    }
    # convert mat to maps & save per subj
    param_maps_images <- matrixToImages(param_mat, antsImageRead(mask, 3))
    # save files
    fnames <- c()
    for (ip in 1:length(param_maps_images)) {
      fname <- paste0(filename, 
                      '_subj-', index[ip],
                      '_fam-', family,
                      '_par-' , toupper(pname),
                      '.nii.gz')
      antsImageWrite(param_maps_images[[ip]], fname)
      fnames[ip] <- fname
    }
  }
  if (return_files) {return(fnames)}
}


map_zscores <- function(zscores, mask, filename, index=NULL, 
                                  return_files=FALSE){
  if (class(zscores) != "vbgamlss.zscores") { stop("zscores must be of class vbgamlss.zscores")}
  if (! is.null(index)) {cat(paste0('Mapping predictions index: ', list(index)), fill=T)}
  
  # prepare usefull info
  nvox <- length(zscores)
  first_vxl <- zscores[[1]]
  # save a subset?
  if (is.null(index)) {
    nsubj <- length(first_vxl)
    subj = 1:nsubj # all
  } else {
    subj = index # subset
    nsubj = length(subj)
  }
  # process
  z_mat <- matrix(nrow=nsubj, ncol=nvox)
  for (ic in 1:nvox) {
    # voxel #subjects
    z_mat[, ic] <- zscores[[ic]][subj]
  }
  # convert mat to maps & save per subj
  z_maps_images <- matrixToImages(z_mat, antsImageRead(mask, 3))
  # save files
  fnames <- c()
  for (ip in 1:length(z_maps_images)) {
    fname <- paste0(filename, 
                    '_subj-', index[ip],
                    '.zscore.nii.gz')
    antsImageWrite(z_maps_images[[ip]], fname)
    fnames[ip] <- fname
  }
  if (return_files) {return(fnames)}
}


################################  Development ##################################
zscore_image.vbgamlss <- function(predictions, yimage, ymask, num_cores=NULL){
  if (missing(predictions)) { stop("vbgamlss.predictions is missing")}
  if (missing(yimage)) { stop("vbgamlss.predictions is missing")}
  if (is.null(num_cores)) {num_cores <- availableCores()}
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  yvoxeldata <- images2matrix(yimage, ymask)
  gc()
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
  Nchunks <- estimate_nchunks(yvoxeldata)
  # parallel call
  zscores = c()
  chunked = as.list(isplitIndices(ncol(yvoxeldata), chunks=Nchunks))
  i = 1
  for (ichunk in chunked){
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    i<- i+1
    # chunk yvoxeldata and predictions
    iterable_chunk <- mapply(list, 
                             pred = predictions[ichunk],
                             yvxldat = yvoxeldata[,ichunk], 
                             SIMPLIFY=F)
    # compute z-scores
    subzs <- pbmclapply(iterable_chunk, 
                        do.zscore, 
                        mc.cores=num_cores)
    zscores <- c(zscores, subzs)
  }
  
  return(structure(zscores, class = "vbgamlss.zscores"))
}
