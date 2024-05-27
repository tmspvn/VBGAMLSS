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


vbgamlss <- function(image, mask, multi.tissue=NULL, 
                     g.formula, train.data, g.family=NO,
                     num_cores=NULL, ...) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  if (missing(train.data)) { stop("subjData is missing")}
  
  if (class(g.formula) != "character") { stop("g.formula class must be character")}
  
  if (class(image) == "character") {image <- oro.nifti::readNIfTI(fname=image)} 
  
  if (class(mask) == "character") {mask <- oro.nifti::readNIfTI(fname=mask)}
  
  if (class(multi.tissue) == "character") {multi.tissue <- oro.nifti::readNIfTI(fname=multi.tissue)} 
  
  if (is.null(num_cores)) {num_cores <- availableCores()}
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  voxeldata <- ts2matrix(image, mask)
  if (!is.null(multi.tissue)){multi.tissue <- ts2matrix(multi.tissue, mask)}
  rm(image)
  rm(mask)
  gc()

  # Coerce update on left hand g.formula for parallel fitting
  g.fo <- update.formula(g.formula, 'Y ~ .')
  
  # parallel settings
  plan(cluster, workers = num_cores)
  handlers(global = TRUE)
  handlers("pbmcapply")
  p <- progressor(ncol(voxeldata))
  future.opt <- list('gamlss2')
  # foreach call --> MAY NEED CHUNKING WITH MORE VOXEL DATA 
  models <- foreach(vxlcol = seq_along(voxeldata),
                    .options.future = future.opt,
                    .combine=c
                    ) %dofuture% {
                               # fit specific voxel
                               vxl_train_data <- train.data
                               vxl_train_data$Y <- voxeldata[,vxlcol]
                               # if multi tissue add
                               if (!is.null(multi.tissue)){
                                 vxl_train_data$tissue <- multi.tissue[,vxlcol]}
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



vbgamlss_chunks <- function(image, mask, multi.tissue=NULL, 
                     g.formula, train.data, g.family=NO,
                     num_cores=NULL, ...) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  if (missing(train.data)) { stop("subjData is missing")}
  
  if (class(g.formula) != "character") { stop("g.formula class must be character")}
  
  if (class(image) == "character") {image <- oro.nifti::readNIfTI(fname=image)} 
  
  if (class(mask) == "character") {mask <- oro.nifti::readNIfTI(fname=mask)}
  
  if (class(multi.tissue) == "character") {multi.tissue <- oro.nifti::readNIfTI(fname=multi.tissue)} 
  
  if (is.null(num_cores)) {num_cores <- availableCores()}
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  voxeldata <- ts2matrix(image, mask)
  if (!is.null(multi.tissue)){multi.tissue <- ts2matrix(multi.tissue, mask)}
  rm(image)
  rm(mask)
  gc()
  
  # Coerce update on left hand g.formula for parallel fitting
  g.fo <- update.formula(g.formula, 'Y ~ .')
  
  # parallel settings
  plan(cluster, workers = num_cores)
  handlers(global = TRUE)
  handlers("pbmcapply")
  future.opt <- list('gamlss2')
  # compute chunk size
  Nchunks <- estimate_nchunks(voxeldata)
  # parallel call
  chunked = as.list(isplitIndices(ncol(voxeldata), chunks=Nchunks))
  i = 1
  models <- list()
  for (ichunk in chunked){
      print(paste0("Chunk: ",i,"/", Nchunks))
      i<- i+1
      voxeldata_chunked <- voxeldata[,ichunk]
      p <- progressor(ncol(voxeldata_chunked))
      # fit vbgamlss 
      models <- foreach(vxlcol = seq_along(voxeldata_chunked),
                        .options.future = future.opt,
                        .combine=c
      ) %dofuture% {
        # fit specific voxel
        vxl_train_data <- train.data
        vxl_train_data$Y <- voxeldata_chunked[,vxlcol]
        # if multi tissue add
        if (!is.null(multi.tissue)){
          vxl_train_data$tissue <- multi.tissue[,vxlcol]}
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
  print(paste0('VBGAMLSS model directory:', directory))
  if (length(missing_numbers) == 0) {
    print(paste0('Found models for ', length(integer_suffixes), ' voxels.'))
  } else {
    print(paste0("Some models are missing. can't find model[s] for voxel: ", missing_numbers))
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
      print(paste0("Chunk: ",i,"/", Nchunks))
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
  
  if (class(yimage) == "character") {yimage <- oro.nifti::readNIfTI(fname=yimage)} 
  
  if (class(ymask) == "character") {ymask <- oro.nifti::readNIfTI(fname=ymask)}

  if (is.null(num_cores)) {num_cores <- availableCores()}
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  yvoxeldata <- ts2matrix(yimage, ymask)
  rm(yimage)
  rm(ymask)
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
    print(paste0("Chunk: ",i,"/", Nchunks))
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


#' Estimate number of chunks based on object size
#'
#' @param object R object for chunk size estimation.
#' @return Estimated number of chunks with max 256MB per job.
#' @examples estimate_nchunks(data)
#' @export
estimate_nchunks <- function(object) {
                          # compute chunk size of max 256 mb per job
                          memory_size_mb <- object.size(object) / (1024^2)
                          Nchunks <- ceiling(memory_size_mb / 256) 
                          return(Nchunks)
                          }

# to be finished:
load_input <- function(image_list, mask){
  if (missing(image_list)) { stop("image_list is missing")}
  if (missing(mask)) { stop("mask is missing")}
  
  if (class(image_list) == "list") {yimage <- oro.nifti::readNIfTI(fname=yimage)} 
  
  if (class(mask) == "character") {ymask <- oro.nifti::readNIfTI(fname=ymask)}
  
  if (is.null(num_cores)) {num_cores <- availableCores()}
  
}













