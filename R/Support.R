################################  Support  #####################################








#' @export
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
#' @import foreach
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

#' @export
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
