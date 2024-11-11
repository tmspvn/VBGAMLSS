################################  Core  ########################################
# install.packages("gamlss2",
#                  repos = c("https://gamlss-dev.R-universe.dev",
#                   "https://cloud.R-project.org"))
# devtools::install_github('ANTsX/ANTsR')
#
# To test if valid response
# fo <- gamlss2:::complete_family('BCPEo')
# fo$valid.response()
#
# TODOs? should save enviroment be added in while logging? even with massive sizes objects?
#
#
#
#' @title Fit Generalized Additive Models for Location Scale and Shape (GAMLSS) with Voxel Data
#'
#' @description
#' The `vbgamlss` function fits Bayesian GAMLSS models on voxel data with optional segmentation,
#' allowing for parallel processing. This function can be used to analyze high-dimensional imaging data.
#'
#' @param imageframe A data frame containing the voxel data. Each column represents a voxel,
#'                   and each row represents a subject.
#' @param g.formula A formula for the GAMLSS model.
#' @param train.data A data frame containing subject-level data. Columns must correspond to covariates in `g.formula`.
#' @param g.family A GAMLSS family object. Defaults to `NO` (Normal).
#' @param segmentation Optional data frame containing segmentation information.
#' @param segmentation_target Optional target for the segmentation data to subset the analysis.
#' @param num_cores Number of cores for parallel processing. If NULL, all available cores are used.
#' @param chunk_max_mb Maximum chunk size in megabytes for processing. Defaults to 64 MB. Increase to 256/512 for HPC.
#' @param afold Optional boolean or integer vector to subset imageframe for cross-validation folds.
#' @param subsample Optional numeric vector specifying indices of columns in imageframe to subset for analysis.
#' @param debug Logical. If TRUE, enables debug mode to log output in `logdir`. Defaults to FALSE.
#' @param logdir Directory path for saving logs if `debug` is TRUE. Defaults to current working directory.
#' @param ... Additional arguments passed to the GAMLSS fitting function.
#'
#' @return A list of GAMLSS models, with class `"vbgamlss"`.
#'
#' @details
#' The function performs the following steps:
#' - Checks input data and prepares voxel and segmentation data as needed.
#' - If parallel processing is enabled, splits voxel data into chunks based on memory constraints.
#' - Fits a GAMLSS model to each voxel using specified covariates and segmentation data, if available.
#' - Returns a list of models, one for each voxel.
#'
#' @export
vbgamlss <- function(imageframe,
                     g.formula,
                     train.data,
                     g.family=NO,
                     segmentation=NULL,
                     segmentation_target=NULL,
                     num_cores=NULL,
                     chunk_max_mb=64,
                     afold=NULL,
                     subsample=NULL,
                     debug=F, # toggle debugging
                     logdir=getwd(), # debug directory
                     chace=F, # save temporary states and force debug=T
                     force_ypositivity=T, # force Y >0
                     ...) {

  # checks
  if (missing(imageframe)) { stop("imageframe is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}


  # Force character columns to factors
  train.data <- as.data.frame(unclass(train.data), stringsAsFactors=TRUE)


  # Cores
  if (is.null(num_cores)) {num_cores <- availableCores()}


  # Input image: subjects (4th dim) x voxels (columns)
  if (!is.data.frame(imageframe)) {stop("Error: imageframe must be a data.frame")}
  voxeldata <- imageframe


  # Segmentation if provided
  if (!is.null(segmentation)){
    if (!is.data.frame(segmentation)) {stop("Error: segmentaion must be a data.frame")}
    }
  gc()


  # subset the imageframe if the input is a fold from CV
  #   a fold must be a boolean vector of length of number of subjects (image 4th dim)
  if (!is.null(afold)){
    if (!is.logical(afold) && !is.integer(afold)) {
      stop("Error: afold must be either logical or integer vector.")
    }
    voxeldata <- voxeldata[afold,]
    segmentation <- segmentation[afold,]
    train.data <- train.data[afold,]
  }


  # subset the imageframe if a subsampling scheme is provided
  if (!is.null(subsample)){
    if (!is.numeric(subsample)) {
      stop("Error: subsample must be a numeric vector of indeces of length sum(mask>0).")
    }
    voxeldata <- voxeldata[,subsample]
    segmentation <- segmentation[,subsample]
  }


  # parallel settings
  plan(stategy="future::cluster", workers=num_cores, rscript_libs=.libPaths())
  options(future.globals.maxSize=20000*1024^2)
  handlers(global = TRUE)
  handlers("pbmcapply")
  p <- with_progress(progressor(ncol(voxeldata)))
  future.opt <- list(packages=c('gamlss2'), seed = TRUE)


  # compute chunk size
  Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
  chunked = as.list(isplitIndices(ncol(voxeldata), chunks=Nchunks))


  # logdir
  if (cache) {debug = T} # force debugging while caching
  if (debug) {
    FID = rand_names(1, l=3)
    # path
    logdir=file.path(logdir, paste0('.voxlog.', FID))
    # create dir
    if (dir.exists(logdir)){
      cat("Existing log dir found, saving and resuming processed chunks if found", fill=T)
    } else {
      cat("No log dir found, creating a new one", fill=T)
      dir.create(logdir, recursive = T, showWarnings = F)
    }
  }


  # loop chunks call
  i = 1
  models <- list()
  for (ichunk in chunked){

    # Counters, names, etc..
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    chunk_id = file.path(logdir, paste0('.voxchunk.', i,))
    i <- i+1


    # subset with chunk indexes
    voxeldata_chunked <- voxeldata[,ichunk]
    if (!is.null(segmentation)){voxelseg_chunked <- segmentation[,ichunk]}


    # set up progress bar
    p <- progressor(ncol(voxeldata_chunked))


    # Continue if submodels already computed
    if (cache){
      if (file.exists(chunk_id)){
        cat("Chunk already processed, skipping", fill=T)
        submodels <- readRDS(chunk_id)
        models <- c(models, submodels)
        next
      }
    }


    # parallel call to fit vbgamlss
    submodels <- foreach(vxlcol = seq_along(voxeldata_chunked),
                         .options.future = future.opt,
                         .combine=c
    ) %dofuture% {
      # fit specific voxel
      vxl_train_data <- train.data
      vxl_train_data$Y <- as.numeric(voxeldata_chunked[,vxlcol])
      # if multi tissue add
      if (!is.null(segmentation)){
        vxl_train_data$tissue <- voxelseg_chunked[,vxlcol]}
      if (! is.null(segmentation_target)){
        vxl_train_data <- vxl_train_data[vxl_train_data$tissue == segmentation_target,]
      }


      # debug
      if (debug) {logfile=file.path(logdir, paste0('log.vxl', vxlcol))} else {logfile=NULL}

      # force Y to be strictly non-negative & !=0
      if (force_ypositivity){
        non_negative = vxl_train_data$Y>0
        vxl_train_data <- vxl_train_data[non_negative,]
      }

      # gamlss
      g <- TRY(gamlss2::gamlss2(formula=as.formula(g.formula),
                                data=vxl_train_data,
                                family=g.family,
                                light=TRUE,
                                trace=FALSE,
                                CG = 100,
                                maxit = c(300, 30),
                                ...),
               logfile, save.env.and.stop=F)
      g$control <- NULL
      g$family <- g$family$family # to re-set via gamlss2:::complete_family()
      g$vxl <- vxlcol
      p() # update progress bar
      list(g)
    }


    # Save chunk
    if (cache){saveRDS(submodels, chunk_id)}


    # Append results to models
    models <- c(models, submodels)

  }
  gc()
  return(structure(models, class = "vbgamlss"))
}





### vbgamlss with chunking ###
#' @export
vbgamlss_chunks <- function(image, mask, g.formula, train.data, g.family=NO,
                            segmentation=NULL,
                            segmentation_target=NULL,
                            num_cores=NULL,
                            chunk_max_mb=64,
                            afold=NULL,
                            subsample=NULL,
                            debug=F,
                            logdir=getwd(),
                            ...) {

  # checks
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}


  # Force character columns to factors
  train.data <- as.data.frame(unclass(train.data), stringsAsFactors=TRUE)


  # Cores
  if (is.null(num_cores)) {num_cores <- availableCores()}


  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  voxeldata <- load_input_image(image, mask)
  # load segmentation if provided
  if (!is.null(segmentation)){segmentation <- images2matrix(segmentation, mask)}
  gc()



  # subset the dataframe if the input is a fold from CV
  #   a fold must be a boolean vector of length of number of subjects (image 4th dim)
  if (!is.null(afold)){
    if (!is.logical(afold) && !is.integer(afold)) {
      stop("Error: afold must be either logical or integer vector.")
    }
    voxeldata <- voxeldata[afold,]
    segmentation <- segmentation[afold,]
    train.data <- train.data[afold,]
  }


  # subset the dataframe if a subsampling scheme is provided
  if (!is.null(subsample)){
    if (!is.numeric(subsample)) {
      stop("Error: subsample must be a numeric vector of indeces of length sum(mask>0).")
    }
    voxeldata <- voxeldata[,subsample]
    segmentation <- segmentation[,subsample]
  }


  # parallel settings
  plan(stategy="future::cluster", workers=num_cores, rscript_libs=.libPaths())
  options(future.globals.maxSize=20000*1024^2)
  handlers(global = TRUE)
  handlers("pbmcapply")
  p <- with_progress(progressor(ncol(voxeldata)))
  future.opt <- list(packages=c('gamlss2'), seed = TRUE)


  # compute chunk size
  Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
  chunked = as.list(isplitIndices(ncol(voxeldata), chunks=Nchunks))


  # logdir
  if (debug) {
    logdir=file.path(logdir, paste0('.voxlog.', rand_names(1, l=3)))
    dir.create(logdir, recursive = T, showWarnings = F)
    }


  # loop chunks call
  i = 1
  models <- list()
  for (ichunk in chunked){
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    i<- i+1
    voxeldata_chunked <- voxeldata[,ichunk]
    if (!is.null(segmentation)){voxelseg_chunked <- segmentation[,ichunk]}
    p <- progressor(ncol(voxeldata_chunked))


    # parallel call to fit vbgamlss
    submodels <- foreach(vxlcol = seq_along(voxeldata_chunked),
                         .options.future = future.opt,
                         .combine=c
    ) %dofuture% {
      # fit specific voxel
      vxl_train_data <- train.data
      vxl_train_data$Y <- as.numeric(voxeldata_chunked[,vxlcol])
      # if multi tissue add
      if (!is.null(segmentation)){
        vxl_train_data$tissue <- voxelseg_chunked[,vxlcol]}
        if (! is.null(segmentation_target)){
          vxl_train_data <- vxl_train_data[vxl_train_data$tissue == segmentation_target,]
        }


      # debug
      if (debug) {logfile=file.path(logdir, paste0('log.vxl', vxlcol))} else {logfile=NULL}
      # Y must be strictly non-negative & !=0
      non_negative = vxl_train_data$Y>0
      vxl_train_data <- vxl_train_data[non_negative,]
      # gamlss
      g <- TRY(gamlss2::gamlss2(formula=as.formula(g.formula),
                         data=vxl_train_data,
                         family=g.family,
                         light=TRUE,
                         trace=FALSE,
                         CG = 100,
                         maxit = c(300, 30),
                         ...),
            logfile, save.env.and.stop=F)
      g$control <- NULL
      g$family <- g$family$family # re-set via gamlss2:::complete_family()
      g$vxl <- vxlcol
      p() # update progressbar
      list(g)
    }

    models <- c(models, submodels)

  }
  gc()
  return(structure(models, class = "vbgamlss"))
}






# ========== #
# DEPRECATED #
# ========== #

### vbgamlss all loaded in mem ###
#
#
# vbgamlss_unchunked <- function(image, mask, g.formula, train.data, g.family=NO,
#                      segmentation=NULL,
#                      num_cores=NULL, ...) {
#
#   if (missing(image)) { stop("image is missing")}
#   if (missing(mask)) { stop("mask is missing")}
#   if (missing(g.formula)) { stop("formula is missing")}
#   if (missing(train.data)) { stop("subjData is missing")}
#
#   if (class(g.formula) != "character") { stop("g.formula class must be character")}
#   if (is.null(num_cores)) {num_cores <- availableCores()}
#
#   # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
#   voxeldata <- images2matrix(image, mask)
#   if (!is.null(segmentation)){segmentation <- images2matrix(segmentation, mask)}
#   gc()
#
#   # Coerce update on left hand g.formula for parallel fitting
#   g.fo <- update.formula(g.formula, 'Y ~ .')
#
#   # parallel settings
#   plan(cluster, workers = num_cores)
#   handlers(global = TRUE)
#   handlers("pbmcapply")
#   p <- with_progress(progressor(ncol(voxeldata)))
#   future.opt <- list(packages=c('gamlss2'), seed = TRUE)
#   # foreach call
#   models <- foreach(vxlcol = seq_along(voxeldata),
#                     .options.future = future.opt,
#                     .combine=c
#   ) %dofuture% {
#     # fit specific voxel
#     vxl_train_data <- train.data
#     vxl_train_data$Y <- voxeldata[,vxlcol]
#     # if multi tissue add
#     if (!is.null(segmentation)){
#       vxl_train_data$tissue <- segmentation[,vxlcol]}
#     # gamlss
#     g <- gamlss2::gamlss2(formula=g.fo,
#                           data=vxl_train_data,
#                           family=g.family,
#                           light=TRUE,
#                           trace=FALSE,
#                           ...)
#     g$control <- NULL
#     g$family <- g$family[1] # reset via gamlss2:::complete_family()
#     g$vxl <- vxlcol
#     p() # update progressbar
#     list(g)
#   }
#   gc()
#   return(structure(models, class = "vbgamlss"))
# }
#
#
#
#
#
#
#
# ### vbgamlss with pbmclapply ###
#
#
# vbgamlss_pbmclapply <- function(image, mask, g.formula, train.data, g.family=NO,
#                              segmentation=NULL,
#                              num_cores=NULL,
#                              chunk_max_mb=258,
#                              afold=NULL,
#                              ...) {
#
#   if (missing(image)) { stop("image is missing")}
#   if (missing(mask)) { stop("mask is missing")}
#   if (missing(g.formula)) { stop("formula is missing")}
#   if (missing(train.data)) { stop("subjData is missing")}
#
#   if (class(g.formula) != "character") { stop("g.formula class must be character")}
#   if (is.null(num_cores)) {num_cores <- availableCores()}
#
#   # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
#   voxeldata <- images2matrix(image, mask)
#   if (!is.null(segmentation)){segmentation <- images2matrix(segmentation, mask)}
#   gc()
#
#   # subset the dataframe if the input is a fold from CV
#   #   a fold must be a boolean vector of length of number of subjects (image 4th dim)
#   if (!is.null(afold)){
#     if (!is.logical(afold) && !is.integer(afold)) {
#       stop("Error: afold must be either logical or integer.")
#     }
#     voxeldata <- voxeldata[afold,]
#     segmentation <- segmentation[afold,]
#     train.data <- train.data[afold,]
#   }
#
#   # Coerce update on left hand g.formula for parallel fitting
#   g.fo <- update.formula(g.formula, 'Y ~ .')
#
#   # compute chunk size
#   Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
#   # parallel call
#   chunked = as.list(isplitIndices(ncol(voxeldata), chunks=Nchunks))
#   i = 1
#   models <- list()
#   for (ichunk in chunked){
#     cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
#     i<- i+1
#     voxeldata_chunked <- voxeldata[,ichunk]
#     # fit vbgamlss
#     submodels <- pbmclapply(seq_along(voxeldata_chunked),
#                             vbgamlss_parallel_call,
#                             # fnc params
#                             voxeldata_chunked_=voxeldata_chunked,
#                             vxl_train_data=train.data,
#                             segmentation_=segmentation,
#                             g.fo_=g.fo,
#                             g.family_=g.family,
#                             # pbmclapply params
#                             mc.cores=num_cores)
#     models <- c(models, submodels)
#   }
#   gc()
#   return(structure(models, class = "vbgamlss"))
# }
#
# vbgamlss_parallel_call <- function(i,
#                                    voxeldata_chunked_=voxeldata_chunked,
#                                    vxl_train_data=train.data,
#                                    segmentation_=NULL,
#                                    g.fo_=g.fo,
#                                    g.family_=g.family,
#                                    ...) {
#   # fit specific voxel i == voxel column
#   vxl_train_data$Y <- voxeldata_chunked_[,i]
#   # if multi tissue add
#   if (!is.null(segmentation_)){
#     vxl_train_data$tissue <- segmentation_[,i]}
#   # gamlss
#   g <- gamlss2::gamlss2(formula=g.fo_,
#                         data=vxl_train_data,
#                         family=g.family_,
#                         light=TRUE,
#                         trace=FALSE,
#                         ...)
#   g$control <- NULL
#   g$family <- g$family[1] # re-set via gamlss2:::complete_family()
#   g$vxl <- i
#   return(list(g))
# }
#
#
#
#
#
#
#
#
#
# ### vbgamlss stored in temp file, unfinished, very slow ###
#
#
# vbgamlss_ <- function(image, mask, g.formula, train.data, g.family=NO,
#                      segmentation=NULL, afold=NULL, subsample=NULL,
#                      num_cores=NULL, ...) {
#
#   if (missing(image)) { stop("image is missing")}
#   if (missing(mask)) { stop("mask is missing")}
#   if (missing(g.formula)) { stop("formula is missing")}
#   if (missing(train.data)) { stop("subjData is missing")}
#   if (is.null(num_cores)) {num_cores <- availableCores()}
#
#   # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
#   voxeldata <- images2matrix(image, mask)
#   voxeldata <- subset_matrix(voxeldata, afold=afold, subsample = subsample)
#   Nvox <- length(voxeldata)
#   # mk registry
#   registry <- vbgamlss_registry(Nvox)
#   # decompose image in timeseries
#   vbgamlss_decompose(voxeldata, registry$vxls_ts)
#   rm(voxeldata)
#
#   if (!is.null(segmentation)){
#     segmentation <- images2matrix(segmentation, mask)
#     voxeldata <- subset_matrix(segmentation, afold=afold, subsample = subsample)
#     vbgamlss_decompose(segmentation, registry$vxls_segts)
#     rm(segmentation)
#   }
#
#   # Save shared train data
#   registry$shared_data <- train.data
#   registry$formula <- g.formula
#   registry$family <- g.family
#
#   # par call
#   submodels <- pbmclapply(registry$ith,
#                           vbgamlss_fit_call,
#                           # extra params
#                           registry=registry,
#                           # pbmclapply params
#                           mc.cores = num_cores)
#   gc()
#   #return(structure(models, class = "vbgamlss"))
# }
#
#
#
# vbgamlss_registry <- function(Nvxls, env=NULL, wd=NULL){
#   registry <- list()
#   if (is.null(wd)) {registry$wd <- getwd()} else {registry$wd <- wd}
#   if (is.null(env)) {registry$env <- NULL} else {registry$env <- env}
#   # generic container
#   registry$main <- file.path(registry$wd, ".vbgamlss.wd")
#   # model specific
#   registry$path <-file.path(registry$main, paste0('.', rand_names(1), '_', Sys.Date()))
#   # model specific but shared by every voxel
#   registry$shared <- file.path(registry$path, paste0('.shared'))
#   # voxel specific
#   registry$ith <- seq(Nvxls)
#   registry$vxls <- paste0('.', rand_names(Nvxls),'.', registry$ith)
#   registry$vxlspaths <- file.path(registry$path, registry$vxls)
#   # others
#   registry$vxls_ts <- file.path(registry$vxlspaths,
#                                 paste0(registry$vxls, '.vxltimeseries.rds'))
#   registry$vxls_results <- file.path(registry$vxlspaths,
#                                      paste0(registry$vxls, '.model.rds'))
#   registry$vxls_segts <- file.path(registry$vxlspaths,
#                                    paste0(registry$vxls, '.segtimeseries.rds'))
#   # mk paths
#   quite(lapply(registry$vxlspaths, dir.create, recursive = TRUE))
#
#   return(registry)
# }
#
#
#
# vbgamlss_decompose <- function(img, registry_paths){
#   # save parallel
#   noout <- pbmclapply(1:dim(img)[2], # iterate vxl
#                       function(i) {
#                         saveRDS(img[,i], file=registry_paths[[i]])
#                         NULL
#                       },
#                       mc.cores=availableCores())
# }
#
#
#
# vbgamlss_fit_call <- function(i, registry=registry, ...){
#   vxl_train_data <- registry$shared_data
#   vxl_train_data$Y <- readRDS(registry$vxls_ts[[i]])
#   if (!is.null(segmentation)){
#     vxl_train_data$tissue <- readRDS(registry$vxls_segts[[i]])}
#
#   g <- gamlss2::gamlss2(formula=registry$formula,
#                         data=vxl_train_data,
#                         family=registry$family,
#                         light=TRUE,
#                         trace=FALSE,
#                         ...)
#   g$control <- NULL
#   g$family <- g$family[1] # reset via gamlss2:::complete_family()
#   saveRDS(g, file=registry$vxls_results[[i]])
#   NULL
# }
#
#
#
# subset_matrix <- function(obj, afold=NULL, subsample=NULL) {
#   if (!is.null(afold)){
#     if (!is.logical(afold) && !is.integer(afold)) {
#       stop("Error: afold must be either logical or integer vector.")
#     }
#     obj <- obj[afold,]
#
#   }
#
#   # subset the dataframe if a subsampling scheme is provided
#   if (!is.null(subsample)){
#     if (!is.numeric(subsample)) {
#       stop("Error: subsample must be a numeric vector of indeces of length sum(mask>0).")
#     }
#     obj <- obj[,subsample]
#   }
#   return(obj)
# }
#



