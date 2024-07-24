################################  Core  ########################################
# install.packages("gamlss2",
#                  repos = c("https://gamlss-dev.R-universe.dev",
#                   "https://cloud.R-project.org"))
# devtools::install_github('ANTsX/ANTsR')




### vbgamlss with chunking ###

#' @export
vbgamlss_chunks <- function(image, mask, g.formula, train.data, g.family=NO,
                            segmentation=NULL,
                            num_cores=NULL,
                            chunk_max_mb=64,
                            afold=NULL,
                            subsample=NULL,
                            ...) {

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
  voxeldata <- images2matrix(image, mask)
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
  # parallel call
  chunked = as.list(isplitIndices(ncol(voxeldata), chunks=Nchunks))
  i = 1
  models <- list()
  for (ichunk in chunked){
    cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
    i<- i+1
    voxeldata_chunked <- voxeldata[,ichunk]
    if (!is.null(segmentation)){voxelseg_chunked <- segmentation[,ichunk]}
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
        vxl_train_data$tissue <- voxelseg_chunked[,vxlcol]}
      # gamlss
      g <- gamlss2::gamlss2(formula=as.formula(g.formula),
                            data=vxl_train_data,
                            family=g.family,
                            light=TRUE,
                            trace=FALSE,
                            ...)
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



### vbgamlss all loaded in mem ###

#' @export
vbgamlss_unchunked <- function(image, mask, g.formula, train.data, g.family=NO,
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
  future.opt <- list(packages=c('gamlss2'), seed = TRUE)
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


### vbgamlss with pbmclapply ###

#' @export
vbgamlss_pbmclapply <- function(image, mask, g.formula, train.data, g.family=NO,
                             segmentation=NULL,
                             num_cores=NULL,
                             chunk_max_mb=258,
                             afold=NULL,
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

  # subset the dataframe if the input is a fold from CV
  #   a fold must be a boolean vector of length of number of subjects (image 4th dim)
  if (!is.null(afold)){
    if (!is.logical(afold) && !is.integer(afold)) {
      stop("Error: afold must be either logical or integer.")
    }
    voxeldata <- voxeldata[afold,]
    segmentation <- segmentation[afold,]
    train.data <- train.data[afold,]
  }

  # Coerce update on left hand g.formula for parallel fitting
  g.fo <- update.formula(g.formula, 'Y ~ .')

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
    # fit vbgamlss
    submodels <- pbmclapply(seq_along(voxeldata_chunked),
                            vbgamlss_parallel_call,
                            # fnc params
                            voxeldata_chunked_=voxeldata_chunked,
                            vxl_train_data=train.data,
                            segmentation_=segmentation,
                            g.fo_=g.fo,
                            g.family_=g.family,
                            # pbmclapply params
                            mc.cores=num_cores)
    models <- c(models, submodels)
  }
  gc()
  return(structure(models, class = "vbgamlss"))
}

vbgamlss_parallel_call <- function(i,
                                   voxeldata_chunked_=voxeldata_chunked,
                                   vxl_train_data=train.data,
                                   segmentation_=NULL,
                                   g.fo_=g.fo,
                                   g.family_=g.family,
                                   ...) {
  # fit specific voxel i == voxel column
  vxl_train_data$Y <- voxeldata_chunked_[,i]
  # if multi tissue add
  if (!is.null(segmentation_)){
    vxl_train_data$tissue <- segmentation_[,i]}
  # gamlss
  g <- gamlss2::gamlss2(formula=g.fo_,
                        data=vxl_train_data,
                        family=g.family_,
                        light=TRUE,
                        trace=FALSE,
                        ...)
  g$control <- NULL
  g$family <- g$family[1] # re-set via gamlss2:::complete_family()
  g$vxl <- i
  return(list(g))
}


### vbgamlss stored in temp file, unfinished, slow ###

#' @export
vbgamlss_ <- function(image, mask, g.formula, train.data, g.family=NO,
                     segmentation=NULL, afold=NULL, subsample=NULL,
                     num_cores=NULL, ...) {

  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  if (missing(train.data)) { stop("subjData is missing")}
  if (is.null(num_cores)) {num_cores <- availableCores()}

  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  voxeldata <- images2matrix(image, mask)
  voxeldata <- subset_matrix(voxeldata, afold=afold, subsample = subsample)
  Nvox <- length(voxeldata)
  # mk registry
  registry <- vbgamlss_registry(Nvox)
  # decompose image in timeseries
  vbgamlss_decompose(voxeldata, registry$vxls_ts)
  rm(voxeldata)

  if (!is.null(segmentation)){
    segmentation <- images2matrix(segmentation, mask)
    voxeldata <- subset_matrix(segmentation, afold=afold, subsample = subsample)
    vbgamlss_decompose(segmentation, registry$vxls_segts)
    rm(segmentation)
  }

  # Save shared train data
  registry$shared_data <- train.data
  registry$formula <- g.formula
  registry$family <- g.family

  # par call
  submodels <- pbmclapply(registry$ith,
                          vbgamlss_fit_call,
                          # extra params
                          registry=registry,
                          # pbmclapply params
                          mc.cores = num_cores)
  gc()
  #return(structure(models, class = "vbgamlss"))
}



vbgamlss_registry <- function(Nvxls, env=NULL, wd=NULL){
  registry <- list()
  if (is.null(wd)) {registry$wd <- getwd()} else {registry$wd <- wd}
  if (is.null(env)) {registry$env <- NULL} else {registry$env <- env}
  # generic container
  registry$main <- file.path(registry$wd, ".vbgamlss.wd")
  # model specific
  registry$path <-file.path(registry$main, paste0('.', rand_names(1), '_', Sys.Date()))
  # model specific but shared by every voxel
  registry$shared <- file.path(registry$path, paste0('.shared'))
  # voxel specific
  registry$ith <- seq(Nvxls)
  registry$vxls <- paste0('.', rand_names(Nvxls),'.', registry$ith)
  registry$vxlspaths <- file.path(registry$path, registry$vxls)
  # others
  registry$vxls_ts <- file.path(registry$vxlspaths,
                                paste0(registry$vxls, '.vxltimeseries.rds'))
  registry$vxls_results <- file.path(registry$vxlspaths,
                                     paste0(registry$vxls, '.model.rds'))
  registry$vxls_segts <- file.path(registry$vxlspaths,
                                   paste0(registry$vxls, '.segtimeseries.rds'))
  # mk paths
  quite(lapply(registry$vxlspaths, dir.create, recursive = TRUE))

  return(registry)
}



vbgamlss_decompose <- function(img, registry_paths){
  # save parallel
  noout <- pbmclapply(1:dim(img)[2], # iterate vxl
                      function(i) {
                        saveRDS(img[,i], file=registry_paths[[i]])
                        NULL
                      },
                      mc.cores=availableCores())
}



vbgamlss_fit_call <- function(i, registry=registry, ...){
  vxl_train_data <- registry$shared_data
  vxl_train_data$Y <- readRDS(registry$vxls_ts[[i]])
  if (!is.null(segmentation)){
    vxl_train_data$tissue <- readRDS(registry$vxls_segts[[i]])}

  g <- gamlss2::gamlss2(formula=registry$formula,
                        data=vxl_train_data,
                        family=registry$family,
                        light=TRUE,
                        trace=FALSE,
                        ...)
  g$control <- NULL
  g$family <- g$family[1] # reset via gamlss2:::complete_family()
  saveRDS(g, file=registry$vxls_results[[i]])
  NULL
}



subset_matrix <- function(obj, afold=NULL, subsample=NULL) {
  if (!is.null(afold)){
    if (!is.logical(afold) && !is.integer(afold)) {
      stop("Error: afold must be either logical or integer vector.")
    }
    obj <- obj[afold,]

  }

  # subset the dataframe if a subsampling scheme is provided
  if (!is.null(subsample)){
    if (!is.numeric(subsample)) {
      stop("Error: subsample must be a numeric vector of indeces of length sum(mask>0).")
    }
    obj <- obj[,subsample]
  }
  return(obj)
}




