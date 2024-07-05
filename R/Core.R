################################  Core  ########################################
# install.packages("gamlss2",
#                  repos = c("https://gamlss-dev.R-universe.dev",
#                   "https://cloud.R-project.org"))
# devtools::install_github('ANTsX/ANTsR')





#' @export
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

  # Coerce update on left hand g.formula for parallel fitting
  g.fo <- update.formula(g.formula, 'Y ~ .')

  # parallel settings
  plan(cluster, workers = num_cores, rscript_libs = .libPaths())
  handlers(global = TRUE)
  handlers("pbmcapply")
  p <- with_progress(progressor(ncol(voxeldata)))
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
