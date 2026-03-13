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
#' @param imageframe A dataframe containing the voxel data. Each column represents a voxel/vertex,
#'                   and each row represents a subject.
#' @param g.formula A formula for the GAMLSS2 model.
#' @param train.data A data frame containing subject-level data. Columns must correspond to covariates in `g.formula`.
#' @param g.family A GAMLSS family object. Defaults to `NO` (Normal).
#' @param segmentation Optional, dataframe generated from a mask containing segmentation information (e.g. 1=GM, 2=WM). It must have the shape as `imageframe`.
#' @param segmentation_target Optional, integer, target/label for the segmentation data to subset the analysis (e.g. 1 for GM).
#' @param num_cores Number of cores for parallel processing. If NULL, all available cores are used.
#' @param chunk_max_mb Maximum chunk size in megabytes for processing. Defaults to 64 MB. Increase to 128/256/512 for HPC.
#' @param afold Optional boolean or integer vector to subset imageframe for cross-validation folds.
#' @param debug Logical. If TRUE, enables debug mode to log output in `logdir`. Defaults to FALSE.
#' @param logdir Directory path for saving logs if `debug` is TRUE. Defaults to current working directory.
#' @param cache Logical, cache intermediate results.
#' @param cachedir Character, path for cached artifacts.
#' @param force_ypositivity Logical, shift/clip y to be positive if the GAMLSS distribution family requires it.
#' @param ... Additional arguments passed to the GAMLSS fitting function.
#'
#' @return A list of GAMLSS models, with class `"vbgamlss"`.
#'
#' @details
#' The function performs the following steps:
#' - Checks input data and prepares voxel and segmentation data as needed.
#' - If parallel processing is enabled, splits voxel data into chunks based on memory constraints.
#'   Each chunks are processed sequentially, but voxels/vertex within each chunk are processed in parallel.
#' - Fits a GAMLSS model to each voxel using specified covariates and segmentation data, if available.
#' - Returns a vbgamlss object containing the list of models, one for each voxel.
#'
#' @import future
#' @import doFuture
#' @import progressr
#' @import gamlss2
#' @import itertools
#' @import RhpcBLASctl
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
                      debug=F, # toggle debugging
                      logdir=getwd(), # debug directory
                      cache=F, # toggle save temporary states and force debug=T
                      cachedir=getwd(), # directory caching
                      force_ypositivity=T, # force Y >0
                      warm_start=NULL, # per voxel by param
                      show_progress=T,
                      eps=1e-5,
                      maxit=c(300, 100),
                      ...) {

  # checks
  if (missing(imageframe)) { stop("imageframe is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  # check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}

  # Force character columns to factors
  train.data <- as.data.frame(train.data, stringsAsFactors=TRUE)

  # Cores
  if (is.null(num_cores)) {num_cores <- future::availableCores()}

  # SPEED OPTIMIZATION: Matrix Conversion
  if (!is.data.frame(imageframe) && !is.matrix(imageframe)) {
    stop("Error: imageframe must be a data.frame or matrix!")
  }
  voxeldata <- as.matrix(imageframe)

  # Segmentation if provided
  if (!is.null(segmentation)){
    if (!is.data.frame(segmentation) && !is.matrix(segmentation)) {
      stop("Error: segmentation must be a data.frame or matrix")
    }
    segmentation <- as.matrix(segmentation)
  }
  gc()

  # Subset the imageframe if the input is a fold from CV
  if (!is.null(afold)){
    if (!is.logical(afold) && !is.integer(afold)) {
      stop("Error: afold must be either logical or integer vector.")
    }
    voxeldata <- voxeldata[afold, , drop=FALSE]
    if (!is.null(segmentation)) segmentation <- segmentation[afold, , drop=FALSE]
    train.data <- train.data[afold, , drop=FALSE]
  }

  # ---------------------------------------------------------
  # SAFE PARALLEL SETUP (Master Node)
  # ---------------------------------------------------------
  # Capture the current future plan and guarantee restoration
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(strategy="future::cluster", workers=num_cores)
  options(future.globals.maxSize=20000*1024^2)
  # ---------------------------------------------------------

  if (show_progress) {
    progressr::handlers(global = TRUE)
    progressr::handlers("pbmcapply")
  }
  future.opt <- list(packages=c('gamlss2'), seed = TRUE)

  # Compute chunk size
  Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
  chunked = as.list(itertools::isplitIndices(ncol(voxeldata), chunks=Nchunks))

  # Cache dir
  if (cache) {
    cat(paste0('Caching\n'))
    if (dir.exists(cachedir)){
      if ('.voxchunk.1' %in% list.files(cachedir, all.files = T)){
        cat('Cache directory found\n')
      } else {
        cachedir=file.path(cachedir, paste0('.voxcache.', rand_names(1, l=4)))
        dir.create(cachedir, recursive = T, showWarnings = F)
        cat(paste0('Cache directory: ', cachedir, '\n'))
      }
    } else {
      stop('ERROR: provide a proper/existing path for the cachedir.')
    }
    debug <- T
    logdir <- cachedir
  }

  # Log dir
  if (debug) {
    logdir=file.path(logdir, paste0('.voxlog.', rand_names(1, l=4)))
    dir.create(logdir, recursive = T, showWarnings = F)
  }

  # Parse formula ONCE outside the loop
  g_form_parsed <- as.formula(g.formula)

  # Pre-allocate models list
  models <- vector("list", length(chunked))

  # Loop chunks call
  for (i in seq_along(chunked)) {
    ichunk <- chunked[[i]]

    cat(paste0("Chunk: ", i, "/", Nchunks, "\n"))
    chunk_id = file.path(cachedir, paste0('.voxchunk.', i))

    # Subset with chunk indexes
    voxeldata_chunked <- voxeldata[, ichunk, drop = FALSE]
    if (!is.null(segmentation)){voxelseg_chunked <- segmentation[, ichunk, drop = FALSE]}
    if (!is.null(warm_start)) {warm_start_chunked <- warm_start[ichunk, , drop = FALSE]}

    if (show_progress) { p <- progressr::progressor(ncol(voxeldata_chunked)) }

    # Check cache, continue if submodels already computed
    if (cache && file.exists(chunk_id)){
      cat("Chunk already processed, skipped\n")
      models[[i]] <- readRDS(chunk_id)
      next
    }

    # Parallel chunk loop
    submodels_raw <- foreach::foreach(vxlcol = seq_along(ichunk),
                                      .options.future = future.opt) %dofuture% {

                                        # SAFE WORKER THREAD CONTROL
                                        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
                                          old_blas <- RhpcBLASctl::blas_get_num_procs()
                                          old_omp  <- RhpcBLASctl::omp_get_max_threads()
                                          RhpcBLASctl::blas_set_num_threads(1L)
                                          RhpcBLASctl::omp_set_num_threads(1L)
                                        }

                                        # Use tryCatch to guarantee restoration even if gamlss2 fails
                                        worker_result <- tryCatch({

                                          # 1. Determine valid row indices based on positivity and segmentation target
                                          Y_vxl <- as.numeric(voxeldata_chunked[, vxlcol])
                                          valid_idx <- rep(TRUE, length(Y_vxl))

                                          if (force_ypositivity) {
                                            valid_idx <- valid_idx & (Y_vxl > 0)
                                          }
                                          if (!is.null(segmentation) && !is.null(segmentation_target)) {
                                            valid_idx <- valid_idx & (voxelseg_chunked[, vxlcol] == segmentation_target)
                                          }

                                          # 2. Subset data vectors matching the exact rows kept
                                          vxl_train_data <- train.data[valid_idx, , drop = FALSE]
                                          vxl_train_data$Y <- Y_vxl[valid_idx]

                                          # 3. Subset warm start safely
                                          vxl_start <- NULL
                                          if (!is.null(warm_start)) {
                                            vxl_start <- warm_start_chunked[vxlcol, ]
                                          }

                                          logfile <- if (debug) file.path(logdir, paste0('log.vxl', vxlcol)) else NULL

                                          # GAMLSS with Warm Start implementation
                                          g <- TRY(gamlss2::gamlss2(formula = g_form_parsed,
                                                                    data = vxl_train_data,
                                                                    family = g.family,
                                                                    light = TRUE,
                                                                    trace = FALSE,
                                                                    start = vxl_start,
                                                                    maxit = maxit,
                                                                    control=gamlss2::gamlss2_control(trace = FALSE, eps = eps),
                                                                    ...),
                                                   logfile, save.env.and.stop = F)

                                          if (show_progress) { p() }

                                          # Error handling and deep environment stripping
                                          if (identical(g, NA)) {
                                            return(list(list(vxl = vxlcol, error = TRUE)))
                                          }

                                          g$control <- NULL
                                          g$converged <- g$iterations < maxit[1L]
                                          g$family <- g$family$family
                                          g$vxl <- vxlcol

                                          # Break environment closures to prevent massive IPC memory bloat
                                          if (!is.null(g$terms)) environment(g$terms) <- globalenv()
                                          if (!is.null(g$formula)) environment(g$formula) <- globalenv()
                                          g$call <- NULL
                                          g$fake_formula <- NULL

                                          list(g)

                                        }, finally = {
                                          # This block ALWAYS runs, guaranteeing the worker resets to normal
                                          if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
                                            RhpcBLASctl::blas_set_num_threads(old_blas)
                                            RhpcBLASctl::omp_set_num_threads(old_omp)
                                          }
                                        })

                                        worker_result # Return the safely computed result
                                      }

    # Flatten chunk results instantly
    submodels <- unlist(submodels_raw, recursive = FALSE)

    if (cache){saveRDS(submodels, chunk_id)}

    # Assign directly to pre-allocated list slot
    models[[i]] <- submodels
  }

  gc()
  # Flatten the pre-allocated chunked list into a single model list
  models <- unlist(models, recursive = FALSE)

  return(structure(models, class = "vbgamlss"))
}



# ================== #
#  TESTING VERSIONS  #
# ================== #

# None

# ========== #
# DEPRECATED #
# ========== #

# vbgamlss <- function(imageframe,
#                      g.formula,
#                      train.data,
#                      g.family=NO,
#                      segmentation=NULL,
#                      segmentation_target=NULL,
#                      num_cores=NULL,
#                      chunk_max_mb=64,
#                      afold=NULL,
#                      debug=F, # toggle debugging
#                      logdir=getwd(), # debug directory
#                      cache=F, # save temporary states and force debug=T
#                      cachedir=getwd(), # debug directory
#                      force_ypositivity=T, # force Y >0
#                      ...) {
#
#   # checks
#   if (missing(imageframe)) { stop("imageframe is missing")}
#   if (missing(g.formula)) { stop("formula is missing")}
#   check_formula_LHS(g.formula)
#   if (missing(train.data)) { stop("subjData is missing")}
#
#
#   # Force character columns to factors
#   train.data <- as.data.frame(train.data, stringsAsFactors=TRUE) # unclass(train.data)
#
#
#   # Cores
#   if (is.null(num_cores)) {num_cores <- availableCores()}
#
#   # Input image: subjects (4th dim) x voxels (columns)
#   if (!is.data.frame(imageframe) && !is.matrix(imageframe)) {
#     stop("Error: imageframe must be a data.frame or matrix! Use images2matrix(image, mask) to convert image to matrix.")
#   }
#   voxeldata <- as.data.frame(imageframe)
#
#
#   # Segmentation if provided
#   if (!is.null(segmentation)){
#     if (!is.data.frame(segmentation)) {stop("Error: segmentaion must be a data.frame")}
#   }
#   gc()
#
#
#   # Subset the imageframe if the input is a fold from CV
#   #   a fold must be a boolean vector of length of number of subjects (image 4th dim)
#   if (!is.null(afold)){
#     if (!is.logical(afold) && !is.integer(afold)) {
#       stop("Error: afold must be either logical or integer vector.")
#     }
#     voxeldata <- voxeldata[afold,]
#     segmentation <- segmentation[afold,]
#     train.data <- train.data[afold,]
#   }
#
#   # Parallel settings
#   future::plan(strategy="future::cluster", workers=num_cores) # ->   with(plan(multisession), local = TRUE)
#   options(future.globals.maxSize=20000*1024^2)
#   progressr::handlers(global = TRUE)
#   progressr::handlers("pbmcapply")
#   p <- with_progress(progressor(ncol(voxeldata)))
#   future.opt <- list(packages=c('gamlss2'), seed = TRUE)
#
#
#   # Compute chunk size
#   Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
#   chunked = as.list(itertools::isplitIndices(ncol(voxeldata), chunks=Nchunks))
#
#
#   # Cache dir
#   if (cache) {
#     cat(paste0('Caching'), fill=T)
#     # exist?
#     if (dir.exists(cachedir)){
#       # is already a cache dir?
#       if ('.voxchunk.1' %in% list.files(cachedir, all.files = T)){
#         cat('Cache directory found', fill=T)
#       } else {
#         # make a new folder to store cache in the cachedir
#         cachedir=file.path(cachedir, paste0('.voxcache.', rand_names(1, l=4)))
#         dir.create(cachedir, recursive = T, showWarnings = F)
#         cat(paste0('Cache directory: ', cachedir), fill=T)
#       }
#     } else {
#       stop('ERROR: provide a proper/existing path for the cachedir.')
#     }
#     # force debugging while in cache mode
#     debug <- T
#     logdir <- cachedir
#   }
#
#
#   # Log dir
#   if (debug) {
#     logdir=file.path(logdir, paste0('.voxlog.', rand_names(1, l=4)))
#     dir.create(logdir, recursive = T, showWarnings = F)
#   }
#
#
#   # Loop chunks call
#   i = 1
#   models <- list()
#   for (ichunk in chunked){
#
#     # Counters, names, etc..
#     cat(paste0("Chunk: ",i,"/", Nchunks), fill=T)
#     chunk_id = file.path(cachedir, paste0('.voxchunk.', i))
#     i <- i+1
#
#
#     # Subset with chunk indexes
#     voxeldata_chunked <- voxeldata[,ichunk]
#     if (!is.null(segmentation)){voxelseg_chunked <- segmentation[,ichunk]}
#
#
#     # Set up progress bar per chunk
#     p <- progressor(ncol(voxeldata_chunked))
#
#
#     # Check cache, continue if submodels already computed
#     if (cache){
#       if (file.exists(chunk_id)){
#         cat("Chunk already processed, skipped", fill=T)
#         submodels <- readRDS(chunk_id)
#         models <- c(models, submodels)
#         next
#       }
#     }
#
#
#     # Parallel call to fit vbgamlss
#     submodels <- foreach(vxlcol = seq_along(voxeldata_chunked),
#                          .options.future = future.opt,
#                          .combine=c
#     ) %dofuture% {
#       # Fit specific voxel
#       vxl_train_data <- train.data
#       vxl_train_data$Y <- as.numeric(voxeldata_chunked[,vxlcol])
#       # If multi tissue fit
#       if (!is.null(segmentation)){
#         vxl_train_data$tissue <- voxelseg_chunked[,vxlcol]}
#       if (! is.null(segmentation_target)){
#         vxl_train_data <- vxl_train_data[vxl_train_data$tissue == segmentation_target,]
#       }
#
#
#       # Debug
#       if (debug) {logfile=file.path(logdir, paste0('log.vxl', vxlcol))} else {logfile=NULL}
#
#
#       # Force Y to be strictly non-negative & !=0
#       if (force_ypositivity){
#         non_negative = vxl_train_data$Y > 0
#         vxl_train_data <- vxl_train_data[non_negative,]
#       }
#
#       # GAMLSS
#       g <- TRY(gamlss2::gamlss2(formula=as.formula(g.formula),
#                                 data=vxl_train_data,
#                                 family=g.family,
#                                 light=TRUE,
#                                 trace=FALSE,
#                                 maxit = c(300, 100),
#                                 ...),
#                logfile, save.env.and.stop=F)
#       g$control <- NULL # remove control to save space
#       g$family <- g$family$family # to re-set via gamlss2:::complete_family()
#       g$vxl <- vxlcol
#       p() # update progress bar
#       list(g)
#     }
#
#
#     # Save chunk
#     if (cache){saveRDS(submodels, chunk_id)}
#
#
#     # Append results to models
#     models <- c(models, submodels)
#
#   }
#   gc()
#   return(structure(models, class = "vbgamlss"))
# }
#
#

