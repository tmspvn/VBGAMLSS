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
#' @param force_constraints list, force voxels outside constraints to be excluded from the fit:  Constr1 <= voxels <= Constr2. Use NULL to exclude constraints. Defaults to c(1e-8, +Inf)
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
                     #afold=NULL,
                     debug=F, # toggle debugging (if T, retains cache dir at the end)
                     cachedir=getwd(), # directory caching
                     force_constraints=c(1e-8, +Inf),
                     warm_start=NULL, # per voxel by param
                     show_progress=T,
                     eps=1e-5,
                     maxit=c(100, 33),
                     future_plan_strategy = "future.mirai::mirai_cluster",
                     save_model=NULL,
                     deep_stripping = TRUE, # remove .Enviroments form fitted models
                     ...) {


  # -----------------------------------------------------------
  # CHECKS

  if (missing(imageframe)) { stop("imageframe is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  # check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}

  if (nrow(imageframe) != nrow(train.data)) {
    stop("Error: imageframe and train.data have different row counts. Subjects must be strictly aligned.")
  }

  # Force character columns to factors
  train.data <- as.data.frame(train.data, stringsAsFactors=TRUE)

  # Cores
  if (is.null(num_cores)) {num_cores <- future::availableCores()}

  # Matrix Conversion for speed
  if (!is.data.frame(imageframe) && !is.matrix(imageframe)) {
    stop("Error: imageframe must be a data.frame or matrix!")
  }
  voxeldata <- as.matrix(imageframe)
  nsub      <- dim(voxeldata)[1]
  nvox      <- dim(voxeldata)[2]
  rm(imageframe)

  # Segmentation if provided
  if (!is.null(segmentation)){
    if (!is.data.frame(segmentation) && !is.matrix(segmentation)) {
      stop("Error: segmentation must be a data.frame or matrix")
    }
    segmentation <- as.matrix(segmentation)
  }

  # Subset the imageframe if the input is a fold from CV
  # if (!is.null(afold)){
  #   if (!is.logical(afold) && !is.integer(afold)) {
  #     stop("Error: afold must be either logical or integer vector.")
  #   }
  #   voxeldata <- voxeldata[afold, , drop=FALSE]
  #   if (!is.null(segmentation)) segmentation <- segmentation[afold, , drop=FALSE]
  #   if (!is.null(warm_start)) warm_start <- warm_start[afold, , drop=FALSE]
  #   train.data <- train.data[afold, , drop=FALSE]
  # }

  gc()



  # -----------------------------------------------------------
  # PARALLEL SETUP

  if (any(grepl("mirai", future_plan_strategy))) {
    mirai::daemons(num_cores)
    future::plan(strategy=future_plan_strategy)
    show_progress <- F
  } else {
    future::plan(strategy=future_plan_strategy, workers=num_cores)
  }
  options(future.globals.maxSize=20000*1024^2)
  # make sure to avoid exporting massive stuff
  future.opt <- list(packages = c('gamlss2'),
                     seed     = TRUE,
                     globals  = structure(TRUE, ignore = "voxeldata")
  )

  # progressr
  if (show_progress) {
    progressr::handlers(global = TRUE)
    progressr::handlers("pbmcapply")
  }

  # get blas omp values?
  master_blas <- RhpcBLASctl::blas_get_num_procs()
  master_omp  <- RhpcBLASctl::omp_get_max_threads()



  # -----------------------------------------------------------
  # CACHING SETUP

  if (debug) {
    cat('Debug=TRUE, cache directory and logs will be retained after execution.\n')
  }
  is_new_cache <- FALSE

  # Define the expected path for the registry
  registry_path <- file.path(cachedir, '.vbgamlss.registry')

  # If the file exists, resume. If not, new cache
  if (file.exists(registry_path)) {
    cat('Cache directory and registry found. Resuming from external/existing cache\n')
    registry <- qs2::qs_read(registry_path)
    logdir <- file.path(cachedir, '.voxlog')
    voxfits_subdir_path <- file.path(cachedir, '.voxfits')

    # Check if registry and fitted models are congruent
    registry$full_paths <- file.path(cachedir, registry$relative_paths)
    cat('Recheck of existing fits convergence\n')
    isconverged <- future.apply::future_vapply(registry$full_paths,
                                               function(x) {
                                                 if (file.exists(x)) {
                                                   m <- qs2::qs_read(x, nthreads = 1)
                                                   if (m$converged) {return(1)} else {return(2)}
                                                 } else {return(404)
                                                 }},
                                               FUN.VALUE = numeric(1))
    cat('Of', length(registry$fitted), 'voxels:',
        sum(registry$fitted),   'were fitted,',
        sum(isconverged == 1),  'converged,',
        sum(isconverged == 2),  'did not, and',
        sum(isconverged == 404),'went missing. \n')
    registry$fitted <- isconverged == 1

  } else {
    # Flag that this is a new cache created in this run
    is_new_cache <- TRUE

    # Make main cache folder
    fit_rand_id <- rand_names(1, l=4)
    cachedir <- file.path(cachedir, paste0('.vbgamlss.cache.', fit_rand_id))
    dir.create(cachedir, recursive = T, showWarnings = F)

    # Make folder for debugging
    logdir <- file.path(cachedir, '.voxlog')
    dir.create(logdir, recursive = T, showWarnings = F)

    # Make cache subfolder for voxfits
    voxfits_subdir_name <- '.voxfits'
    voxfits_subdir_path <- file.path(cachedir, voxfits_subdir_name)
    dir.create(voxfits_subdir_path, recursive = T, showWarnings = F)

    # Make registry
    registry <- data.frame(voxel            = 1:nvox,
                           fitted           = logical(nvox),
                           converged        = logical(nvox),
                           relative_paths   = file.path(voxfits_subdir_name,
                                                        paste0(".vbgamlss.voxel.", 1:nvox)),
                           full_paths       = file.path(voxfits_subdir_path,
                                                        paste0(".vbgamlss.voxel.", 1:nvox)),
                           stringsAsFactors = FALSE)

    # Update registry_path to the new subfolder and save
    registry_path <- file.path(cachedir, '.vbgamlss.registry')
    qs2::qs_save(registry, registry_path)
    cat(paste0('Cache directory created: ', cachedir, '\n'))
  }



  # ---------------------------------------------------------
  # LARGE IMAGE CHUNKING & PREPPING ROUTINE

  # Parse formula ONCE outside the loop
  g_form_parsed <- as.formula(g.formula)

  # Compute chunk size
  Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
  chunked <- as.list(itertools::isplitIndices(ncol(voxeldata), chunks=Nchunks))

  # Loop chunks call
  start.time   <- Sys.time()
  current.time <- start.time
  for (i in seq_along(chunked)) {

    # Tracking
    ichunk <- chunked[[i]]
    cat(paste0("Chunk: ", i, "/", Nchunks, format(Sys.time(), " (started: %X, %d %b %Y)"),
               ' [elapsed: ', format(difftime(current.time, start.time, units = "hours"),
                                     digits = 3), ']',"\n"))

    # Explicitly match the subsetting
    registry_rows <- match(ichunk, registry$voxel)
    is_unfitted   <- registry$fitted[registry_rows] %in% FALSE
    ichunk        <- ichunk[is_unfitted]

    # If the entire chunk is already fitted, skip to the next chunk immediately
    if (length(ichunk) == 0) {
      cat("All voxels in this chunk are already fitted. Skipping...\n")
      next
    }

    # Subset datasets with chunk indexes
    voxeldata_chunked <- voxeldata[, ichunk, drop = FALSE] # sub x vxl
    registry_chunked  <- registry[ichunk, , drop = FALSE] # vxl x [voxel, fitted, relative_paths, full_paths]
    if (!is.null(segmentation))
    {voxelseg_chunked <- segmentation[, ichunk, drop = FALSE]} # 1 x vxl


    # ---------------------------------------------------------
    # VOXEL PREPROCESSING

    # Internal function for preparing voxel for fitting
    prep_func <- function(idvxl) {
      # Registry for voxel
      vxl_reg_entry <- as.list(registry_chunked[idvxl, ])

      # Prepare voxel data
      Y_vxl     <- as.numeric(voxeldata_chunked[, idvxl])
      valid_idx <- rep(TRUE, length(Y_vxl))

      if (!is.null(force_constraints)) {
        valid_idx <- valid_idx & !is.na(Y_vxl) &
          (Y_vxl >= force_constraints[1]) &
          (Y_vxl <= force_constraints[2])
        # disrespectfully crash everything if constrains are exaggerated
        if (sum(valid_idx) < 5) {
          stop(paste0('force_constraints left less than 5 observations for voxel ',
                      vxl_reg_entry$voxel))
        }
      }

      # propagate indexes to segmentation
      if (!is.null(segmentation) && !is.null(segmentation_target)) {
        valid_idx <- valid_idx & (voxelseg_chunked[, idvxl] == segmentation_target)
      }

      # Subset data vectors matching the exact rows kept
      vxl_train_data   <- train.data[valid_idx, , drop = FALSE]
      vxl_train_data$Y <- Y_vxl[valid_idx]
      vxl_warm_start <- NULL
      if (!is.null(warm_start)) {
        vxl_warm_start <- warm_start[valid_idx, ]
      }

      # Return only what the worker needs
      list(
        vxlcol         = vxl_reg_entry$voxel,
        data           = vxl_train_data,
        start_params   = vxl_warm_start,
        registry_entry = vxl_reg_entry
      )
    }

    # Prepare input list for each voxel
    prepared_voxel_data <- future.apply::future_lapply(seq_along(ichunk), prep_func)

    # Clean up chunk matrices to free memory before parallel execution
    rm(voxeldata_chunked)
    rm(registry_chunked)
    if (!is.null(segmentation)) rm(voxelseg_chunked)
    gc()


    # ---------------------------------------------------------
    # PARALLEL PROCESSING ROUTINE

    # Track progress per chunk
    if (show_progress) { p <- progressr::progressor(length(prepared_voxel_data)) }

    # Parallel chunk loop iterating over the prepared list
    cat('Processing ', length(ichunk), ' voxels of', nvox,'\n')
    submodels <- foreach::foreach(vxl_item = prepared_voxel_data,
                                  .options.future = future.opt)  %dofuture% {

                                    # WORKER THREAD CONTROL
                                    if (RhpcBLASctl::blas_get_num_procs() > 1L)
                                          {RhpcBLASctl::blas_set_num_threads(1L)}
                                    if (RhpcBLASctl::omp_get_num_procs() > 1L)
                                          {RhpcBLASctl::omp_set_num_threads(1L)}

                                    registry_entry <- vxl_item$registry_entry
                                    vxlcol         <- registry_entry$voxel

                                    logfile <- NULL
                                    if (!is.null(logdir))
                                          {logfile <- file.path(logdir, paste0('log.vxl', vxlcol))}

                                    # GAMLSS fit using pre-assembled data
                                    g <- TRY(gamlss2::gamlss2(formula = g_form_parsed,
                                                              data    = vxl_item$data,
                                                              family  = g.family,
                                                              start   = vxl_item$start_params,
                                                              maxit   = maxit,
                                                              control = gamlss2::gamlss2_control(trace = FALSE,
                                                                                                 light = TRUE,
                                                                                                 eps  = eps),
                                                              ...),
                                             logfile, save.env.and.stop = F)

                                    if (show_progress) { p() }

                                    # Error handling and deep environment stripping
                                    error = FALSE
                                    if (identical(g, NA)) {
                                      error <- TRUE
                                      g     <- list(vxl = vxlcol, error = TRUE, converged = F)

                                    } else {
                                      # Good fit, strip extras
                                      g$control      <- NULL
                                      g$converged    <- g$iterations < maxit[1L]
                                      g$family       <- g$family$family
                                      g$vxl          <- vxlcol
                                      if (deep_stripping){
                                        g <- deep_env_stripping(g)}
                                    }

                                    # Save on disk to unload master
                                    qs2::qs_save(g, registry_entry$full_paths)
                                    registry_entry$converged <- g$converged
                                    registry_entry$fitted    <- T

                                    # But if error
                                    if (error){
                                      registry_entry$converged <- F
                                      registry_entry$fitted    <- F
                                    }

                                    # Return registry entry
                                    return(registry_entry)

                                    # End parallel loop
                                  }
    # End chunking routine
    gc()

    # Update registry
    worker_returned_entry_voxels    <- vapply(submodels, function(x) x$voxel, numeric(1))
    worker_returned_entry_fitted    <- vapply(submodels, function(x) x$fitted, logical(1))
    worker_returned_entry_converged <- vapply(submodels, function(x) x$converged, logical(1))

    # Update the master registry in place
    match_idx <- match(worker_returned_entry_voxels, registry$voxel)
    registry$fitted[match_idx]    <- worker_returned_entry_fitted
    registry$converged[match_idx] <- worker_returned_entry_converged

    # Save the updated master registry
    qs2::qs_save(registry, registry_path)

    gc()
    current.time <- Sys.time()
  }


  # ---------------------------------------------------------
  # CLOSING
  gc()

  # Reverse blas and openmp threads control, probably useless
  if (RhpcBLASctl::blas_get_num_procs() != master_blas)
        {RhpcBLASctl::blas_set_num_threads(master_blas)}
  if (RhpcBLASctl::omp_get_num_procs() != master_omp)
        {RhpcBLASctl::omp_set_num_threads(master_omp)}


  # Report of success
  isconverged <- future.apply::future_vapply(registry$full_paths,
                                             function(x) {
                                               if (file.exists(x)) {
                                                 m <- qs2::qs_read(x, nthreads = 1)
                                                 if (m$converged) {return(1)} else {return(2)}
                                               } else {return(404)
                                               }},
                                             FUN.VALUE = numeric(1))
  cat('\nOf', length(registry$fitted), 'voxels:',
      sum(registry$fitted),   'were fitted,',
      sum(isconverged == 1),  'converged,',
      sum(isconverged == 2),  'did not, and',
      sum(isconverged == 404),'went missing. \n')


  # Aggregating
  cat("Aggregating individual voxel models\n")
  agg.start.time <- Sys.time()
  models <- vector("list", nvox)

  for (i in seq_along(registry$full_paths)) {
    # read as raw bytes
    f_size <- file.info(registry$full_paths[i])$size
    models[[i]] <- readBin(registry$full_paths[i], what = "raw", n = f_size)
    # Release memory or OOM
    if (i %% 1e4 == 0 || i == nvox) {gc(verbose = FALSE)}
  }
  models <- structure(models, class = "vbgamlss")


  # Save
  if (! is.null(save_model)) {
    qs2::qs_save(models,
                 file = paste0(save_model, '.vbgamlss'),
                 compress_level = 0L) # uncompressed
    cat('Model saved: ', paste0(save_model, '.vbgamlss'), '\n')
  }


  # Cleanup
  if (!debug) {
    if (is_new_cache &&                                                         # is newly made, so not passed?
        grepl("\\.vbgamlss\\.cache", cachedir) &&                               # is names .vbgamlss.cache
        normalizePath(cachedir, mustWork = FALSE) != normalizePath(getwd(),     # is the path != from pwd?
                                                                   mustWork = FALSE)) {
      # Remove cache
      unlink(cachedir, recursive = TRUE)
      }}

  # Bye bye
  cat(paste0('Completed in ', format(difftime(current.time, start.time, units = "hours"), digits = 3),"\n"))
  return(models)
}



#' Define the S3 subsetting method, here it would cool to include restore_family..
#'
#' @rawNamespace S3method("[[", vbgamlss)
`[[.vbgamlss` <- function(x, i, ...) {
  raw_bytes <- unclass(x)[[i, ...]]
  if (is.null(raw_bytes)) return(NULL)
  qs2::qs_deserialize(raw_bytes)
}


# Deep environment stripping
deep_env_stripping <- function(model) {
  # rm redundant backup
  attr(model$xterms, "terms") <- NULL
  # rm top-level environment tethers
  attr(model$terms, ".Environment") <- globalenv()
  attr(model$formula, ".Environment") <- globalenv()
  attr(model$xterms, ".Environment") <- globalenv()
  # rm environments inside individual parameter terms (mu, sigma, etc.)
  for (param in names(model$terms)) {
    environment(model$terms[[param]]) <- globalenv()
    attr(model$terms[[param]], ".Environment") <- globalenv()
  }
  # do the exact same for the formulas to prevent mirroring
  for (param in names(model$formula)) {
    environment(model$formula[[param]]) <- globalenv()
    attr(model$formula[[param]], ".Environment") <- globalenv()
  }
  return(model)
}



# ================== #
#  TESTING VERSIONS  #
# ================== #


# ========== #
# DEPRECATED #
# ========== #
















