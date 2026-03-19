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
                     debug=F, # toggle debugging and force cache=T
                     cache=F, # toggle save temporary states
                     cachedir=getwd(), # directory caching
                     force_ypositivity=T, # force Y >0
                     warm_start=NULL, # per voxel by param
                     show_progress=T,
                     eps=1e-5,
                     maxit=c(100, 33),
                     future_plan_strategy = "future.mirai::mirai_cluster",
                     save_model=NULL,
                     remove_voxel_fits_from_cache = T, # if chaching=T, specify if keepign the individual voxels fitsa in memory
                     recheck_existing_fits=F,
                     ...) {


  # -----------------------------------------------------------
  # CHECKS

  if (missing(imageframe)) { stop("imageframe is missing")}
  if (missing(g.formula)) { stop("formula is missing")}
  # check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}

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
  if (!is.null(afold)){
    if (!is.logical(afold) && !is.integer(afold)) {
      stop("Error: afold must be either logical or integer vector.")
    }
    voxeldata <- voxeldata[afold, , drop=FALSE]
    if (!is.null(segmentation)) segmentation <- segmentation[afold, , drop=FALSE]
    train.data <- train.data[afold, , drop=FALSE]
  }

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
  # CHACHING SETUP

  # If debug, force caching and log dir with caching directory
  if (debug) {
    cache <- T
    remove_voxel_fits_from_cache <- F
    warning('Debug=TRUE, forcing chace=TRUE and remove_voxel_fits_from_cache=FALSE as well.\n')
    if (remove_voxel_fits_from_cache) {
      warning('remove_voxel_fits_from_cache = TRUE, cached individual voxel fits will be removed. Set to FALSE if you would like to keep them. \n')
    }
  }

  # Cache dir
  if (cache) {
    cat(paste0('Caching\n'))
    if (dir.exists(cachedir)){
      if ('.vbgamlss.registry' %in% list.files(cachedir, all.files = T)){
        cat('Cache directory found\n')
        registry_path       <- file.path(cachedir, '.vbgamlss.registry')
        registry            <- qs2::qs_read(registry_path)

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
            sum(registry$fitted), 'were fitted,',
            sum(isconverged == 1), 'converged,',
            sum(isconverged == 2), 'did not, and',
            sum(isconverged == 404),'went missing. \n')
        registry$fitted <- isconverged == 1
        }
      } else {
        # make main cache fold
        fit_rand_id <- rand_names(1, l=4)
        cachedir <- file.path(cachedir, paste0('.vbgamlss.cache.', fit_rand_id))
        dir.create(cachedir, recursive = T, showWarnings = F)
        # make cache subfold for voxfits
        voxfits_subdir_name <- '.voxfits'
        voxfits_subdir_path <- file.path(cachedir, voxfits_subdir_name)
        dir.create(voxfits_subdir_path, recursive = T, showWarnings = F)
        # make registry
        registry <- data.frame(voxel            = 1:nvox,
                               fitted           = logical(nvox),
                               converged        = logical(nvox),
                               relative_paths   = file.path(voxfits_subdir_name, paste0(".vbgamlss.voxel.", 1:nvox)),
                               full_paths       = file.path(voxfits_subdir_path, paste0(".vbgamlss.voxel.", 1:nvox)),
                               stringsAsFactors = FALSE)
        registry_path <- file.path(cachedir, '.vbgamlss.registry')
        qs2::qs_save(registry, registry_path)
        cat(paste0('Cache directory: ', cachedir, '\n'))
      }
    } else {
      stop('ERROR: provide a proper/existing path for the cachedir.')
    }


  # Finish setting up debugging path
  if (debug && ! dir.exists(cachedir)) {
    logdir <- cachedir
    logdir <- file.path(logdir, '.voxlog')
    dir.create(logdir, recursive = T, showWarnings = F)
  } else {
    logdir <- file.path(cachedir, '.voxlog')
  }


  # ---------------------------------------------------------
  # LARGE IMAGE CHUNKING & PREPPING ROUTINE

  # Pre-allocate models list if not caching
  if (! cache){models <- vector("list", length(chunked))}

  # Parse formula ONCE outside the loop
  g_form_parsed <- as.formula(g.formula)

  # Compute chunk size
  Nchunks <- estimate_nchunks(voxeldata, chunk_max_Mb=chunk_max_mb)
  chunked <- as.list(itertools::isplitIndices(ncol(voxeldata), chunks=Nchunks))

  # Loop chunks call
  start.time <- Sys.time()
  current.time <- start.time
  for (i in seq_along(chunked)) {

    # Tracking
    ichunk <- chunked[[i]]
    cat(paste0("Chunk: ", i, "/", Nchunks, format(Sys.time(), " (started: %X, %d %b %Y)"),
               ' [elapsed: ', format(current.time - start.time, digits = 3), ']',"\n"))


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
    if (!is.null(warm_start))
      {warm_start_chunked <- warm_start[ichunk, , drop = FALSE]} # vxl x [mu, sigma, nu, tau]

    # ---------------------------------------------------------
    # VOXEL PREPROCESSING

    # Internal function for preparing voxel for fitting
    prep_func <- function(idvxl) {
      # Registry for voxel
      vxl_reg_entry <- as.list(registry_chunked[idvxl, ])

      # Prepare voxel data
      Y_vxl     <- as.numeric(voxeldata_chunked[, idvxl])
      valid_idx <- rep(TRUE, length(Y_vxl))

      if (force_ypositivity) {
        valid_idx <- valid_idx & (Y_vxl > 0)
      }
      if (!is.null(segmentation) && !is.null(segmentation_target)) {
        valid_idx <- valid_idx & (voxelseg_chunked[, idvxl] == segmentation_target)
      }

      # Subset data vectors matching the exact rows kept
      vxl_train_data   <- train.data[valid_idx, , drop = FALSE]
      vxl_train_data$Y <- Y_vxl[valid_idx]

      # Subset warm start safely
      vxl_start <- NULL
      if (!is.null(warm_start)) {
        vxl_start <- warm_start_chunked[idvxl, ]
      }

      # Return only what the worker needs
      list(
        vxlcol         = vxl_reg_entry$voxel,
        data           = vxl_train_data,
        start          = vxl_start,
        registry_entry = vxl_reg_entry
      )
    }

    # Prepare input list for each voxel
    prepared_voxel_data <- future.apply::future_lapply(seq_along(ichunk), prep_func)

    # Clean up chunk matrices to free memory before parallel execution
    rm(voxeldata_chunked)
    rm(registry_chunked)
    if (!is.null(segmentation)) rm(voxelseg_chunked)
    if (!is.null(warm_start)) rm(warm_start_chunked)
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
                                    logfile        <- if (debug) file.path(logdir, paste0('log.vxl', vxlcol)) else NULL

                                    # GAMLSS fit using pre-assembled data
                                    g <- TRY(gamlss2::gamlss2(formula = g_form_parsed,
                                                              data    = vxl_item$data,
                                                              family  = g.family,
                                                              start   = vxl_item$start,
                                                              maxit   = maxit,
                                                              control =gamlss2::gamlss2_control(trace = FALSE,
                                                                                                light = TRUE,
                                                                                                 eps  = eps),
                                                              ...),
                                             logfile, save.env.and.stop = F)

                                    if (show_progress) { p() }

                                    # Error handling and deep environment stripping
                                    error = FALSE
                                    if (identical(g[1], NA)) {
                                      error <- TRUE
                                      g     <- list(vxl = vxlcol, error = TRUE, errout = g[2])

                                    } else {
                                      # Good fit, strip extras
                                      g$control      <- NULL
                                      g$converged    <- g$iterations < maxit[1L]
                                      g$family       <- g$family$family
                                      g$vxl          <- vxlcol
                                      g$call         <- NULL
                                      g$fake_formula <- NULL
                                    }

                                    if (cache){
                                      # Save on disk to unload master
                                      #  note use qs2 package instead of RDS for speed
                                      #saveRDS(g, registry_entry$full_paths, compress = FALSE)
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

                                    } else {
                                      # Just return the model
                                      return(g)
                                    }

                                    # End parallel loop
                                  }
    # End chunking routine
    gc()

    if (cache){
      # Update registry

      # Extract the updated values from the worker returns
      worker_returned_entry_voxels    <- vapply(submodels, function(x) x$voxel, numeric(1))
      worker_returned_entry_fitted    <- vapply(submodels, function(x) x$fitted, logical(1))
      worker_returned_entry_converged <- vapply(submodels, function(x) x$converged, logical(1))

      # Update the master registry in place
      match_idx <- match(worker_returned_entry_voxels, registry$voxel)
      registry$fitted[match_idx]    <- worker_returned_entry_fitted
      registry$converged[match_idx] <- worker_returned_entry_converged

      # Save the updated master registry
      qs2::qs_save(registry, registry_path)

    } else {
      # Assign directly to pre-allocated list slot
      models[[i]] <- submodels
    }

    gc()
    current.time <- Sys.time()
    }




  # ---------------------------------------------------------
  # CLOSING
  gc()

  if (! cache) {
    # Flatten the pre-allocated chunked list into a single model list
    models <- do.call(c, models)
    cat(paste0('Completed in ', format(current.time - start.time, digits = 3),"\n"))
    return(structure(models, class = "vbgamlss"))
  }

  # Collect and merge all model saved locally into one
  cat("Aggregating individual voxel models\n")
  models <- future.apply::future_lapply(registry$full_paths,
                                        function(x) {qs2::qs_read(x, nthreads = 1)})
  models <- structure(models, class = "vbgamlss")

  # Save
  if (! is.null(save_model)) {
    qs2::qs_save(models, paste0(save_model, '.vbgamlss'), nthreads = num_cores)
    cat('Model saved: ', paste0(save_model, '.vbgamlss'))
    }

  # Then remove the smaller temp RDS
  if (remove_voxel_fits_from_cache) {file.remove(registry$full_paths)}

  # Bye bye
  cat(paste0('Completed in ', format(current.time - start.time, digits = 3),"\n"))
  return(models)
}




# ================== #
#  TESTING VERSIONS  #
# ================== #


# ========== #
# DEPRECATED #
# ========== #


