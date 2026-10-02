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
                     cache_id = NULL, # give the cachedir a specific name ie {cachedir}/{<chache id>}.vbgamlss.cache.{<rand string>}
                     force_constraints=c(1e-8, +Inf),
                     warm_start=NULL, # per voxel by param
                     show_progress=T,
                     eps=1e-5,
                     maxit=c(100, 33),
                     future_plan_strategy = "future.mirai::mirai_cluster",
                     save_model=NULL,
                     deep_stripping = TRUE, # remove .Enviroments form fitted models
                     retry_failed = FALSE, # on resume, refit voxels whose status is error/nonconverged
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
  owned_cache <- FALSE # created (or found via cache_id) by vbgamlss, so safe to delete
  cache_root  <- cachedir

  # Resume: cachedir is either a cache itself, or holds a previous "<cache_id>.vbgamlss.cache.*"
  registry_path <- file.path(cachedir, '.vbgamlss.registry')
  if (!file.exists(registry_path) && !is.null(cache_id)) {
    prev <- list.dirs(cachedir, recursive = FALSE)
    prev <- prev[startsWith(basename(prev), paste0(cache_id, '.vbgamlss.cache.')) &
                   file.exists(file.path(prev, '.vbgamlss.registry'))]
    if (length(prev) > 0) {
      cachedir      <- prev[1]
      registry_path <- file.path(cachedir, '.vbgamlss.registry')
      owned_cache   <- TRUE
    }
  }

  registry <- NULL
  if (file.exists(registry_path)) {
    registry <- qs2::qs_read(registry_path)
    if (is.null(registry$shard)) {
      if (!owned_cache) { stop('Cache at ', cachedir, ' was written by an older vbgamlss (one file per voxel), delete it or pass another cachedir') }
      cat('Ignoring old-format cache', cachedir, '\n')
      registry <- NULL
      cachedir <- cache_root
    }
  }

  if (!is.null(registry)) {
    cat('Cache directory and registry found. Resuming from', cachedir, '\n')
    logdir      <- file.path(cachedir, '.voxlog')
    voxfits_dir <- file.path(cachedir, '.voxfits')
    registry$fitted <- registry$fitted & file.exists(file.path(cachedir, registry$shard))
    registry <- recover_shards(registry, cachedir) # voxels finished after the last registry save
    if (retry_failed) { registry$fitted[registry$status %in% c('error', 'nonconverged')] <- FALSE }
    report_registry(registry)

  } else {
    owned_cache <- TRUE

    # Make main cache folder
    fit_rand_id <- rand_names(1, l=4)
    cachedir <- file.path(cachedir, paste0(cache_id, '.vbgamlss.cache.', fit_rand_id))
    dir.create(cachedir, recursive = T, showWarnings = F)

    # Make folder for debugging
    logdir <- file.path(cachedir, '.voxlog')
    dir.create(logdir, recursive = T, showWarnings = F)

    # Shards: one append-only file per worker process, see append_shard()
    voxfits_dir <- file.path(cachedir, '.voxfits')
    dir.create(voxfits_dir, recursive = T, showWarnings = F)

    # Make registry. fitted = attempted; status: ok / warned / nonconverged / error / too_few_obs
    # shard/offset/nbytes: where the voxel model sits in its shard (shard relative to cachedir)
    registry <- data.frame(voxel     = 1:nvox,
                           fitted    = logical(nvox),
                           converged = logical(nvox),
                           status    = NA_character_,
                           message   = NA_character_,
                           shard     = NA_character_,
                           offset    = NA_real_,
                           nbytes    = NA_real_,
                           stringsAsFactors = FALSE)

    registry_path <- file.path(cachedir, '.vbgamlss.registry')
    qs2::qs_save(registry, registry_path)
    cat(paste0('Cache directory created: ', cachedir, '\n'))
  }

  # New shard files on every call: never append after a tail a killed job left half-written
  run_token <- paste0(Sys.getpid(), '.', as.integer(Sys.time()))



  # ---------------------------------------------------------
  # LARGE IMAGE CHUNKING & PREPPING ROUTINE

  # Parse formula ONCE. Its environment would otherwise be this whole frame (voxeldata included),
  # serialized to every worker together with the formula.
  g_form_parsed <- as.formula(g.formula)
  environment(g_form_parsed) <- globalenv()

  # Workers get the covariates once, and only those the formula uses
  tdata <- train.data[, intersect(all.vars(g_form_parsed), names(train.data)), drop = FALSE]
  rm(train.data)

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
    registry_rows <- registry_rows[is_unfitted]

    # If the entire chunk is already fitted, skip to the next chunk immediately
    if (length(ichunk) == 0) {
      cat("All voxels in this chunk are already fitted. Skipping...\n")
      next
    }

    # One small item per voxel, so each worker receives only its own voxels
    items <- lapply(seq_along(ichunk), function(j) {
      list(voxel = ichunk[j],
           y     = voxeldata[, ichunk[j]],
           seg   = if (!is.null(segmentation)) segmentation[, ichunk[j]])
    })


    # ---------------------------------------------------------
    # PARALLEL PROCESSING ROUTINE

    # Track progress per chunk
    if (show_progress) { p <- progressr::progressor(length(items)) }

    cat('Processing ', length(ichunk), ' voxels of', nvox,'\n')
    submodels <- foreach::foreach(vxl_item = items,
                                  .options.future = future.opt)  %dofuture% {

                                    # WORKER THREAD CONTROL
                                    if (RhpcBLASctl::blas_get_num_procs() > 1L)
                                          {RhpcBLASctl::blas_set_num_threads(1L)}
                                    if (RhpcBLASctl::omp_get_num_procs() > 1L)
                                          {RhpcBLASctl::omp_set_num_threads(1L)}

                                    vxlcol  <- vxl_item$voxel
                                    logfile <- NULL
                                    if (!is.null(logdir))
                                          {logfile <- file.path(logdir, paste0('log.vxl', vxlcol))}

                                    # Observations kept for this voxel
                                    Y_vxl <- as.numeric(vxl_item$y)
                                    valid <- !is.na(Y_vxl)
                                    if (!is.null(force_constraints)) {
                                      valid <- valid & Y_vxl >= force_constraints[1] & Y_vxl <= force_constraints[2]
                                    }
                                    if (!is.null(vxl_item$seg) && !is.null(segmentation_target)) {
                                      valid <- valid & (vxl_item$seg %in% segmentation_target)
                                    }

                                    # GAMLSS fit
                                    if (sum(valid) < 5) {
                                      fit    <- list(value = NA, warnings = character(0),
                                                     error = paste0('fewer than 5 valid observations (', sum(valid), ')'))
                                      status <- 'too_few_obs'
                                    } else {
                                      vxl_data   <- tdata[valid, , drop = FALSE]
                                      vxl_data$Y <- Y_vxl[valid]
                                      fit <- TRY(gamlss2::gamlss2(formula = g_form_parsed,
                                                                  data    = vxl_data,
                                                                  family  = g.family,
                                                                  start   = if (!is.null(warm_start)) warm_start[valid, ],
                                                                  maxit   = maxit,
                                                                  control = gamlss2::gamlss2_control(trace = FALSE,
                                                                                                     light = TRUE,
                                                                                                     eps  = eps),
                                                                  ...),
                                                 logfile)
                                      status <- if (!is.null(fit$error)) 'error'
                                                else if (fit$value$iterations >= maxit[1L]) 'nonconverged'
                                                else if (length(fit$warnings) > 0) 'warned'
                                                else 'ok'
                                    }
                                    msg <- if (!is.null(fit$error)) fit$error
                                           else if (length(fit$warnings) > 0) paste(fit$warnings, collapse = ' | ')
                                           else NA_character_

                                    if (show_progress) { p() }

                                    if (status %in% c('error', 'too_few_obs')) {
                                      g <- list(vxl = vxlcol, error = TRUE, converged = F, status = status, message = msg)
                                    } else {
                                      # Good fit, strip extras
                                      g <- fit$value
                                      # SEs of the linear coefficients for eb_shrink(): vcov() needs the data, only available here
                                      g$se_linear <- tryCatch(suppressWarnings(stats::vcov(g, type = "se")), error = function(e) NULL)
                                      g$control      <- NULL
                                      g$converged    <- g$iterations < maxit[1L]
                                      g$family       <- g$family$family
                                      g$vxl          <- vxlcol
                                      g$status       <- status
                                      g$message      <- msg
                                      if (deep_stripping){
                                        g <- deep_env_stripping(g)}
                                    }
                                    # Everything the fit had, so a fit that was not clean can be replayed alone, see refit_voxel()
                                    # ponytail: adds ~nsub x covariates doubles per such voxel; make optional if many warn
                                    if (status %in% c('warned', 'nonconverged', 'error')) {
                                      g$debug_call <- list(formula = g_form_parsed, data = vxl_data, family = g.family,
                                                           start = if (!is.null(warm_start)) warm_start[valid, ],
                                                           maxit = maxit,
                                                           control = list(trace = FALSE, light = TRUE, eps = eps),
                                                           dots = list(...))
                                    }

                                    # Append to this worker's shard to unload master
                                    shard <- paste0('shard.', run_token, '.', Sys.getpid(), '.bin')
                                    loc   <- append_shard(file.path(voxfits_dir, shard), vxlcol, qs2::qs_serialize(g))

                                    list(voxel = vxlcol, converged = isTRUE(g$converged),
                                         status = status, message = msg,
                                         shard = file.path('.voxfits', shard), offset = loc[1], nbytes = loc[2])
                                  }
    rm(items)
    gc()

    # Update the master registry in place
    match_idx <- match(vapply(submodels, function(x) x$voxel, numeric(1)), registry$voxel)
    registry$fitted[match_idx]    <- TRUE
    registry$converged[match_idx] <- vapply(submodels, function(x) x$converged, logical(1))
    registry$status[match_idx]    <- vapply(submodels, function(x) x$status, character(1))
    registry$message[match_idx]   <- vapply(submodels, function(x) x$message, character(1))
    registry$shard[match_idx]     <- vapply(submodels, function(x) x$shard, character(1))
    registry$offset[match_idx]    <- vapply(submodels, function(x) x$offset, numeric(1))
    registry$nbytes[match_idx]    <- vapply(submodels, function(x) x$nbytes, numeric(1))
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

  cat('\n'); report_registry(registry)

  # Aggregating
  cat("Aggregating individual voxel models\n")
  models <- vector("list", nvox)

  # raw bytes, read shard by shard in file order
  done <- which(registry$fitted)
  for (sh in unique(registry$shard[done])) {
    rows <- done[registry$shard[done] == sh]
    rows <- rows[order(registry$offset[rows])]
    con  <- file(file.path(cachedir, sh), "rb")
    for (r in rows) {
      seek(con, registry$offset[r] + SHARD_HEADER_BYTES)
      models[[r]] <- readBin(con, what = "raw", n = registry$nbytes[r])
    }
    close(con)
    gc(verbose = FALSE)
  }
  models <- structure(models, class = "vbgamlss", cachedir = cachedir)


  # Save
  if (! is.null(save_model)) {
    qs2::qs_save(models,
                 file = paste0(save_model, '.vbgamlss'),
                 compress_level = 0L) # uncompressed
    cat('Model saved: ', paste0(save_model, '.vbgamlss'), '\n')
  }


  # Cleanup
  if (!debug) {
    if (owned_cache &&                                                          # made or found by vbgamlss, not passed in
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



#' Define the S3 subsetting method
#'
#' @rawNamespace S3method("[[", vbgamlss)
`[[.vbgamlss` <- function(x, i, ...) {
  raw_bytes <- unclass(x)[[i, ...]]
  if (is.null(raw_bytes)) return(NULL)
  restore_family(qs2::qs_deserialize(raw_bytes))
}



# Deep environment stripping
deep_env_stripping <- function(model) {
  target_env <- environment(model$terms$mu)
  # explode passed enviroment, too many links to it to find them all manually
  if (!is.null(target_env) && !identical(target_env, globalenv()) && !identical(target_env, baseenv())) {
    rm(list = ls(envir = target_env, all.names = TRUE), envir = target_env)
  }
  # gamlss2 also pins the frame it was called from (the worker's) on the terms list
  if (!is.null(attr(model$terms, ".Environment"))) {
    attr(model$terms, ".Environment") <- globalenv()
  }

  return(model)
}



# Per-status voxel counts of a vbgamlss registry
report_registry <- function(registry) {
  st <- ifelse(registry$fitted, registry$status, 'pending')
  st <- table(factor(st, levels = c('ok', 'warned', 'nonconverged', 'error', 'too_few_obs', 'pending')))
  cat('Of', nrow(registry), 'voxels:', paste(names(st), st, sep = '=', collapse = ', '), '\n')
  errs <- registry$message[registry$fitted & registry$status %in% c('error', 'too_few_obs')]
  if (length(errs) > 0) {
    top <- sort(table(errs), decreasing = TRUE)[1:min(3, length(unique(errs)))]
    cat('\t most frequent failures:\n', paste0('\t   ', names(top), ' (', top, ')\n'), sep = '')
  }
}



# Shard record: [magic | voxel | nbytes] as 4-byte integers, then nbytes of qs2-serialized model.
# Each worker process appends to its own file, so a file never has two writers.
SHARD_MAGIC <- 1447249713L # "VBG1"
SHARD_HEADER_BYTES <- 12

# Append one voxel model; returns c(offset of the record, nbytes of the model)
append_shard <- function(path, voxel, bytes) {
  offset <- if (file.exists(path)) file.size(path) else 0
  con <- file(path, "ab")
  on.exit(close(con))
  writeBin(c(SHARD_MAGIC, as.integer(voxel), length(bytes)), con, size = 4)
  writeBin(bytes, con)
  c(offset, length(bytes))
}

# Index of the complete records in a shard; stops at a truncated or corrupt tail
scan_shard <- function(path) {
  size <- file.size(path)
  con  <- file(path, "rb")
  on.exit(close(con))
  vox <- integer(0); off <- numeric(0); nb <- numeric(0); pos <- 0
  while (pos + SHARD_HEADER_BYTES <= size) {
    h <- readBin(con, "integer", n = 3, size = 4)
    if (length(h) < 3 || h[1] != SHARD_MAGIC || pos + SHARD_HEADER_BYTES + h[3] > size) break
    vox <- c(vox, h[2]); off <- c(off, pos); nb <- c(nb, h[3])
    pos <- pos + SHARD_HEADER_BYTES + h[3]
    seek(con, pos)
  }
  data.frame(voxel = vox, offset = off, nbytes = nb)
}

# Register voxels that are complete in the shards but missing from the registry
# (fitted after its last save, e.g. the job was killed mid-chunk)
recover_shards <- function(registry, cachedir) {
  shards <- list.files(file.path(cachedir, '.voxfits'), pattern = '^shard\\..*\\.bin$')
  n_rec <- 0
  for (sh in shards) {
    idx <- scan_shard(file.path(cachedir, '.voxfits', sh))
    idx <- idx[idx$voxel %in% registry$voxel[!registry$fitted], , drop = FALSE]
    if (nrow(idx) == 0) next
    con <- file(file.path(cachedir, '.voxfits', sh), "rb")
    for (k in seq_len(nrow(idx))) {
      seek(con, idx$offset[k] + SHARD_HEADER_BYTES)
      g <- qs2::qs_deserialize(readBin(con, what = "raw", n = idx$nbytes[k]))
      r <- match(idx$voxel[k], registry$voxel)
      registry[r, c('shard', 'offset', 'nbytes')] <- list(file.path('.voxfits', sh), idx$offset[k], idx$nbytes[k])
      registry$fitted[r]    <- TRUE
      registry$converged[r] <- isTRUE(g$converged)
      registry$status[r]    <- if (is.null(g$status)) NA_character_ else g$status
      registry$message[r]   <- if (is.null(g$message)) NA_character_ else g$message
      n_rec <- n_rec + 1
    }
    close(con)
  }
  if (n_rec > 0) { cat('Recovered', n_rec, 'voxels from shards that the registry had not recorded yet\n') }
  registry
}



#' Replay the fit of one voxel alone, e.g. to debug a fit that failed, warned or did not converge
#'
#' Every voxel whose status is warned/nonconverged/error stores `$debug_call`: the exact data,
#' formula, family, start, maxit, control and extra arguments its gamlss2 fit received.
#' @param models A vbgamlss model.
#' @param i Voxel index.
#' @param ... Override any stored argument, e.g. `control = list(trace = TRUE)` (merged into the
#'   stored control), `maxit = c(300, 50)`, `family = SHASH`.
#' @return The gamlss2 fit (errors are not caught, so traceback() works).
#' @export
refit_voxel <- function(models, i, ...) {
  a <- models[[i]]$debug_call
  if (is.null(a)) { stop('voxel ', i, ' (status: ', models[[i]]$status, ') has no debug_call, only warned/nonconverged/error fits keep it') }
  over <- list(...)
  if (!is.null(over$control)) { a$control <- utils::modifyList(a$control, over$control); over$control <- NULL }
  args <- c(list(formula = a$formula, data = a$data, family = a$family, start = a$start, maxit = a$maxit,
                 control = do.call(gamlss2::gamlss2_control, a$control)), a$dots)
  do.call(gamlss2::gamlss2, utils::modifyList(args, over, keep.null = TRUE))
}



# Delete a model's voxel cache once the aggregated model is safely on disk
drop_vbgamlss_cache <- function(model) {
  cdir <- attr(model, "cachedir")
  if (!is.null(cdir) && grepl("\\.vbgamlss\\.cache", basename(cdir))) { unlink(cdir, recursive = TRUE) }
}



# ================== #
#  TESTING VERSIONS  #
# ================== #


# # CONSIDER REPLACING list of binary with binary of binary stored locally
#
# cat("Consolidating models into a single binary archive\n")
#
# archive_file <- "vbgamlss_archive.bin"
# con <- file(archive_file, "wb") # Open for binary writing
#
# offsets <- numeric(nvox)
# lengths <- numeric(nvox)
# current_byte_pos <- 0
#
# for (i in seq_along(registry$full_paths)) {
#   # 1. Read the temporary individual file (or generate the raw bytes directly)
#   f_size <- file.info(registry$full_paths[i])$size
#   raw_bytes <- readBin(registry$full_paths[i], what = "raw", n = f_size)
#
#   # 2. Record where this model lives in the giant file
#   offsets[i] <- current_byte_pos
#   lengths[i] <- length(raw_bytes)
#
#   # 3. Write bytes to the archive and advance the position counter
#   writeBin(raw_bytes, con)
#   current_byte_pos <- current_byte_pos + lengths[i]
# }
#
# close(con)


# # Your object now uses almost zero RAM
# models <- structure(
#   list(
#     archive = archive_file,
#     offsets = offsets,
#     lengths = lengths
#   ),
#   class = "vbgamlss"
# )
#

# # S3 method class to load it
# `[[.vbgamlss` <- function(x, i, ...) {
#   # Open connection to the single large file
#   con <- file(x$archive, "rb")
#
#   # Ensure the connection closes even if an error occurs
#   on.exit(close(con))
#
#   # Jump straight to the model's location on the disk
#   seek(con, where = x$offsets[i], origin = "start")
#
#   # Read only the bytes for this specific model
#   raw_bytes <- readBin(con, what = "raw", n = x$lengths[i])
#
#   restore_family(qs2::qs_deserialize(raw_bytes))
# }














# ========== #
# DEPRECATED #
# ========== #










