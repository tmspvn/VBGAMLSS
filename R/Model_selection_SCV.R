








#' @export
vbgamlss.model_selection_NCV <- function(# model selection commands
  result_file,  # where to save
  run_name = "ncv_run",
  resume_registry=NULL,  # pass registry file, works as toggle for resume
  reset_enviroment=T,   # reset passed global environment instead of resuming
  sbatch_resources=list(t='71:59:00', m='40G', c='12'),
  # gamlss nested cv commands
  mu_formulas,
  sigma_formulas,
  nu_formulas,
  tau_formulas,
  fold.var, # outer LOBO variable (e.g. batch/site), see LOSOfolds()
  k_inner = 5, # number of inner stratified CV folds run within each outer (LOBO) training set
  drop_re = T, # drop random/batch effects when predicting the held-out outer batch
  images,  # pass named list of paths
  constraints, # pass named list of numeric vectors
  mask,
  train.data,
  families = c('NO'),
  segmentation = NULL,
  chunk_max_mb = 256,
  k.penalty=NULL,
  verbose=F,
  return_all_GD=T,
  ...){

  if (!is.character(train.data)) { stop("train.data class must be a path") }

  # ---------------------------------------------------------
  # PREPARE INPUTS

  # Combine formulas for gamlss2
  all_formulas <- combine_formulas_gamlss2(mu_formulas, sigma_formulas,
                                           nu_formulas, tau_formulas)

  # Extract inputs safely
  img_names <- names(images)
  images_vec <- unlist(images, use.names = FALSE)
  families_vec <- unlist(families, use.names = FALSE)

  # Match constraints to images by name (safest), or fallback to positional mapping
  if (!is.null(img_names) && !is.null(names(constraints))) {
    matched_constraints <- constraints[img_names]
  } else if (length(constraints) == 1 && length(images_vec) > 1) {
    matched_constraints <- rep(constraints, length(images_vec))
  } else {
    matched_constraints <- constraints
  }

  # Generate all combinations using an index to avoid messing up the list
  grid_combinations <- expand.grid(
    img_idx = seq_along(images_vec),
    formula = all_formulas,
    family  = families_vec,
    stringsAsFactors = FALSE
  )

  # Build the final structured grid
  job_grid <- data.frame(
    image      = images_vec[grid_combinations$img_idx],
    formula    = grid_combinations$formula,
    family     = grid_combinations$family,
    stringsAsFactors = FALSE
  )
  # Append the constraint as a list-column to preserve the numeric vectors
  job_grid$constraint <- matched_constraints[grid_combinations$img_idx]
  nfm <- nrow(job_grid)


  # ---------------------------------------------------------
  # PREPARE REGISTRY

  ### Resume registry? ###
  if (is.null(resume_registry)) {
    cat(paste0('Generating *new* registry\n'))
    # Prepare registry
    registry <- slurm_registry(nfm, run_name = run_name)

    # Merge grid columns directly into the registry list
    registry$image      <- job_grid$image
    registry$formula    <- job_grid$formula
    registry$family     <- job_grid$family
    registry$constraint <- job_grid$constraint
    registry$mask       <- mask

    registry$pkgs <- loadedNamespaces()

    # Save session info
    writeLines(capture.output(sessionInfo()),
               file(paste0(registry$path, '/.SessionInfo.txt')))

    # Save the environment & registry to file
    registry$genv <- paste0(registry$path, '/.Global.Enviroment.RData')
    save.image(file = registry$genv)

    # Save registry using qs2
    qs2::qs_save(registry, registry$reg_path)

  } else {
    cat(paste0('Loading *past* registry\n'))
    registry <- qs2::qs_read(resume_registry)
    # reset global environment, useful in case of changes to core functions
    if (reset_enviroment) {
      cat(paste0('Reset global enviroment\n'))
      save.image(file = registry$genv)
    }
  }


  # ---------------------------------------------------------
  # SLURM JOBS PREPARATION

  ### Prepare template ###
  template <- slurm_template()

  ### Generate jobs to the HPC ###
  if (is.null(resume_registry)) {
    cat(paste0('Generating ', nfm, ' jobs\n'))
    for (i in 1:nfm) {
      # assign correct paths
      slurm <- slurm_resources(n=registry$jobs[i],
                               o=registry$jobs_sout[i],
                               e=registry$jobs_sout[i],
                               r=registry$jobs_results[i],
                               t=sbatch_resources$t,
                               m=sbatch_resources$m,
                               c=sbatch_resources$c)
      slurm$wd <- registry$jobspaths[i]
      slurm$jenv <- registry$jobs_env[i]

      # Extract cleanly directly from the registry list
      image         <- registry$image[i]
      g.formula     <- registry$formula[i]
      g.family      <- registry$family[i]
      g.mask        <- registry$mask

      # Use [[i]] to extract the actual numeric vector from the list element
      g.constraints <- registry$constraint[[i]]

      # save job local environment
      save(list = ls(all.names = TRUE), file = slurm$jenv)

      # Populate call
      CALL <- glue("
                   load('{registry$genv}')
                   load('{slurm$jenv}')
                   suppressWarnings(
                     suppressPackageStartupMessages(
                       invisible(
                         lapply(registry$pkgs, require, character.only = TRUE, quietly = TRUE))))
                   set.seed(04281945)
                   imageframe <- images2matrix(image, g.mask)
                   dtfr <- read.csv(train.data, stringsAsFactors = T)
                   out <- vbgamlss.nested_cv( g.formula         = g.formula,
                                              imageframe        = imageframe,
                                              train.data        = dtfr,
                                              fold.var          = fold.var,
                                              k_inner           = k_inner,
                                              drop_re           = drop_re,
                                              g.family          = g.family,
                                              segmentation      = segmentation,
                                              chunk_max_mb      = chunk_max_mb,
                                              force_constraints = g.constraints,
                                              k.penalty         = k.penalty,
                                              verbose           = verbose,
                                              return_all_GD     = return_all_GD,
                                              num_cores         = NULL,
                                              debug             = T,
                                              save_states       = T,
                                              resume            = T,
                                              logdir            = slurm$wd)
                   warnings()
                   qs2::qs_save(out, slurm$rdsout)
                   ")

      SCRIPT <- registry$jobs_script[i]
      JOB <- registry$jobs_sbatch[i]

      # Write call to file
      writeLines(CALL, SCRIPT)
      # Populate template
      brew(text = template, output = JOB)
    }
  } else {
    cat(paste0('Using the pre-generated ', nfm, ' jobs\n'))
  }

  sanity_check(c(registry$jobs_script, registry$jobs_sbatch))


  # ---------------------------------------------------------
  # SLURM JOBS SBATCHING and monitor jobs to HPC

  #resbatch (max 4 times) if timeout or other
  RESBATCH = 0
  user_interrupted <- FALSE
  while (TRUE){
    ## Sbatch jobs to the HPC ##
    cat(paste0('Sbatching ', nfm, ' jobs\n'))
    registry <- sbatch_jobs(registry)

    ## MONITOR ##
    cat(paste0('Monitoring ', nfm, ' jobs\n'))
    tryCatch({
      registry <- monitor_jobs(registry, sleep=5, resbatch=RESBATCH)
    },
    interrupt = function(e) {user_interrupted <<- TRUE})

    if (user_interrupted) {
      system(registry$killall)
      stop("Script interrupted by the user!\nKilling all the sbatched jobs...\n")
    }

    # Check if all jobs are terminated
    if (all(registry$status %in% c("FAILED", "COMPLETED"))) {break}

    # max cluster time limit
    if (RESBATCH > 4) {break}
    RESBATCH <- RESBATCH + 1
  }


  # ---------------------------------------------------------
  # Gather the results
  cat(paste0('gathering results\n'))
  results <- gather_jobs_outputs(registry)

  # Extra, disposable: a rank-based composite ranking across candidate formulas,
  # combining the outer LOBO (between-batch) and inner CV (within-distribution)
  # metrics stored in `results`. Stored alongside, not in place of, the raw
  # per-formula results -- rank_composite_model_selection() can be re-run later
  # with different weights/metrics without recomputing anything.
  model_ranking <- tryCatch(rank_composite_model_selection(results),
                            error = function(e) {
                              warning("rank_composite_model_selection() failed: ", conditionMessage(e))
                              NULL
                            })

  qs2::qs_save(list(results = results, model_ranking = model_ranking), result_file)
  cat(paste0('Done.\n\n\n\n'))
  warnings()
}




# -----------------------------------------------------------
# Delete fully-converged .voxfits caches from a "vbgamlss.cv.states"-style
# directory (shared by vbgamlss.lobocv, vbgamlss.stratified_cv and each inner
# fold of vbgamlss.nested_cv). Returns the number of files (dry-run: that could
# be) deleted.
.cleanup_vbgamlss_cache_dir <- function(state_dir, label = "", dry_run = TRUE) {
  total_freed <- 0
  if (!dir.exists(state_dir)) { return(total_freed) }

  cache_dirs <- list.dirs(state_dir, recursive = FALSE, full.names = TRUE)
  cache_dirs <- cache_dirs[grepl("\\.vbgamlss\\.cache", basename(cache_dirs))]
  if (length(cache_dirs) == 0) { return(total_freed) }

  for (j in seq_along(cache_dirs)) {
    cdir <- cache_dirs[j]
    cache_name <- basename(cdir)
    local_reg_path <- file.path(cdir, ".vbgamlss.registry")
    voxfits_dir <- file.path(cdir, ".voxfits")

    if (!file.exists(local_reg_path)) {
      cat(sprintf("\t\t| %s cache %d/%d (%s): Missing internal registry, skipping.\n",
                  label, j, length(cache_dirs), cache_name))
      next
    }

    local_reg <- tryCatch({ qs2::qs_read(local_reg_path) }, error = function(e) NULL)
    if (is.null(local_reg)) {
      cat(sprintf("\t\t| %s cache %d/%d (%s): Corrupted internal registry, skipping.\n",
                  label, j, length(cache_dirs), cache_name))
      next
    }

    total_voxels     <- nrow(local_reg)
    fitted_voxels    <- sum(local_reg$fitted == TRUE, na.rm = TRUE)
    converged_voxels <- sum(local_reg$converged == TRUE, na.rm = TRUE)

    if (fitted_voxels < total_voxels) {
      cat(sprintf("\t\t| %s cache %d/%d (%s): Incomplete (%d/%d fitted), skipping.\n",
                  label, j, length(cache_dirs), cache_name, fitted_voxels, total_voxels))
      next
    }
    if (converged_voxels < total_voxels) {
      cat(sprintf("\t\t| %s cache %d/%d (%s): %d voxels failed to converge, skipping.\n",
                  label, j, length(cache_dirs), cache_name, total_voxels - converged_voxels))
      next
    }

    if (dir.exists(voxfits_dir)) {
      files_to_delete <- list.files(voxfits_dir, full.names = TRUE, all.files = TRUE)
      n_files <- length(files_to_delete)
      if (n_files > 0) {
        if (dry_run) {
          cat(sprintf("\t\t| %s cache %d/%d (%s): dry_run=TRUE, %d files could be deleted.\n",
                      label, j, length(cache_dirs), cache_name, n_files))
        } else {
          cat(sprintf("\t\t| %s cache %d/%d (%s): Deleting %d files.\n",
                      label, j, length(cache_dirs), cache_name, n_files))
          unlink(files_to_delete, force = TRUE)
        }
        total_freed <- total_freed + n_files
      }
    }
  }
  return(total_freed)
}




#' @export
cleanup_ncv_cached_voxels <- function(resume_registry, dry_run = TRUE) {

  if (!file.exists(resume_registry)) {
    stop("Error: SLURM registry file not found at ", resume_registry)
  }

  cat("Loading SLURM registry:", resume_registry, "\n")
  slurm_reg <- qs2::qs_read(resume_registry)
  n_jobs <- length(slurm_reg$jobspaths)
  cat("Found", n_jobs, "jobs (formulas) in the registry.\n\n")

  total_freed <- 0

  for (i in seq_len(n_jobs)) {
    job_dir <- slurm_reg$jobspaths[i]

    f_str <- "Unknown Formula"
    if (!is.null(slurm_reg$formula) && length(slurm_reg$formula) >= i) {
      f_str <- paste(deparse(slurm_reg$formula[[i]]), collapse = " ")
    }
    cat(sprintf("Job %d/%d;", i, n_jobs), f_str, "\n")

    ncv_state_dir <- file.path(job_dir, "vbgamlss.ncv.states")
    if (!dir.exists(ncv_state_dir)) {
      cat("\t| No nested-CV states directory found, skipping.\n\n")
      next
    }

    # Outer-level caches (shared across outer folds, cachedir = ncv_state_dir)
    total_freed <- total_freed + .cleanup_vbgamlss_cache_dir(ncv_state_dir, label = "outer", dry_run = dry_run)

    # Inner-level caches, one "vbgamlss.cv.states" dir per outer fold
    outer_dirs <- list.dirs(ncv_state_dir, recursive = FALSE, full.names = TRUE)
    outer_dirs <- outer_dirs[grepl("^outer_fold_", basename(outer_dirs))]
    for (od in outer_dirs) {
      inner_state_dir <- file.path(od, "inner_cv", "vbgamlss.cv.states")
      total_freed <- total_freed + .cleanup_vbgamlss_cache_dir(inner_state_dir,
                                                               label = paste0(basename(od), "/inner"),
                                                               dry_run = dry_run)
    }
    cat("\n")
  }

  cat("Done.\n")
  return(NULL)
}
