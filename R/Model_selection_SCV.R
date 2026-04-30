








#' @export
vbgamlss.model_selection_SCV <- function(# model selection commands
  result_file,  # where to save
  run_name = "cv_run",
  resume_registry=NULL,  # pass registry file, works as toggle for resume
  reset_enviroment=T,   # reset passed global environment instead of resuming
  sbatch_resources=list(t='71:59:00', m='40G', c='12'),
  # gamlss cv commands
  mu_formulas,
  sigma_formulas,
  nu_formulas,
  tau_formulas,
  fold.var,
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
                   out <- vbgamlss.stratified_cv( g.formula         = g.formula,
                                                  imageframe        = imageframe,
                                                  train.data        = dtfr,
                                                  fold.var          = fold.var,
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
  qs2::qs_save(results, result_file)
  cat(paste0('Done.\n\n\n\n'))
  warnings()
}
