#############################  model selection  ################################
# this set of functions are hard coded to work with our local SLURM HPC systems







#' @export
vbgamlss.model_selection_LOBO <- function(# model selection commands
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
                    out <- vbgamlss.cv( g.formula    = g.formula,
                                        imageframe   = imageframe,
                                        train.data   = dtfr,
                                        fold.var     = fold.var,
                                        g.family     = g.family,
                                        segmentation = segmentation,
                                        chunk_max_mb = chunk_max_mb,
                                        force_constraints = g.constraints,
                                        k.penalty    = k.penalty,
                                        verbose      = verbose,
                                        return_all_GD = return_all_GD,
                                        num_cores    = NULL,
                                        debug        = T,
                                        save_states  = T,
                                        resume       = T,
                                        drop_re      = T,
                                        logdir     = slurm$wd)
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





# -----------------------------------------------------------
slurm_template <- function(){
  return('#!/bin/bash
#SBATCH --job-name=<%= slurm$jobname %>
#SBATCH --output=<%= slurm$output %>
#SBATCH --error=<%= slurm$error %>
#SBATCH --time=<%= slurm$time %>
#SBATCH --mem-per-cpu=<%= slurm$mempercpu %>
#SBATCH --cpus-per-task=<%= slurm$ncpu %>
#SBATCH --nodes=1

# Enviroment variables
module load singularityce/4.1.0
export SINGULARITY_BINDPATH="/users,/tmp,/reference,/dcsrsoft,/work,/scratch,/data"
containeR=/data/PRTNR/CHUV/RADMED/ijelescu/micmap/tom/norming/containers/NormMod_R_v1.3.sif
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=<%= slurm$ncpu %>

# call
singularity exec $containeR Rscript <%= SCRIPT %>
  ')}






# -----------------------------------------------------------
slurm_resources <- function(n='VBGAMLSS', o=NULL, e=NULL,
                            t='71:59:00', m='40G', c='12', a=1, r=NULL) {
  slurm           <- list()
  slurm$jobname   <- as.character(n)     # name of the slurm job
  slurm$output    <- as.character(o)     # path to the .out slurm file
  slurm$error     <- as.character(e)     # path to the .err slurm file
  slurm$time      <- as.character(t)     # max slurm time
  slurm$mempercpu <- as.character(m)     # memory per cpu
  slurm$ncpu      <- as.character(c)     # number of cpus
  slurm$rdsout    <- as.character(r)     # path to the results qs file (kept var name to avoid upstream breaks)
  slurm$array     <- a                   # array jobs
  return(slurm)
}






# -----------------------------------------------------------
slurm_registry <- function(Njobs, env=NULL, run_name="cv_run"){
  registry              <- list()
  registry$wd           <- getwd()
  registry$main         <- file.path(registry$wd, "vbgamlss.slurm")

  registry$path         <- file.path(registry$main,
                                     paste0('.', run_name, '_', rand_names(1), '_', Sys.Date()))
  registry$reg_path     <- file.path(registry$path, '.Slurm.Registry')

  registry$jobs         <- paste0('.', run_name, '_job', seq(Njobs),'_', rand_names(Njobs))

  registry$jobspaths    <- file.path(registry$path, registry$jobs)
  registry$jobs_sout    <- file.path(registry$jobspaths,
                                     paste0(registry$jobs, '.stdout'))
  registry$jobs_results <- file.path(registry$jobspaths,
                                     paste0(registry$jobs, '.results.qs'))
  registry$jobs_script  <- file.path(registry$jobspaths,
                                     paste0(registry$jobs, '.R'))
  registry$jobs_sbatch  <- file.path(registry$jobspaths,
                                     paste0(registry$jobs, '.job'))
  registry$jobs_env     <- file.path(registry$jobspaths,
                                     paste0(registry$jobs, '.local.env'))
  quite(lapply(registry$jobspaths, dir.create, recursive = TRUE))
  registry$env <- NULL
  return(registry)
}









# -----------------------------------------------------------
sanity_check <- function(file_list){
  ## safety check ##
  for (file_path in file_list) {
    # Check if the file exists
    if (!file.exists(file_path)) {
      # Stop execution and return an error message if the file is not found
      stop("ERROR: can't find file but it should have been generated:", file_path)
    }
  }
}







# -----------------------------------------------------------
sbatch_jobs <- function(registry) {
  # example: 'Submitted batch job 22119441'
  jobs_id <- c()
  i <- 0
  for (job in registry$jobs_sbatch){

    # sbatch
    job_status <- system(glue("sbatch {job}"), intern = TRUE)
    i <- i + 1
    cat(glue('{i}  '), fill=T)

    # sanity stop
    if (! startsWith(job_status, 'Submitted batch job')) {
      stop(paste0("SBATCHING went wrong, sbatch call returned: ", job_status))}

    # parse
    id <- as.numeric(sub("^Submitted batch job ", "", job_status))
    jobs_id <- c(jobs_id, id)
    Sys.sleep(0.5)
  }
  registry$jobs_id <- jobs_id
  registry$killall <- paste('scancel', paste(registry$jobs_id, collapse=" "))
  return(registry)
}






# -----------------------------------------------------------
# To be changed, SLURM has output in machine readable format
jobs_status <- function(registry, param='status') {
  # sacct -j 22119441 --format="JobID,State,Elapsed,ExitCode,End"
  # JobID             State    Elapsed ExitCode                 End
  # ------------ ---------- ---------- -------- -------------------
  #   22119441         FAILED   00:00:02      1:0 2024-07-12T11:54:51

  jobs_status <- c()
  jobs_elapsed <- c()
  jobs_exitcode <- c()
  jobs_end <- c()
  for (job in registry$jobs_id){
    status_ <- system(glue('sacct -j {job} --format="JobID,State,ExitCode,Elapsed,End"'),
                      intern = TRUE)
    jobs_status <- c(jobs_status, strsplit(status_[3], '\\s+')[[1]][2])
    jobs_exitcode <- c(jobs_exitcode, strsplit(status_[3], '\\s+')[[1]][3])
    jobs_elapsed <- c(jobs_elapsed, strsplit(status_[3], '\\s+')[[1]][4])
    jobs_end <- c(jobs_end, strsplit(status_[3], '\\s+')[[1]][5])
  }
  if (param=='status'){return(jobs_status)}
  if (param=='elapsed'){return(jobs_elapsed)}
  if (param=='exitcode'){return(jobs_exitcode)}
  if (param=='end'){return(jobs_end)}
}






# -----------------------------------------------------------
monitor_jobs <- function(registry, sleep=10, resbatch=NULL) {
  RUNNING = TRUE
  start_time = Sys.time()
  finished <- F
  ct=0
  while (RUNNING) {
    Sys.sleep(1)

    if ((ct %% sleep) == 0) {
      ct <- 0
      # Format status into a data frame
      registry$status <- jobs_status(registry, param='status')
      status_df <- data.frame(JobID = registry$jobs_id,
                              JobCode = registry$jobs,
                              Status = registry$status,
                              Time = jobs_status(registry, param='elapsed'),
                              ExitCode = jobs_status(registry, param='exitcode'),
                              End = jobs_status(registry, param='end')
      )
    }

    # Text printing
    system("clear")
    cat(' ', fill=T)
    cat(paste('Sbatched time        ', format(start_time)), fill=T)
    cat(' ', fill=T)
    cat(paste('Last squeue refresh  ', format(Sys.time())), fill=T)
    cat(' ', fill=T)
    cat(paste('Elapsed              ', format(Sys.time()-start_time, digits=3)), fill=T)
    cat(' ', fill=T)
    print(status_df)
    cat(' ', fill=T)
    cat(paste('Refresh in ', sleep - ct - 1, ' sec'), fill=T)
    if (! is.null(resbatch)) {cat(paste('Resbatch ', resbatch), fill=T)}


    # Check if all jobs are terminated
    if ((ct %% sleep) == 0) {
      if (all(registry$status %in% c("COMPLETED", "FAILED",
                                     "CANCELLED", "BOOT_FAIL",
                                     "DEADLINE", "OUT_OF_MEMORY",
                                     "OUT_OF_ME+",
                                     "NODE_FAIL", "PREEMPTED",
                                     "REVOKED", "TIMEOUT"))) {
        finished <- T
        break
      } }

    # arbitrary max time limit (include pending jobs)
    if (difftime(Sys.time(), start_time, units='hours') > 24*15) {
      # kill and resbatch
      system(registry$killall)
      break
    }

    # Time counter increment
    ct <- ct + 1
  }

  if (! finished) {system("clear")} # skip clear if it is done
  return(registry)
}



# print(status_df)
# cat(paste0('\nAll Jobs terminated:', Sys.time(), ''), fill=T)
# cat(paste0('Registry path:', registry$path), fill=T)
# cat("Done", fill = T)




# -----------------------------------------------------------
gather_jobs_outputs <- function(registry){
  final <- setNames(lapply(seq_along(registry$formula),
                           function(i) qs2::qs_read(registry$jobs_results[[i]])),
                    registry$formula)
  return(final)
}






# -----------------------------------------------------------
#' @export
cleanup_cv_cached_voxels <- function(resume_registry, dry_run = TRUE) {

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

    # Format formula nicely for the log
    f_str <- "Unknown Formula"
    if (!is.null(slurm_reg$formula) && length(slurm_reg$formula) >= i) {
      f_str <- paste(deparse(slurm_reg$formula[[i]]), collapse = " ")
    }

    cat(sprintf("Job %d/%d;", i, n_jobs), f_str, "\n")

    state_dir <- file.path(job_dir, "vbgamlss.cv.states")

    if (!dir.exists(state_dir)) {
      cat("\t| No CV states directory found, skipping. \n")
      next
    }

    # Locate all random cache directories for this specific job's folds
    cache_dirs <- list.dirs(state_dir, recursive = FALSE, full.names = TRUE)
    cache_dirs <- cache_dirs[grepl("\\.vbgamlss\\.cache", basename(cache_dirs))]

    if (length(cache_dirs) == 0) {
      cat("\t| No vbgamlss caches found in states directory, skipping.\n")
      next
    }

    # Iterate through the caches generated for this job's folds
    for (j in seq_along(cache_dirs)) {
      cdir <- cache_dirs[j]
      cache_name <- basename(cdir)
      local_reg_path <- file.path(cdir, ".vbgamlss.registry")
      voxfits_dir <- file.path(cdir, ".voxfits")

      if (!file.exists(local_reg_path)) {
        cat(sprintf("\t\t| Cache %d/%d (%s): Missing internal registry, skipping.\n", j, length(cache_dirs), cache_name))
        next
      }

      # Safely read the local voxel registry
      local_reg <- tryCatch({ qs2::qs_read(local_reg_path) }, error = function(e) NULL)

      if (is.null(local_reg)) {
        cat(sprintf("\t\t| Cache %d/%d (%s): Corrupted internal registry, skipping.\n", j, length(cache_dirs), cache_name))
        next
      }

      # Check convergence status
      total_voxels     <- nrow(local_reg)
      fitted_voxels    <- sum(local_reg$fitted == TRUE, na.rm = TRUE)
      converged_voxels <- sum(local_reg$converged == TRUE, na.rm = TRUE)

      if (fitted_voxels < total_voxels) {
        cat(sprintf("\t\t| Cache %d/%d (%s): Incomplete (%d/%d fitted), skipping.\n",
                    j, length(cache_dirs), cache_name, fitted_voxels, total_voxels))
        next
      }

      if (converged_voxels < total_voxels) {
        cat(sprintf("\t\t| Cache %d/%d (%s): %d voxels failed to converge, skipping.\n",
                    j, length(cache_dirs), cache_name, total_voxels - converged_voxels))
        next
      }

      # If we reach here, this fold's cache is fully fitted and 100% converged. Nuke the temporary files.
      if (dir.exists(voxfits_dir)) {
        files_to_delete <- list.files(voxfits_dir, full.names = TRUE, all.files = TRUE)
        n_files <- length(files_to_delete)

        if (n_files > 0) {
          if (dry_run) {
            cat(sprintf("\t\t| Cache %d/%d (%s): dry_run=TRUE, %d files could be deleted.\n",
                        j, length(cache_dirs), cache_name, n_files))
          } else {
            cat(sprintf("\t\t| Cache %d/%d (%s): Deleting %d files.\n",
                        j, length(cache_dirs), cache_name, n_files))
            unlink(files_to_delete, force = TRUE)
          }
          total_freed <- total_freed + n_files
        } else {
          cat(sprintf("\t\t| Cache %d/%d (%s): .voxfits already empty, skipping.\n", j, length(cache_dirs), cache_name))
        }
      } else {
        cat(sprintf("\t\t| Cache %d/%d (%s): No .voxfits directory, skipping.\n", j, length(cache_dirs), cache_name))
      }
    }
    cat("\n")
  }

  cat("Done.\n")
  return(NULL)
}

























