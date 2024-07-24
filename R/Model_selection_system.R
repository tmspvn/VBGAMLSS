#############################  model selection  ################################








vbgamlss.model_selection <- function(result_file,
                                      mu_formulas,
                                      sigma_formulas,
                                      nu_formulas,
                                      tau_formulas,
                                      fold.var,
                                      image,
                                      mask,
                                      train.data,
                                      families =  c('NO'),
                                      segmentation = NULL,
                                      chunk_max_mb = 128,
                                      n_folds = 10,
                                      k.penalty=NULL,
                                      verbose=F,
                                      return_all_GD=T,
                                      subsample=NULL,
                                      subsample.type=c('regular', 'random'),
                                      ...){

  if (class(train.data) != "character") { stop("train.data class must be a path")}

  # Combine formulas for gamlss2:
  all_formulas <- combine_formulas_gamlss2(mu_formulas, sigma_formulas,
                                           nu_formulas, tau_formulas,
                                           families)
  nfm <- length(all_formulas)

  # Prepare registry
  template <- slurm_template()
  registry <- slurm_registry(nfm)
  registry$formulas <- all_formulas
  registry$pkgs <- loadedNamespaces()
  # Save the environment to the file
  registry$env <- paste0(registry$path, '/.Enviroment.RData')
  save.image(file = registry$env)
  # Save session info
  writeLines(capture.output(sessionInfo()),
             file(paste0(registry$path, '/.SessionInfo.txt')))

  #### Generate jobs to the HPC ####
  cat(paste0('Generating ', nfm, ' jobs'), fill = TRUE)

  for (i in 1:nfm) {
    # assign correct paths
    slurm <- slurm_resources(n=registry$jobs[i],
                             o=registry$jobs_sout[i],
                             e=registry$jobs_sout[i],
                             r=registry$jobs_results[i])

    # Parse formula for slurm call
    parsed <- strsplit(registry$formulas[i], " :: ", fixed = TRUE)
    g.family <- parsed[[1]][1]
    g.formula <- parsed[[1]][2]
    # Populate call
    CALL <- glue("
    load('{registry['env']}')
    quite(lapply({registry['pkgs']}, require, character.only = TRUE))
    out <- vbgamlss.cv(g.formula = '{g.formula}',
                      image = '{image}',
                      mask = '{mask}',
                      train.data = read.csv('{train.data}'),
                      fold.var = {fold.var},
                      g.family = '{g.family}',
                      segmentation = '{segmentation}',
                      chunk_max_mb = {chunk_max_mb},
                      n_folds = {n_folds},
                      k.penalty = {if (is.null(k.penalty)) 'NULL' else k.penalty},
                      verbose = {verbose},
                      return_all_GD = {return_all_GD},
                      subsample = {if (is.null(subsample)) 'NULL' else subsample},
                      subsample.type = '{subsample.type}',
                      num_cores = NULL)
    saveRDS(out, '{slurm$rdsout}')
    ")

    SCRIPT <- registry$jobs_script[i]
    JOB <- registry$jobs_sbatch[i]
    # Write call to fine
    writeLines(CALL, SCRIPT)
    # Populate template
    brew(text = template,
         output = JOB)
  }

  sanity_check(c(registry$jobs_script, registry$jobs_sbatch))


  #### Sbatch jobs to the HPC ####
  cat(paste0('Sbatching ', nfm, ' jobs'), fill = TRUE)
  registry$jobs_id <- sbatch_jobs(registry)

  ### MONITOR ###
  cat(paste0('Monitoring ', nfm, ' jobs'), fill = TRUE)
  tryCatch({
    monitor_jobs(registry, sleep=5)
  }, interrupt = function(e) {
    cat("Script interrupted by the user!", fill=T)
    cat("Killing the sbatched jobs...", fill=T)
    registry$killall <- paste('scancel', paste(registry$jobs_id, collapse=" "))
    system(registry$killall)
  })

  ### Gather the results ###
  results <- gather_jobs_outputs(registry) # to finish
  #saveRDS(results, RESULTS_FILE)
}


slurm_template <- function(){
return(
'#!/bin/bash
#SBATCH --job-name=<%= slurm$jobname %>
#SBATCH --output=<%= slurm$output %>
#SBATCH --error=<%= slurm$error %>
#SBATCH --time=<%= slurm$time %>
#SBATCH --mem-per-cpu=<%= slurm$mempercpu %>
#SBATCH --cpus-per-task=<%= slurm$ncpu %>
#SBATCH --nodes=1

# Enviroment variables
module load singularityce/3.11.3
containeR=/scratch/tpavan1/scripts_tommaso/test_batchtools/NormMod_R.sif
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=<%= slurm$ncpu %>

# call
singularity exec $containeR Rscript <%= SCRIPT %>

' # end template
)}


slurm_resources <- function(n='VBGAMLSS', o=NULL, e=NULL,
                            t='71:59:00',m='30G', c='12', a=1, r=NULL) {
  slurm <- list()
  slurm$jobname <- as.character(n)
  slurm$output <- as.character(o)
  slurm$error <- as.character(e)
  slurm$time <- as.character(t)
  slurm$mempercpu <- as.character(m)
  slurm$ncpu <- as.character(c)
  slurm$rdsout <- as.character(r)
  slurm$array <- a
  return(slurm)
}

slurm_registry <- function(Njobs, env){
  registry <- list()
  registry$wd <- getwd()
  registry$main <- file.path(registry$wd, ".vbgamlss.slurm")
  registry$path <-file.path(registry$main, paste0('.', rand_names(1), '_', Sys.Date()))
  registry$jobs <- paste0('.', rand_names(Njobs),'_', seq(Njobs))
  registry$jobspaths <- file.path(registry$path, registry$jobs)
  registry$jobs_sout <- file.path(registry$jobspaths,
                                  paste0(registry$jobs, '.stdout'))
  registry$jobs_results <- file.path(registry$jobspaths,
                                     paste0(registry$jobs, '.RDS'))
  registry$jobs_script <- file.path(registry$jobspaths,
                                    paste0(registry$jobs, '.R'))
  registry$jobs_sbatch <- file.path(registry$jobspaths,
                                     paste0(registry$jobs, '.job'))
  quite(lapply(registry$jobspaths, dir.create, recursive = TRUE))
  registry$env <- NULL
  return(registry)
}


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
  return(jobs_id)
}

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

monitor_jobs <- function(registry, sleep=10) {
  RUNNING = TRUE
  while (RUNNING) {
    Sys.sleep(sleep)
    # Format status into a data frame
    registry$status <- jobs_status(registry, param='status')
    status_df <- data.frame(JobID = registry$jobs_id,
                            JobCode = registry$jobs,
                            Status = registry$status,
                            Time = jobs_status(registry, param='elapsed'),
                            ExitCode = jobs_status(registry, param='exitcode'),
                            End = jobs_status(registry, param='end')
    )
    system("clear")
    cat(' ', fill=T)
    cat(paste('Last squeue refresh', Sys.time()), fill=T)
    cat(' ', fill=T)
    print(status_df)

    # Check if all jobs are terminated
    if (all(registry$status %in% c("COMPLETED", "FAILED",
                                 "CANCELLED", "BOOT_FAIL",
                                 "DEADLINE", "OUT_OF_MEMORY",
                                 "NODE_FAIL", "PREEMPTED",
                                 "REVOKED", "TIMEOUT"))) {
      break
    }
  }
  system("clear")
  print(status_df)
  cat(paste0('\nAll Jobs terminated:', Sys.time(), ''), fill=T)
  cat(paste0('Registry path:', registry$path), fill=T)
  cat("Done", fill = T)
}

gather_jobs_outputs <- function(registry){
  NULL
}


