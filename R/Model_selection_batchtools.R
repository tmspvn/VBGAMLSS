#############################  model selection  ################################








vbgamlss.model_selection <- function(RESULTS_FILE,
                                      SLURM_JOB_TMPL,
                                      RESOURCES,
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
                                     DEBUG_HPC=F,
                                     DELAYS=NULL,
                                      ...){

  if (is.null(DELAYS)) {DELAYS=c(30, 2)}

  # Combine formulas for gamlss2:
  all_formulas <- combine_formulas_gamlss2(mu_formulas, sigma_formulas,
                                           nu_formulas, tau_formulas,
                                           families)

  # Setup the future plan to use the batchtools back-end with the Slurm template
  if (DEBUG_HPC) {options(future.debug = TRUE, future.delete=FALSE)}
  plan(batchtools_slurm,
       template = SLURM_JOB_TMPL,
       resources = RESOURCES,
       scheduler.latency=DELAYS[1],
       fs.latency=DELAYS[2])

  # Submit jobs to the HPC
  cat(paste0('Sbatching ', length(all_formulas), ' jobs'), fill = TRUE)
  results <- foreach(f = all_formulas) %dofuture% {
                      # parse formula
                       print(f)
                      parsed <- strsplit(f, " :: ", fixed = TRUE)
                      g.family <- parsed[[1]][1]
                      g.formula <- parsed[[1]][2]
                      # fit
                      future::future({
                        vbgamlss.cv(g.formula = g.formula,
                                    image = image,
                                    mask = mask,
                                    train.data = train.data,
                                    fold.var = fold.var,
                                    g.family = g.family,
                                    segmentation = segmentation,
                                    chunk_max_mb = chunk_max_mb,
                                    n_folds = n_folds,
                                    k.penalty = k.penalty,
                                    verbose = verbose,
                                    return_all_GD = return_all_GD,
                                    subsample = subsample,
                                    subsample.type = subsample.type,
                                    num_cores = NULL)
                        })
  }

  # Gather the results
  results <- lapply(results, value)
  saveRDS(results, RESULTS_FILE)
}



