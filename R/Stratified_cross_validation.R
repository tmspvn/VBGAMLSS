








#' @export
vbgamlss.stratified_cv <- function(imageframe,
                                   g.formula,
                                   train.data,
                                   fold.var, # vector of length Nsubjects with integers indicating fold di appartenenza
                                   g.family = NO,
                                   segmentation = NULL,
                                   segmentation_target=NULL,
                                   num_cores = NULL,
                                   chunk_max_mb = 64,
                                   force_constraints=c(1e-8, +Inf),
                                   k.penalty=NULL,
                                   verbose=F,
                                   debug=T,
                                   logdir=getwd(),
                                   return_all_GD=T,
                                   save_states=T,
                                   resume=T,
                                   ...) {

  cat(paste0("Starting cross validation"), fill=T)

  if (missing(imageframe)) { stop("imageframe is missing") }
  if (missing(fold.var)) { stop("vector of length N subjects with integers indicating the fold is missing") }
  if (missing(g.formula)) { stop("formula is missing")}
  check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}

  # get length fold var
  n_folds <- length(unique(fold.var))

  # states
  if (save_states) {
    state.dir=file.path(logdir, paste0('vbgamlss.cv.states'))
    dir.create(state.dir, recursive = T, showWarnings = F)
    reg.file=file.path(state.dir, paste0('.registry.qs'))
    # Load or
    if (!file.exists(reg.file)) {
      registry = matrix(data = F, nrow=4, ncol=n_folds)
      colnames(registry) <- paste0('fold', 1:n_folds)
      rownames(registry) <- c('modfit', 'GD', 'stat', 'done')
      if (save_states) {
        save.image(file.path(state.dir, '.env.Rdata'))
      }
    } else {
      registry = qs2::qs_read(reg.file)
    }
  }

  # check registry and verify from where to skip
  fold_completed <- apply(registry, 2, all)
  resume_from <- 1
  for (i in 1:n_folds) {
    if (!fold_completed[i]) {
      resume_from = i
      break # first non completed fold
    }
  }

  # Force character columns to factors
  train.data <- as.data.frame(unclass(train.data), stringsAsFactors=TRUE)

  # Create or load stratified folds
  fold.file = file.path(state.dir, '.folds.variable.qs')
  if (resume & file.exists(fold.file)){
    cat(paste0('\t| Found fold information, loading it'), fill=T)
    train.data$folds <- qs2::qs_read(fold.file)
  } else {
    train.data$folds <- fold.var
    if (save_states){qs2::qs_save(train.data$folds, fold.file)}
  }

  # Prepare storage for results
  cvresults <- c()

  # Perform stratified cross-validation
  for (fold in resume_from:n_folds) {
    cat(paste0("\nProcessing fold ", fold, " of ", n_folds), fill=T)

    # Split data into training and validation sets
    # train
    training_fold <- train.data$folds != fold
    train_fold_data <- train.data[training_fold, ]
    train_fold_data <- droplevels(train_fold_data) # drop levels
    # test
    test_indices <- which(train.data$folds == fold)
    test_fold_data <- train.data[test_indices, ]

    # Fit or resume the model on the training fold
    fold.model.file = file.path(state.dir, paste0('.fold.', fold, '.model.qs'))
    if (resume & file.exists(fold.model.file)){
      cat(paste0("\t| Training: found existing train fold models, loading them"), fill=T)
      model <- qs2::qs_read(fold.model.file)
    } else {
      cat('\t| Training: fitting train fold', fill=T)
      cat("\033[34m")
      model <- quite(
        vbgamlss(imageframe=imageframe[training_fold,],
                 g.formula=g.formula, # parsing handled in Core.R
                 train.data=train_fold_data,
                 g.family=g.family,
                 segmentation=segmentation,
                 segmentation_target=segmentation_target,
                 force_constraints = force_constraints,
                 num_cores=num_cores,
                 chunk_max_mb=chunk_max_mb,
                 debug=debug,
                 cachedir=state.dir,
                 ...
        ),
        skip=verbose)
      cat("\033[0m")
      if (save_states){qs2::qs_save(model, fold.model.file)}

      # update and save status on the registry
      registry[1, fold] = TRUE
      if (save_states){qs2::qs_save(registry, reg.file)}
    }

    # Predict new fold
    fold.gd.file = file.path(state.dir, paste0('.fold.', fold, '.GD.qs'))
    if (resume & file.exists(fold.gd.file)){
      cat(paste0("\t| Test: found existing fold GDs, loading them"), fill=T)
      GDs <- qs2::qs_read(fold.gd.file)
    } else {
      cat('\t| Test: predict test metrics', fill=T)
      GDs <- predictGD.strat(model,
                       test_imageframe     = imageframe[test_indices,],
                       newdata             = test_fold_data,
                       verbose             = verbose,
                       segmentation        = segmentation,
                       segmentation_target = segmentation_target,
                       resume              = resume,
                       save_states         = save_states,
                       loginfo             = c(fold, state.dir))

      if (save_states){qs2::qs_save(GDs, fold.gd.file)}

      # update and save status on the registry
      registry[2, fold] = TRUE
      if (save_states){qs2::qs_save(registry, reg.file)}
    }

    # get statistics for validation Global Deviance brain-wide
    cat('\t| Summarizing statistics', fill=T)
    stats <- statGD(GDs,
                    k.penalty,
                    deg.fre=model[[1]]$df,
                    return_all_GD=return_all_GD)

    # update and save status on the registry
    registry[3, fold] = TRUE

    # Store the model and validation set
    cvresults[[fold]] <- stats
    fold.resuts.file = file.path(state.dir, paste0('.fold.', fold, '.results.qs'))
    if (! file.exists(fold.resuts.file)){
      if (save_states){qs2::qs_save(cvresults, fold.resuts.file)}

      # update and save status on the registry
      registry[4, fold] = TRUE
      if (save_states){qs2::qs_save(registry, reg.file)}
    }

    # clean
    rm(model)
    rm(GDs)
    gc()
  }


  # --------------------------------------------------------
  # Gather previously completed folds from disk
  if (resume && resume_from > 1) {
    cat('\t| Compiling final results: loading previously completed folds', fill = T)
    for (past_fold in 1:(resume_from - 1)) {

      # Check the registry matrix
      if (registry[4, past_fold] == TRUE) {
        past_res_file <- file.path(state.dir, paste0('.fold.', past_fold, '.results.qs'))

        if (file.exists(past_res_file)) {
          past_data <- tryCatch(qs2::qs_read(past_res_file), error = function(e) NULL)

          # Safely extract the target fold from the past list
          if (!is.null(past_data) && !is.null(past_data[[past_fold]])) {
            cvresults[[past_fold]] <- past_data[[past_fold]]
          } else {
            warning(paste("Could not extract data for fold", past_fold, "from the saved state file."))
          }
        }
      }
    }
  }

  # --------------------------------------------------------
  # Give a name to each fold
  names(cvresults) <- paste0("fold", seq_along(cvresults))

  # save the final result
  if (save_states){
    cvres_filepath <- file.path(state.dir, 'vbgamlss.cvresults')
    qs2::qs_save(cvresults, file = cvres_filepath)
    cat('Saved CV results: ', cvres_filepath, '\n')
  }

  return(cvresults)
}



# --------------------------------
# Predict new fold Global Deviance
predictGD.strat <- function (object,
                       test_imageframe,
                       newdata = NULL,
                       verbose=F,
                       segmentation=NULL,
                       segmentation_target=NULL,
                       loginfo=c(fold, state.dir),
                       resume=T,
                       save_states=T,
                       ...) {

  if (is.null(newdata)){stop("newdata is not set")}
  .fold = loginfo[1]
  .state.dir = loginfo[2]

  ## predict GD ##
  familyobj <- restore_family(object[[1]])$family

  # .predicted.parameters
  fold.P.file = file.path(.state.dir, paste0('.fold.', .fold, '.predicted.parameters.qs'))
  if (resume & file.exists(fold.P.file)){
    cat(paste0("\t\t| found existing test fold predicted parameters, loading them"), fill=T)
    nfitted <- qs2::qs_read(fold.P.file)
  } else {
    quite(cat('\t| predicting test fold parameters', fill=T), skip=verbose)

    # Just predict, do not drop re
    cat("\033[34m")
    nfitted <- quite(
      predict.vbgamlss(object,
                       newdata             = newdata,
                       ptype               = 'parameter',
                       segmentation        = segmentation,
                       segmentation_target = segmentation_target
      ),
      skip=verbose)
    cat("\033[0m")

    # save
    if (save_states){qs2::qs_save(nfitted, fold.P.file)}
  }

  # .predicted.response
  fold.R.file = file.path(.state.dir, paste0('.fold.', .fold, '.predicted.response.qs'))
  if (resume & file.exists(fold.R.file)){
    cat(paste0("\t| Found existing fold predicted response, loading them"), fill=T)
    resp <- qs2::qs_read(fold.R.file)
  } else {
    quite(cat('\t| Predicting fold response', fill=T), skip=verbose)

    cat("\033[34m")
    resp <- quite(
      predict.vbgamlss(object,
                       newdata             = newdata,
                       ptype               = 'response',
                       segmentation        = segmentation,
                       segmentation_target = segmentation_target
      ),
      skip=verbose)
    cat("\033[0m")
    if (save_states){qs2::qs_save(resp, fold.R.file)}
  }

  # add response to nfitted (from gamlss2.predict)
  nsub <- dim(newdata)[1] # == length(nfitted[[1]]$mu)
  for(i in 1:length(nfitted)){
    condA = ! "try-error" %in% class(nfitted[[i]])
    condB = ! all(is.na(nfitted[[i]]))
    if (condA & condB) {
      nfitted[[i]][['y']] <- test_imageframe[,i] # subj x vox mat
      nfitted[[i]][['yhat']] <- as.numeric(unlist(resp[[i]])[1:nsub])
    } else {
      nfitted[[i]] <- NA
    }
  }
  rm(resp)
  rm(object)

  # missfits
  not_missfits <- ! is.na(nfitted)

  # test GD
  quite(cat('\t| Evaluating test fold GD ', fill=T), skip=verbose)
  plan(strategy="future::cluster", workers=availableCores())
  GDs <- foreach(i=seq_along(nfitted)) %dofuture% {
    if (not_missfits[i]) {

      # Wrapped in tryCatch to catch uniroot failures
      tryCatch({
        vxlGD <- testGD(nfitted[[i]], familyobj)
        vxlGD$vxl <- nfitted[[i]]$vxl
        vxlGD
      }, error = function(e) {
        print(e)
        NA # Return NA if the math fails
      })

    } else {
      NA
    }
  }

  class(GDs) <- "vbgamlss.predictions.GD"
  GDs
}




#############################  nested cross validation  ##########################
# Outer loop = leave-one-batch-out (between-batch generalization)
# Inner loop = stratified k-fold CV run on each outer training set
#              (within-distribution predictive power, uncontaminated by the
#              held-out batch)
#
#' @export
vbgamlss.nested_cv <- function(imageframe,
                               g.formula,
                               train.data,
                               fold.var, # outer LOBO variable: integer/factor per subject indicating batch
                               g.family = NO,
                               segmentation = NULL,
                               segmentation_target = NULL,
                               num_cores = NULL,
                               chunk_max_mb = 64,
                               force_constraints = c(1e-8, +Inf),
                               k.penalty = NULL,
                               k_inner = 5, # number of inner stratified CV folds run within each outer training set
                               verbose = F,
                               debug = T,
                               logdir = getwd(),
                               return_all_GD = T,
                               save_states = T,
                               resume = T,
                               drop_re = T, # drop random/batch effects when predicting the held-out outer batch
                               ...) {

  cat(paste0("Starting nested cross validation (outer LOBO x inner ", k_inner, "-fold CV)"), fill=T)

  if (missing(imageframe)) { stop("imageframe is missing") }
  if (missing(fold.var)) { stop("vector of length N subjects with integers indicating the outer (batch) fold is missing") }
  if (missing(g.formula)) { stop("formula is missing") }
  check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing") }

  n_outer <- length(unique(fold.var))

  # states
  state.dir <- file.path(logdir, 'vbgamlss.ncv.states')
  dir.create(state.dir, recursive = T, showWarnings = F)
  reg.file <- file.path(state.dir, '.registry.qs')
  if (!file.exists(reg.file)) {
    registry <- matrix(data = F, nrow = 4, ncol = n_outer)
    colnames(registry) <- paste0('outerfold', 1:n_outer)
    rownames(registry) <- c('outer_modfit', 'outer_pred', 'outer_stat', 'inner_cv')
    if (save_states) { save.image(file.path(state.dir, '.env.Rdata')) }
  } else {
    registry <- qs2::qs_read(reg.file)
  }

  fold_completed <- apply(registry, 2, all)
  resume_from <- 1
  for (i in 1:n_outer) {
    if (!fold_completed[i]) {
      resume_from <- i
      break
    }
  }

  # Force character columns to factors
  train.data <- as.data.frame(unclass(train.data), stringsAsFactors = TRUE)

  # Create or load the outer (LOBO) fold variable
  fold.file <- file.path(state.dir, '.outer.folds.variable.qs')
  if (resume & file.exists(fold.file)) {
    cat(paste0('\t| Found outer fold information, loading it'), fill=T)
    train.data$outerfolds <- qs2::qs_read(fold.file)
  } else {
    train.data$outerfolds <- fold.var
    if (save_states) { qs2::qs_save(train.data$outerfolds, fold.file) }
  }

  ncvresults <- list(outer = list(), inner = list())

  # ---------------------------------------------------------
  # Outer loop: leave-one-batch-out
  for (fold in resume_from:n_outer) {
    cat(paste0("\nProcessing OUTER fold ", fold, " of ", n_outer), fill=T)

    outer.dir <- file.path(state.dir, paste0('outer_fold_', fold))
    dir.create(outer.dir, recursive = T, showWarnings = F)

    training_fold <- train.data$outerfolds != fold
    train_fold_data <- train.data[training_fold, ]
    train_fold_data <- droplevels(train_fold_data)
    test_indices <- which(train.data$outerfolds == fold)
    test_fold_data <- train.data[test_indices, ]

    # ---- Outer model fit (trained on all-but-one batch) ----
    outer.model.file <- file.path(outer.dir, '.outer.model.qs')
    if (resume & file.exists(outer.model.file)) {
      cat(paste0("\t| Outer training: found existing model, loading it"), fill=T)
      outer.model <- qs2::qs_read(outer.model.file)
    } else {
      cat('\t| Outer training: fitting on all-but-one-batch', fill=T)
      cat("\033[34m")
      outer.model <- quite(
        vbgamlss(imageframe = imageframe[training_fold, ],
                 g.formula = g.formula,
                 train.data = train_fold_data,
                 g.family = g.family,
                 segmentation = segmentation,
                 segmentation_target = segmentation_target,
                 force_constraints = force_constraints,
                 num_cores = num_cores,
                 chunk_max_mb = chunk_max_mb,
                 debug = debug,
                 cachedir = state.dir,
                 ...
        ),
        skip = verbose)
      cat("\033[0m")
      if (save_states) { qs2::qs_save(outer.model, outer.model.file) }
      registry[1, fold] <- TRUE
      if (save_states) { qs2::qs_save(registry, reg.file) }
    }

    # ---- Outer LOBO prediction on the held-out batch (unseen level -> drop_re) ----
    outer.gd.file <- file.path(outer.dir, '.outer.GD.qs')
    if (resume & file.exists(outer.gd.file)) {
      cat(paste0("\t| Outer test: found existing GDs, loading them"), fill=T)
      outer.GDs <- qs2::qs_read(outer.gd.file)
    } else {
      cat('\t| Outer test: predicting held-out batch (LOBO)', fill=T)
      outer.GDs <- predictGD(outer.model,
                             test_imageframe     = imageframe[test_indices, ],
                             newdata             = test_fold_data,
                             verbose             = verbose,
                             segmentation        = segmentation,
                             segmentation_target = segmentation_target,
                             resume              = resume,
                             save_states         = save_states,
                             drop_re             = drop_re,
                             loginfo             = c('outer', outer.dir))
      if (save_states) { qs2::qs_save(outer.GDs, outer.gd.file) }
      registry[2, fold] <- TRUE
      if (save_states) { qs2::qs_save(registry, reg.file) }
    }

    cat('\t| Summarizing outer (LOBO) statistics', fill=T)
    outer_stats <- statGD(outer.GDs,
                          k.penalty,
                          deg.fre = outer.model[[1]]$df,
                          return_all_GD = return_all_GD)
    ncvresults$outer[[paste0('outerfold', fold)]] <- outer_stats
    registry[3, fold] <- TRUE
    if (save_states) { qs2::qs_save(registry, reg.file) }

    rm(outer.model, outer.GDs)
    gc()

    # ---- Inner loop: stratified k-fold CV within this outer training set ----
    inner.dir <- file.path(outer.dir, 'inner_cv')
    dir.create(inner.dir, recursive = T, showWarnings = F)

    # Stratify inner folds by the (remaining) outer batch variable, so each
    # inner fold keeps the batch composition of the outer training set balanced.
    inner_fold_var <- stratCVfolds(as.factor(train.data$outerfolds[training_fold]), k.fold = k_inner)

    cat('\t| Inner CV: running ', k_inner, '-fold stratified CV on outer training set', fill=T)
    inner_results <- vbgamlss.stratified_cv(imageframe          = imageframe[training_fold, ],
                                            g.formula          = g.formula,
                                            train.data         = train_fold_data,
                                            fold.var           = inner_fold_var,
                                            g.family           = g.family,
                                            segmentation       = segmentation,
                                            segmentation_target = segmentation_target,
                                            num_cores          = num_cores,
                                            chunk_max_mb       = chunk_max_mb,
                                            force_constraints  = force_constraints,
                                            k.penalty          = k.penalty,
                                            verbose            = verbose,
                                            debug              = debug,
                                            logdir             = inner.dir,
                                            return_all_GD      = return_all_GD,
                                            save_states        = save_states,
                                            resume             = resume,
                                            ...)
    ncvresults$inner[[paste0('outerfold', fold)]] <- inner_results
    registry[4, fold] <- TRUE
    if (save_states) { qs2::qs_save(registry, reg.file) }

    gc()
  }

  # --------------------------------------------------------
  # Gather previously completed outer folds from disk (resume)
  if (resume && resume_from > 1) {
    cat('\t| Compiling final results: loading previously completed outer folds', fill = T)
    for (past_fold in 1:(resume_from - 1)) {
      if (registry[4, past_fold] == TRUE) {
        outer.dir <- file.path(state.dir, paste0('outer_fold_', past_fold))
        outer.gd.file <- file.path(outer.dir, '.outer.GD.qs')
        outer.model.file <- file.path(outer.dir, '.outer.model.qs')

        if (file.exists(outer.gd.file)) {
          past_GDs <- tryCatch(qs2::qs_read(outer.gd.file), error = function(e) NULL)
          deg.fre <- NA
          if (file.exists(outer.model.file)) {
            m <- tryCatch(qs2::qs_read(outer.model.file), error = function(e) NULL)
            if (!is.null(m)) { deg.fre <- m[[1]]$df; rm(m) }
          }
          if (!is.null(past_GDs)) {
            ncvresults$outer[[paste0('outerfold', past_fold)]] <-
              statGD(past_GDs, k.penalty, deg.fre = deg.fre, return_all_GD = return_all_GD)
          }
        }

        inner.dir <- file.path(outer.dir, 'inner_cv')
        inner.cvres.file <- file.path(inner.dir, 'vbgamlss.cv.states', 'vbgamlss.cvresults')
        if (file.exists(inner.cvres.file)) {
          past_inner <- tryCatch(qs2::qs_read(inner.cvres.file), error = function(e) NULL)
          if (!is.null(past_inner)) {
            ncvresults$inner[[paste0('outerfold', past_fold)]] <- past_inner
          }
        }
      }
    }
  }

  # --------------------------------------------------------
  # Aggregate summary statistics
  cat('\t| Aggregating nested CV summary', fill=T)
  outer_summary <- getCV_All_Metrics(ncvresults$outer)
  inner_flat <- unlist(ncvresults$inner, recursive = FALSE)
  inner_summary <- getCV_All_Metrics(inner_flat)

  ncvresults$summary <- list(LOBO = outer_summary, innerCV = inner_summary)
  ncvresults$n_outer <- n_outer
  ncvresults$k_inner <- k_inner

  if (save_states) {
    ncvres_filepath <- file.path(state.dir, 'vbgamlss.ncvresults')
    qs2::qs_save(ncvresults, file = ncvres_filepath)
    cat('Saved nested CV results: ', ncvres_filepath, '\n')
  }

  return(ncvresults)
}




# --------------------------------
# Rank-based composite model selection across candidate formulas.
#
# This is a *post-hoc* helper, kept deliberately separate from the CV engine
# above: it only reads the metrics already stored by vbgamlss.nested_cv, so the
# selection rule can be swapped/re-run at any time without recomputing any fit.
#
# `results` is the named list produced by vbgamlss.model_selection_NCV / gathered
# via gather_jobs_outputs(): one vbgamlss.nested_cv() output per candidate formula.
#' @export
rank_composite_model_selection <- function(results,
                                           metrics = c('GD', 'MAE', 'LL', 'CLL'),
                                           stat = 'mean',
                                           lobo_weight = 0.5,
                                           inner_weight = 0.5) {

  formulas <- names(results)
  nfm <- length(formulas)

  get_val <- function(res, source, metric) {
    v <- tryCatch(res$summary[[source]][[metric]][[stat]], error = function(e) NA)
    if (is.null(v)) NA else v
  }

  lobo_tab <- sapply(metrics, function(m) sapply(results, get_val, source = 'LOBO', metric = m))
  inner_tab <- sapply(metrics, function(m) sapply(results, get_val, source = 'innerCV', metric = m))
  colnames(lobo_tab) <- paste0('LOBO_', metrics)
  colnames(inner_tab) <- paste0('innerCV_', metrics)

  # lower is better for all of GD/MAE/LL/CLL -> rank ascending
  rank_lobo <- apply(lobo_tab, 2, rank, na.last = 'keep', ties.method = 'average')
  rank_inner <- apply(inner_tab, 2, rank, na.last = 'keep', ties.method = 'average')
  colnames(rank_lobo) <- paste0('rank_LOBO_', metrics)
  colnames(rank_inner) <- paste0('rank_innerCV_', metrics)

  mean_rank_lobo <- rowMeans(rank_lobo, na.rm = TRUE)
  mean_rank_inner <- rowMeans(rank_inner, na.rm = TRUE)
  composite_rank <- (lobo_weight * mean_rank_lobo) + (inner_weight * mean_rank_inner)

  out <- data.frame(formula = formulas,
                    lobo_tab, inner_tab,
                    rank_lobo, rank_inner,
                    mean_rank_LOBO = mean_rank_lobo,
                    mean_rank_innerCV = mean_rank_inner,
                    composite_rank = composite_rank,
                    stringsAsFactors = FALSE)

  out <- out[order(out$composite_rank), ]
  rownames(out) <- NULL
  return(out)
}
