








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
