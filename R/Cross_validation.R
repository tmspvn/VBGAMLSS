#############################  cross validation  ###############################
#
#
#
#
#
#
#
#
#
#
#



#' @export
vbgamlss.cv <- function(image,
                        mask,
                        g.formula,
                        train.data,
                        fold.var, # vector of length Nsubjects with integers indicating fold di appartenenza
                        g.family = NO,
                        segmentation = NULL,
                        segmentation_target=NULL,
                        num_cores = NULL,
                        chunk_max_mb = 64,
                        k.penalty=NULL,
                        verbose=F,
                        debug=T,
                        logdir=getwd(),
                        return_all_GD=T,
                        save_states=T,
                        resume=T,
                        ...) {

  cat(paste0("Starting"), fill=T)

  if (! is.null(segmentation)){
    if (segmentation == 'NULL') {segmentation=NULL}}
  if (missing(image)) { stop("image is missing") }
  if (missing(mask)) { stop("mask is missing") }
  if (missing(fold.var)) { stop("vector of legnth Nsubjects with integers indicating the fold is missing") }
  if (missing(g.formula)) { stop("formula is missing")}
  check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}

  # get length fold var
  n_folds <- length(unique(fold.var))

  # states
  if (save_states) {
    state.dir=file.path(logdir, paste0('.cvstates'))
    dir.create(state.dir, recursive = T, showWarnings =F)
    reg.file=file.path(state.dir, paste0('.registry'))
    if (!dir.exists(reg.file)) {
      registry = matrix(data = F, nrow=4, ncol=n_folds)
      colnames(registry) <- paste0('fold', 1:n_folds)
      rownames(registry) <- c('modfit', 'GD', 'stat', 'done')
      if (save_states) {
        save.image(file.path(state.dir, '.env.Rdata'))
      }
    } else {
      registry = readRDS(reg.file)
      }
  }
  # check registry and verify from where to skip
  fold_completed <- apply(registry, 2, all)
  for (i in 1:n_folds) {
    if (!fold_completed[i]) {
      resume_from = i
      break # first non completed fold
      }
  }


  # Force character columns to factors
  train.data <- as.data.frame(unclass(train.data), stringsAsFactors=TRUE)

  # Create or load stratified folds
  fold.file = file.path(state.dir, '.fold.rds')
  if (resume & file.exists(fold.file)){
    cat(paste0("\nFound fold information, loading it"), fill=T)
    train.data$folds <- readRDS(fold.file)
  } else {
    train.data$folds <- fold.var
    if (save_states){saveRDS(train.data$folds, fold.file)}
  }

  # Prepare storage for results
  cvresults <- c()

  # Perform stratified cross-validation
  for (fold in resume_from:n_folds) {
    cat(paste0("\n\n Processing fold ", fold, " of ", n_folds), fill=T)


    # Split data into training and validation sets
    training_fold = train.data$folds != fold
    test_indices <- which(train.data$folds == fold)
    test_fold_data <- train.data[test_indices, ]
    cat('-Fitting fold \n', fill=T)


    # Fit or resume the model on the training fold
    fold.model.file = file.path(state.dir, paste0('.fold.', fold, '.model'))
    if (resume & file.exists(fold.model.file)){
      cat(paste0("Found fold models, loading them"), fill=T)
      model <- readRDS(fold.model.file)
    } else {
      model <- quite(
                      vbgamlss(image=image,
                               mask=mask,
                               g.formula=as.formula(g.formula),
                               train.data=train.data,
                               g.family=g.family,
                               segmentation=segmentation,
                               segmentation_target=segmentation_target,
                               num_cores=num_cores,
                               chunk_max_mb=chunk_max_mb,
                               afold=training_fold, # pass fold indexes
                               debug=T,
                               logdir=logdir
                               ),
                     skip=verbose)
      if (save_states){saveRDS(model, fold.model.file)}
      # update and save status on the registry
      registry[1, fold] = TRUE
      if (save_states){saveRDS(registry, reg.file)}
    }


    # Predict new fold
    cat('-Estimating GD \n', fill=T)
    fold.gd.file = file.path(state.dir, paste0('.fold.', fold, '.GD'))
    if (resume & file.exists(fold.gd.file)){
      cat(paste0("\n\n Found fold GDs, loading them\n"), fill=T)
      GDs <- readRDS(fold.gd.file)
    } else {
      # options(future.globals.maxSize=2000*1024^2) # 2000mb limit may be needed
      GDs <- predictGD(model,
                       newdata = test_fold_data,
                       verbose=verbose,
                       segmentation=segmentation,
                       segmentation_target=segmentation_target,
                       mask=mask,
                       afold=test_indices,
                       resume=resume,
                       save_states=save_states,
                       loginfo=c(fold, state.dir))

      if (save_states){saveRDS(GDs, fold.gd.file)}
      # update and save status on the registry
      registry[2, fold] = TRUE
      if (save_states){saveRDS(registry, reg.file)}
    }


    # get statistics for validation Global Deviance
    cat('-Summarizing statistics \n', fill=T)
    stats <- statGD(GDs,
                    k.penalty,
                    deg.fre=model[[1]]$df,
                    return_all_GD=return_all_GD)
    # update and save status on the registry
    registry[3, fold] = TRUE


    # Store the model and validation set
    cvresults[[fold]] <- stats
    fold.resuts.file = file.path(state.dir, paste0('.fold.', fold, '.results'))
    if (! file.exists(fold.resuts.file)){
      if (save_states){saveRDS(cvresults, fold.resuts.file)}
      # update and save status on the registry
      registry[4, fold] = TRUE
      if (save_states){saveRDS(registry, reg.file)}
    }


    # clean
    rm(model)
    rm(GDs)
    gc()
  }
  return(cvresults)
}








# Predict new fold Global Deviance
predictGD <- function (object, newdata = NULL, verbose=F,
                       segmentation=NULL, segmentation_target=NULL,
                       mask=NULL, afold=NULL,
                       loginfo=c(fold, state.dir),
                       resume=T, save_states=T, ...) {

  if (is.null(newdata)){stop("newdata is not set")}
  .fold = loginfo[1]
  .state.dir = loginfo[2]

  ## predict GD ##
  fname <- as.character(object[[1]]$family)
  familyobj <- gamlss2:::complete_family(get(fname))

  # .predicted.parameters
  quite(cat('Predicting fold parameters', fill=T), skip=verbose)
  fold.P.file = file.path(.state.dir, paste0('.fold.', .fold, '.predicted.parameters'))
  if (resume & file.exists(fold.P.file)){
    cat(paste0(" Found fold predicted parameters, loading them"), fill=T)
    nfitted <- readRDS(fold.P.file)
  } else {
    nfitted <- quite(
                     predict.vbgamlss(object,
                                      newdata = newdata,
                                      ptype='parameter',
                                      segmentation=segmentation,
                                      segmentation_target=segmentation_target,
                                      mask=mask,
                                      afold=afold
                                      ),
                     skip=verbose)
    if (save_states){saveRDS(nfitted, fold.P.file)}
  }


  # .predicted.response
  quite(cat('Predicting fold response', fill=T), skip=verbose)
  fold.R.file = file.path(.state.dir, paste0('.fold.', .fold, '.predicted.response'))
  if (resume & file.exists(fold.R.file)){
    cat(paste0(" Found fold predicted parameters, loading them"), fill=T)
    resp <- readRDS(fold.R.file)
  } else {
  resp <- quite(
                predict.vbgamlss(object,
                                 newdata = newdata,
                                 ptype='response',
                                 segmentation=segmentation,
                                 segmentation_target=segmentation_target,
                                 mask=mask,
                                 afold=afold
                                 ),
                skip=verbose)
  # save state if needed
  if (save_states){saveRDS(resp, fold.R.file)}
  }

  # add response to nfitted (from gamlss2.predict)
  nsub <- dim(newdata)[1] # == length(nfitted[[1]]$mu)
  for(i in 1:length(nfitted)){
    condA = ! "try-error" %in% class(nfitted[[i]])
    condB = ! all(is.na(nfitted[[i]]))
    if (condA & condB) {
        nfitted[[i]][['y']] <- as.numeric(unlist(resp[[i]])[1:nsub])
      } else {
      nfitted[[i]] <- NA
      }
  }
  rm(resp)
  rm(object)

  # missfits
  not_missfits <- ! is.na(nfitted)

  # test GD
  quite(cat('Predicting test fold GD ', fill=T), skip=verbose)
  plan(strategy="future::cluster", workers=availableCores())
  GDs <- foreach(i=seq_along(nfitted)) %dofuture% {
    #TRY(
    if (not_missfits[i]) {
        vxlGD <- testGD(nfitted[[i]], familyobj)
        vxlGD$vxl <- nfitted[[i]]$vxl
        vxlGD
      } else {
        NA
      }#,
    #logfile="/home/localadmin/urblauna/tpavan1/scripts_tommaso/test_batchtools/.vbgamlss.slurm/.bQCfBFCbek_2024-07-30/.1_ayIYAlblha/gdhelp/err.txt",
    #save.env.and.stop=T)
  }
  class(GDs) <- "vbgamlss.predictions.GD"
  GDs
}









# Test new fold Global Deviance
testGD <- function(nfit, familyobj){

  ## Internal of predictGD ##
  dfun <- paste("d", familyobj$family, sep = "")
  pfun <- paste("p", familyobj$family, sep = "")
  lpar <- length(familyobj$names)

  # NaN check for params missfits (ADDED MAE, LL, CLL)
  if (all(is.nan(nfit$y))) {
    out <- list(TGD = NA, predictError = NA, resid = nfit$mu*NA, MAE = NA, LL = NA, CLL = NA)
    return(out)
  }
  # Na check for params missfits (ADDED MAE, LL, CLL)
  if (all(is.na(nfit$y))) {
    out <- list(TGD = NA, predictError = NA, resid = nfit$mu*NA, MAE = NA, LL = NA, CLL = NA)
    return(out)
  }

  # subset NA
  yisnan <- is.nan(nfit$y)
  for (mtd in c('y', 'mu', 'sigma', 'nu', 'tau')[1:(lpar+1)]) {
    eval(parse(text=paste0('nfit$',mtd,' <- nfit$',mtd,'[! yisnan]')))
  }

  # Compute TGD on cleaned data (code from GAMLSS)
  if (is.null(nfit$y))
    stop("the response variables is missing in the newdata")

  if (familyobj$family %in% .gamlss.bi.list) {
    if (NCOL(nfit$y) == 1) {
      y1 <- nfit$y
      bd <- nfit$bd
    }
    else {
      bd <- nfit$y[, 1] + nfit$y[, 2]
      y1 <- nfit$y[, 1]
    }
  } else {
    y1 <- nfit$y
  }

  if (lpar == 1) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfit$mu, bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfit$mu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu)
    }
  } else if (lpar == 2) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$nmu, sigma = nfit$sigma, bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma)
    }
  } else if (lpar == 3) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu)
    }
  } else {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, tau = nfit$tau, bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, tau = nfit$tau, bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, tau = nfit$tau, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, tau = nfit$tau)
    }
  }

  Vresid <- qNO(eval(ures))
  dev <- -2 * sum(eval(devi))

  # Mean Absolute Error (MAE)
  vxl_mae <- mean(abs(y1 - nfit$mu), na.rm = TRUE)

  # Extract observation-wise log-loss (negative log-likelihood) and CDF
  ll_obs_raw <- -eval(devi)
  cdf_vals <- eval(ures)

  # Uncensored Mean Log-Loss (LL)
  vxl_ll <- mean(ll_obs_raw, na.rm = TRUE)

  # Censored Mean Log-Loss (CLL)
  ll_obs_censored <- ll_obs_raw
  ll_obs_censored[cdf_vals < 0.01 | cdf_vals > 0.99] <- -log(0.02)
  vxl_cll <- mean(ll_obs_censored, na.rm = TRUE)

  # Recompose subsetted NaNs
  Vresid_na <- yisnan * NA
  Vresid_na[! yisnan] <- Vresid

  # output
  out <- list(TGD = dev,
              predictError = dev/length(nfit$mu),
              resid = Vresid_na,
              MAE = vxl_mae,
              LL = vxl_ll,
              CLL = vxl_cll)
  return(out)
}







# Summarize fold Global Deviance statistics
statGD <- function(GDs, k.penalty=NULL, deg.fre=1, return_all_GD=F) {
  missfits <- sum(is.na(GDs))
  nsub <- length(GDs[[1]]$resid)
  nvxl <- length(GDs)

  TGDs <- numeric(nvxl) * NA
  predictErrors <- numeric(nvxl) * NA
  MAEs <- numeric(nvxl) * NA
  LLs <- numeric(nvxl) * NA
  CLLs <- numeric(nvxl) * NA

  resids <- matrix(data=NA, nrow=nsub, ncol=nvxl)
  resids_4math <- c()

  for (i in seq_len(nvxl)) {
    if (! is.na(GDs)[i]) {
      TGDs[i] <- GDs[[i]]$TGD
      predictErrors[i] <- GDs[[i]]$predictError
      MAEs[i] <- GDs[[i]]$MAE
      LLs[i] <- GDs[[i]]$LL
      CLLs[i] <- GDs[[i]]$CLL
      resids[, i] <- GDs[[i]]$resid
      resids_4math <- c(resids_4math, GDs[[i]]$resid)
    } else {
      TGDs[i] <- NA
      predictErrors[i] <- NA
      MAEs[i] <- NA
      LLs[i] <- NA
      CLLs[i] <- NA
      resids[, i] <- NA
      resids_4math <- NA
    }
  }

  # penalize TGDs
  if (! is.null(k.penalty)){
    TGDs_pen <- (k.penalty * deg.fre) - TGDs
  }

  # statistics
  GD = describe_stats(TGDs, '')
  if (! is.null(k.penalty)){
    GDp = describe_stats(TGDs_pen, '')
  }
  PE = describe_stats(predictErrors, '')

  # Summarize metrics across all voxels
  MAE_stat = describe_stats(MAEs, '')
  LL_stat = describe_stats(LLs, '')
  CLL_stat = describe_stats(CLLs, '')

  resids_4math[is.infinite(resids_4math)] <- NA
  RS = describe_stats(resids_4math, '')

  # return
  to_return = list(GD=GD, GDpen=GDp, predErr=PE, Resids=RS,
                   MAE=MAE_stat, LL=LL_stat, CLL=CLL_stat,
                   missFits=missfits, missFitsPerc=missfits/length(GDs))

  if (return_all_GD) {
    to_return <- append(to_return, list(allGD=TGDs,
                                        allPredErr=predictErrors,
                                        allMAE=MAEs,
                                        allLL=LLs,
                                        allCLL=CLLs,
                                        allResid=resids))
  }
  return(to_return)
}










# Descriptive statistics helper
describe_stats <- function(x, vname) {
  out <- c()
  x <- na.omit(x)
  for (fn in c('mean', 'sd', 'quantile', 'min', 'max')) {
    eval(parse(text=paste0(fn, vname,'<- ', fn,'(x, na.rm=T)')))
    eval(parse(text=paste0('out <- c(out,', fn, vname, ' = ', fn, vname, ')')))
  }
  return(out)
}







# Create leave-one-site/dataset-out cross-validation folds
#' @export
LOSOfolds <- function(df){
  if (is.factor(df)){
    return(as.integer(df))
  } else {
     errorCondition('input must be a factor.')
    }
}








# Create stratified cross-validation folds
#' @export
stratCVfolds <- function(df, k.fold=10){
  folds_row <- rep(0, length(df))
  folds_col <- create_folds(df,
                            k = k.fold,
                            type='stratified',
                            invert=T)
  for (i in seq(length(folds_col))) {
    folds_row[unlist(folds_col[i])] <- i
  }
  return(folds_row)
}








# Lazy get cross-validated Global Deviance from cvresults
getCVGD <- function(cvresults, term='mean') {
  # mean, sd,
  # quantile.0%, quantile.25%, quantile.50%, quantile.75%, quantile.100%,
  # min, max
  CVGD = 0
  for (fold in cvresults) {CVGD <- CVGD + fold$GD[[term]]}
  return(CVGD)
}








# Lazy get cross-validated penalized Global Deviance from cvresults
getCVGD.pen <- function(cvresults, term='quantile.50%') {
  CVGD = 0
  for (fold in cvresults) {CVGD <- CVGD + fold$GDpen[[term]]}
  return(CVGD)
}








# Lazy get all (pen and non) cross-validated Global Deviance from cvresults
getCVGD.all <- function(cvresults, penalized=F) {
  out <- c()
  for (t in c('mean', 'sd', 'quantile.0%', 'quantile.25%', 'quantile.50%',
              'quantile.75%', 'quantile.100%', 'min', 'max')) {
    if (penalized) {
    out[[t]] <- getCVGD.pen(cvresults, term=t)
    } else {
      out[[t]] <- getCVGD(cvresults, term=t)
    }
  }
  return(out)
}






# Lazy get cross-validated Mean Absolute Error (MAE) from cvresults
getCVMAE <- function(cvresults, term='mean') {
  CVMAE = 0
  for (fold in cvresults) {
    # If the fold failed entirely, skip it to prevent crashing
    if (!is.null(fold$MAE[[term]])) {
      CVMAE <- CVMAE + fold$MAE[[term]]
    }
  }
  # Note: A true cross-validated average MAE requires dividing by the number of folds
  return(CVMAE / length(cvresults))
}






# Lazy get cross-validated Uncensored Log-Loss (LL) from cvresults
getCVLL <- function(cvresults, term='mean') {
  CVLL = 0
  for (fold in cvresults) {
    if (!is.null(fold$LL[[term]])) {
      CVLL <- CVLL + fold$LL[[term]]
    }
  }
  return(CVLL / length(cvresults))
}






# Lazy get cross-validated Censored Log-Loss (CLL) from cvresults
getCVCLL <- function(cvresults, term='mean') {
  CVCLL = 0
  for (fold in cvresults) {
    if (!is.null(fold$CLL[[term]])) {
      CVCLL <- CVCLL + fold$CLL[[term]]
    }
  }
  return(CVCLL / length(cvresults))
}





# (Optional) Lazy get ALL cross-validated metrics at once
getCV_All_Metrics <- function(cvresults) {
  out <- list()
  for (t in c('mean', 'sd', 'quantile.0%', 'quantile.25%', 'quantile.50%',
              'quantile.75%', 'quantile.100%', 'min', 'max')) {
    out$GD[[t]] <- getCVGD(cvresults, term=t)
    out$MAE[[t]] <- getCVMAE(cvresults, term=t)
    out$LL[[t]] <- getCVLL(cvresults, term=t)
    out$CLL[[t]] <- getCVCLL(cvresults, term=t)
  }
  return(out)
}





# Akaike weights from a vector of AICc or similar values
#' @export
akaike_weights <- function(v){
  # https://link.springer.com/content/pdf/10.3758/BF03206482.pdf
  # Dinga et al 2021
  v <- v - min(v)
  d <- exp(-0.5 * v)
  return(d/sum(d))
}










#===============================================================================
#                                 TESTING                                      #
#===============================================================================


#' @export
vbgamlss.nested_cv <- function(image,
                               mask,
                               g.formula,
                               full.data,
                               site.var, # <- dataset variables
                               fold.var,
                               g.family = NO,
                               segmentation = NULL,
                               segmentation_target=NULL,
                               num_cores = NULL,
                               chunk_max_mb = 64,
                               n_folds = 10,
                               k.penalty=NULL,
                               verbose=F,
                               debug=T,
                               logdir=getwd(),
                               return_all_GD=T,
                               save_states=T,
                               resume=T,
                               ...) {

  # Identify all unique sites for the Outer Loop (LOSO)
  sites <- unique(full.data[[site.var]])
  nested_results <- list()

  for (site in sites) {
    cat(paste0("\n=======================================================\n"), fill=T)
    cat(paste0(" OUTER LOOP: Leaving out Dataset/Site '", site, "' for Testing"), fill=T)
    cat(paste0("=======================================================\n"), fill=T)

    # Split Data: Outer Train vs. Outer Test
    train_indices <- full.data[[site.var]] != site
    test_indices <- full.data[[site.var]] == site

    outer_train_data <- full.data[train_indices, ]
    outer_test_data <- full.data[test_indices, ]

    # Create a site-specific directory to prevent state overwriting across outer loops
    site_logdir <- file.path(logdir, paste0("Nested_LOSO_", site))
    dir.create(site_logdir, recursive = TRUE, showWarnings = FALSE)

    # ---------------------------------------------------------
    # INNER LOOP: Run your existing k-fold CV on the outer training set
    # ---------------------------------------------------------
    cat(paste0("\n--- Starting Inner k-fold CV (Training Set Only) ---\n"), fill=T)
    inner_cv_results <- vbgamlss.cv(image = image,
                                    mask = mask,
                                    g.formula = g.formula,
                                    train.data = outer_train_data,
                                    fold.var = fold.var,
                                    logdir = site_logdir, # Saves .cvstates strictly inside this site's folder
                                    save_states = save_states,
                                    g.family = g.family,
                                    segmentation = segmentation,
                                    segmentation_target=segmentation_target,
                                    num_cores = num_cores,
                                    chunk_max_mb = chunk_max_mb,
                                    n_folds = n_folds,
                                    k.penalty = k.penalty,
                                    verbose=F,
                                    debug=T,
                                    return_all_GD=T,
                                    resume=T,
                                    ...
    )

    # ---------------------------------------------------------
    # OUTER FIT: Train a final model on the ENTIRE outer_train_data
    # ---------------------------------------------------------
    cat(paste0("\n--- Fitting Final Outer Model (All Training Sites) ---\n"), fill=T)
    outer_model_file <- file.path(site_logdir, "final_outer_model.rds")

    if (save_states && file.exists(outer_model_file)) {
      cat("Found existing outer model, loading it...", fill=T)
      outer_model <- readRDS(outer_model_file)
    } else {
      outer_model <- vbgamlss(
        image = image,
        mask = mask,
        g.formula = as.formula(g.formula),
        train.data = outer_train_data,
        logdir = site_logdir,
        afold = rep(TRUE, nrow(outer_train_data)), # Use all outer_train_data
        ...
      )
      if (save_states) saveRDS(outer_model, outer_model_file)
    }

    # ---------------------------------------------------------
    # OUTER TEST: Evaluate the final model on the Held-Out Site
    # ---------------------------------------------------------
    cat(paste0("\n--- Evaluating Model on Held-Out Site '", site, "' ---\n"), fill=T)
    outer_gd_file <- file.path(site_logdir, "final_outer_GDs.rds")

    if (save_states && file.exists(outer_gd_file)) {
      cat("Found existing outer test statistics, loading them...", fill=T)
      outer_GDs <- readRDS(outer_gd_file)
    } else {
      outer_GDs <- predictGD(
        object = outer_model,
        newdata = outer_test_data,
        mask = mask,
        afold = which(test_indices), # The row indices in the full dataset
        loginfo = c(paste0("TestSite_", site), site_logdir),
        save_states = save_states,
        ...
      )
      if (save_states) saveRDS(outer_GDs, outer_gd_file)
    }

    # Summarize the outer test performance using your updated statGD function
    outer_test_stats <- statGD(outer_GDs)

    # ---------------------------------------------------------
    # Store Results
    # ---------------------------------------------------------
    nested_results[[as.character(site)]] <- list(
      Site = site,
      N_Train = nrow(outer_train_data),
      N_Test = nrow(outer_test_data),
      Inner_CV = inner_cv_results,
      Outer_Test = outer_test_stats
    )

    # Clean up heavy objects before next iteration
    rm(outer_model, outer_GDs)
    gc()
  }

  cat("\n=======================================================\n", fill=T)
  cat(" Nested Cross-Validation Complete!", fill=T)
  cat("=======================================================\n", fill=T)

  return(nested_results)
}




#===============================================================================
#                                 DEPRECATED                                   #
#===============================================================================

# Test new fold Global Deviance
# testGD <- function(nfit, familyobj){
#
#   ## Internal of predictGD ##
#   dfun <- paste("d", familyobj$family, sep = "")
#   pfun <- paste("p", familyobj$family, sep = "")
#   lpar <- length(familyobj$names)
#
#   # NaN check for params missfits
#   if (all(is.nan(nfit$y))) {
#     out <- list(TGD = NA, predictError = NA, resid = nfit$mu*NA)
#     return(out)
#   }
#   # Na check for params missfits
#   if (all(is.na(nfit$y))) {
#     out <- list(TGD = NA, predictError = NA, resid = nfit$mu*NA)
#     return(out)
#   }
#
#   # subset NA
#   yisnan <- is.nan(nfit$y)
#   for (mtd in c('y', 'mu', 'sigma', 'nu', 'tau')[1:(lpar+1)]) {
#     eval(parse(text=paste0('nfit$',mtd,' <- nfit$',mtd,'[! yisnan]')))
#   }
#
#   # Compute TGD on cleaned data (code from GAMLSS)
#   if (is.null(nfit$y))
#     stop("the response variables is missing in the newdata")
#
#   if (familyobj$family %in% .gamlss.bi.list) {
#     if (NCOL(nfit$y) == 1) {
#       y1 <- nfit$y
#       bd <- nfit$bd
#     }
#     else {
#       bd <- nfit$y[, 1] + nfit$y[, 2]
#       y1 <- nfit$y[, 1]
#     }
#   } else {
#     y1 <- nfit$y
#   }
#
#   if (lpar == 1) {
#     if (familyobj$family %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfit$mu, bd = bd,
#                    log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$mu, bd = bd)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfit$mu, log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$mu)
#     }
#   }
#   else if (lpar == 2) {
#     if (familyobj$family %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    bd = bd, log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$nmu, sigma = nfit$sigma,
#                    bd = bd)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma)
#     }
#   }
#   else if (lpar == 3) {
#     if (familyobj$family %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu, bd = bd, log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu, bd = bd)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu, log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu)
#     }
#   }
#   else {
#     if (familyobj$family %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu, tau = nfit$tau, bd = bd,
#                    log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu, tau = nfit$tau, bd = bd)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu, tau = nfit$tau, log = TRUE)
#       ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
#                    nu = nfit$nu, tau = nfit$tau)
#     }
#   }
#
#   Vresid <- qNO(eval(ures))
#   dev <- -2 * sum(eval(devi))
#
#   # Recompose subsetted NaNs
#   Vresid_na <- yisnan * NA
#   Vresid_na[! yisnan] <- Vresid
#
#   # output
#   out <- list(TGD = dev,
#               predictError = dev/length(nfit$mu),
#               resid = Vresid_na)
#   return(out)
#
# }


# # Summarize fold Global Deviance statistics
# statGD <- function(GDs, k.penalty=NULL, deg.fre=1, return_all_GD=F) {
#   missfits <- sum(is.na(GDs))
#   nsub <- length(GDs[[1]]$resid)
#   nvxl <- length(GDs)
#   TGDs <- numeric(nvxl) * NA
#   predictErrors <- numeric(nvxl) * NA
#   resids <- matrix(data=NA, nrow=nsub, ncol=nvxl)
#   resids_4math <- c()
#
#   for (i in seq_len(nvxl)) {
#     if (! is.na(GDs)[i]) {
#       TGDs[i] <- GDs[[i]]$TGD
#       predictErrors[i] <- GDs[[i]]$predictError
#       resids[, i] <- GDs[[i]]$resid
#       resids_4math <- c(resids_4math, GDs[[i]]$resid)
#     } else {
#       TGDs[i] <- NA
#       predictErrors[i] <- NA
#       resids[, i] <- NA
#       resids_4math <- NA
#     }
#   }
#
#   # penalize TGDs
#   if (! is.null(k.penalty)){
#     TGDs_pen <- (k.penalty * deg.fre) - TGDs
#   }
#
#   # statistics
#   GD = describe_stats(TGDs, '')
#   if (! is.null(k.penalty)){
#     GDp = describe_stats(TGDs_pen, '')
#   }
#   PE = describe_stats(predictErrors, '')
#   resids_4math[is.infinite(resids_4math)] <- NA
#   RS = describe_stats(resids_4math, '')
#
#   to_return = list(GD=GD, GDpen=GDp, predErr=PE, Resids=RS,
#                    missFits=missfits, missFitsPerc=missfits/length(GDs))
#
#   if (return_all_GD) {to_return <- append(to_return, list(allGD=TGDs,
#                                                           allPredErr=predictErrors,
#                                                           allResid=resids))}
#   return(to_return)
# }

