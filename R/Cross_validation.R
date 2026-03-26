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
#' @export
vbgamlss.cv <- function(imageframe,
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
                        drop_re=T,
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
    test_fold_data <- droplevels(test_fold_data)

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
                 num_cores=num_cores,
                 chunk_max_mb=chunk_max_mb,
                 #afold=training_fold,
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
      GDs <- predictGD(model,
                       test_imageframe     = imageframe[test_indices,],
                       newdata             = test_fold_data,
                       verbose             = verbose,
                       segmentation        = segmentation,
                       segmentation_target = segmentation_target,
                       #afold               = test_indices,
                       resume              = resume,
                       save_states         = save_states,
                       drop_re             = drop_re,
                       loginfo             = c(fold, state.dir))

      if (save_states){qs2::qs_save(GDs, fold.gd.file)}

      # update and save status on the registry
      registry[2, fold] = TRUE
      if (save_states){qs2::qs_save(registry, reg.file)}
    }

    # get statistics for validation Global Deviance
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

  # Give a name to each fold
  names(cvresults) <- paste0("fold", seq_along(cvresults))

  # save the final result
  if (save_states){
    cvres_filepath <- paste0(state.dir, 'vbgamlss.cvresults')
    qs2::qs_save(cvresults, file = cvres_filepath)
    cat('Saved CV results: ', cvres_filepath, '\n')
  }

  return(cvresults)
}




# --------------------------------
# Predict new fold Global Deviance
predictGD <- function (object,
                       test_imageframe,
                       newdata = NULL,
                       verbose=F,
                       segmentation=NULL,
                       segmentation_target=NULL,
                       #afold=NULL,
                       loginfo=c(fold, state.dir),
                       resume=T,
                       save_states=T,
                       drop_re=T,
                       ...) {

  if (is.null(newdata)){stop("newdata is not set")}
  .fold = loginfo[1]
  .state.dir = loginfo[2]

  ## predict GD ##
  fname <- as.character(object[[1]]$family)
  familyobj <- gamlss2:::complete_family(get(fname))



  # .predicted.parameters
  fold.P.file = file.path(.state.dir, paste0('.fold.', .fold, '.predicted.parameters.qs'))
  if (resume & file.exists(fold.P.file)){
    cat(paste0("- - found existing test fold predicted parameters, loading them"), fill=T)
    nfitted <- qs2::qs_read(fold.P.file)
  } else {
    quite(cat('\t| predicting test fold parameters', fill=T), skip=verbose)
    # Remove random effects if unseen factor level is predicted like in LOSO ---
    if (drop_re) {

      # use first model to get formula and exclude RE
      terms_exre <- get_fixed_terms(object[[1]])
      m <- restore_family(object[[1]])
      parnames <- m$family$names

      # make output df
      list_nfitted <- list()

      # Predict per parameter feeding terms without RE.
      for (param in parnames){
        quite(cat('\t\t| ', param, fill=T), skip=verbose)
        cat("\033[34m")
        pfit <- quite(
          predict.vbgamlss(object,
                           newdata             = newdata,
                           ptype               = 'parameter',
                           segmentation        = segmentation,
                           segmentation_target = segmentation_target,
                           #afold               = afold,
                           what                = param,
                           terms               = terms_exre
          ),
          skip=verbose)
        cat("\033[0m")
        # store it
        list_nfitted[[param]] <- pfit
      }
      # store them properly
      nfitted <- list_nfitted[['mu']]
      for (n in parnames) {
        if (n == 'mu') {next}
        nfitted <- mapply(c, nfitted, list_nfitted[[n]], SIMPLIFY = FALSE)
      }

    } else {
      # Just predict, do not drop re if e.g. there are no unseen levels ---
      cat("\033[34m")
      nfitted <- quite(
        predict.vbgamlss(object,
                         newdata             = newdata,
                         ptype               = 'parameter',
                         segmentation        = segmentation,
                         segmentation_target = segmentation_target,
                         #afold               = afold
        ),
        skip=verbose)
      cat("\033[0m")
    }
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

    RE_to_drop=NULL
    if (drop_re){
      terms_exre <- get_fixed_terms(object[[1]])
    }
    cat("\033[34m")
    resp <- quite(
      predict.vbgamlss(object,
                       newdata             = newdata,
                       ptype               = 'response',
                       segmentation        = segmentation,
                       segmentation_target = segmentation_target,
                       #afold               = afold,
                       terms               = terms_exre
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
  quite(cat('\t| Predicting test fold GD ', fill=T), skip=verbose)
  plan(strategy="future::cluster", workers=availableCores())
  GDs <- foreach(i=seq_along(nfitted)) %dofuture% {
    if (not_missfits[i]) {
      vxlGD <- testGD(nfitted[[i]], familyobj)
      vxlGD$vxl <- nfitted[[i]]$vxl
      vxlGD
    } else {
      NA
    }
  }
  class(GDs) <- "vbgamlss.predictions.GD"
  GDs
}



# --------------------------------
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
  # https://github.com/gamlss-dev/gamlss/blob/6a1e3c34ebf55825fb4774c60206c777442da21c/R/gamlssVGD_23_12_21.R#L235
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

  # Mean Absolute Error (MAE) using the true predicted mean
  vxl_mae <- mean(abs(y1 - nfit$yhat), na.rm = TRUE)

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
  out <- list(TGD          = dev,
              predictError = dev/length(nfit$mu),
              resid        = Vresid_na,
              MAE          = vxl_mae,
              LL           = vxl_ll,
              CLL          = vxl_cll)
  return(out)
}



# --------------------------------
# Summarize fold Global Deviance statistics
statGD <- function(GDs, k.penalty=NULL, deg.fre=1, return_all_GD=F) {
  missfits <- sum(is.na(GDs))
  nsub <- length(GDs[[1]]$resid)
  nvxl <- length(GDs)

  TGDs <- numeric(nvxl) * NA
  predictErrors <- numeric(nvxl) * NA
  MAEs <- numeric(nvxl) * NA
  LLs <- numeric(nvxl)  * NA
  CLLs <- numeric(nvxl) * NA

  resids <- matrix(data=NA, nrow=nsub, ncol=nvxl)
  resids_4math <- c()

  for (i in seq_len(nvxl)) {
    if (! is.na(GDs)[i]) {
      TGDs[i]          <- GDs[[i]]$TGD
      predictErrors[i] <- GDs[[i]]$predictError
      MAEs[i]          <- GDs[[i]]$MAE
      LLs[i]           <- GDs[[i]]$LL
      CLLs[i]          <- GDs[[i]]$CLL
      resids[, i]      <- GDs[[i]]$resid
      resids_4math     <- c(resids_4math, GDs[[i]]$resid)
    } else {
      TGDs[i]          <- NA
      predictErrors[i] <- NA
      MAEs[i]          <- NA
      LLs[i]           <- NA
      CLLs[i]          <- NA
      resids[, i]      <- NA
      resids_4math     <- NA
    }
  }

  # penalize TGDs
  if (! is.null(k.penalty)){
    TGDs_pen <- (k.penalty * deg.fre) - TGDs
  }

  # statistics
  GD = describe_stats(TGDs, '')
  GDp = NULL
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



# --------------------------------
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





# --------------------------------
# Return list of term for each gamlss2 formula without random effects in them
#' @export
get_fixed_terms <- function(m) {
  m <- restore_family(m)

  fixed_terms <- c("(Intercept)")

  for (param in m$family$names) {
    # Add all linear terms
    if (length(m$xterms[[param]]) > 0) {
      fixed_terms <- c(fixed_terms, m$xterms[[param]])
    }

    # Evaluate special/smooth terms
    s_terms <- m$sterms[[param]]
    if (length(s_terms) > 0) {
      for (t in s_terms) {
        is_re <- FALSE

        # mgcv random effect (bs="re")?
        if (!is.null(m$specials[[t]]) && inherits(m$specials[[t]], "random.effect")) {
          is_re <- TRUE}

        # legacy gamlss random coefficient?
        if (!is_re && !is.null(m$fitted.specials[[param]][[t]]$coefficients)) {
          if (inherits(m$fitted.specials[[param]][[t]]$coefficients, "random")) {
            is_re <- TRUE}
        }

        # string fallbacks (just in case)
        is_random_effect <- grepl("bs\\s*=\\s*['\"]re['\"]", t) |
          grepl("^random\\(", t) |
          grepl("^re\\(", t)
        if (!is_re && is_random_effect) {
          is_re <- TRUE }

        # If it passed all checks, it's a fixed smooth/effect. Keep it.
        if (!is_re) {
          fixed_terms <- c(fixed_terms, t)
        }
      }
    }
  }

  # Return a unique, flat vector of all safe terms
  return(unique(fixed_terms))
}





# --------------------------------
# Create leave-one-site/dataset-out cross-validation folds
#' @export
LOSOfolds <- function(v){
  if (is.factor(v)){
    return(as.integer(v))
  } else {
    errorCondition('input must be a factor.')
  }
}



# --------------------------------
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



# --------------------------------
# Lazy get cross-validated Global Deviance from cvresults
getCVGD <- function(cvresults, term='mean') {
  CVGD = 0
  for (fold in cvresults) {CVGD <- CVGD + fold$GD[[term]]}
  return(CVGD)
}



# --------------------------------
# Lazy get cross-validated penalized Global Deviance from cvresults
getCVGD.pen <- function(cvresults, term='quantile.50%') {
  CVGD = 0
  for (fold in cvresults) {CVGD <- CVGD + fold$GDpen[[term]]}
  return(CVGD)
}



# --------------------------------
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



# --------------------------------
# Lazy get cross-validated Mean Absolute Error (MAE) from cvresults
getCVMAE <- function(cvresults, term='mean') {
  CVMAE = 0
  for (fold in cvresults) {
    # If the fold failed entirely, skip it to prevent crashing
    if (!is.null(fold$MAE[[term]])) {
      CVMAE <- CVMAE + fold$MAE[[term]]
    }
  }
  return(CVMAE / length(cvresults))
}



# --------------------------------
# Lazy get cross-validated Uncensored Log-Loss (LL) or Negative Log-Predictive Density (NLPD) from cvresults
getCVLL <- function(cvresults, term='mean') {
  CVLL = 0
  for (fold in cvresults) {
    if (!is.null(fold$LL[[term]])) {
      CVLL <- CVLL + fold$LL[[term]]
    }
  }
  return(CVLL / length(cvresults))
}



# --------------------------------
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



# --------------------------------
# Lazy get ALL cross-validated metrics at once
getCV_All_Metrics <- function(cvresults) {
  out <- list(GD = list(), MAE = list(), LL = list(), CLL = list())
  for (t in c('mean', 'sd', 'quantile.0%', 'quantile.25%', 'quantile.50%',
              'quantile.75%', 'quantile.100%', 'min', 'max')) {
    out$GD[[t]] <- getCVGD(cvresults, term=t)
    out$MAE[[t]] <- getCVMAE(cvresults, term=t)
    out$LL[[t]] <- getCVLL(cvresults, term=t)
    out$CLL[[t]] <- getCVCLL(cvresults, term=t)
  }
  return(out)
}



# --------------------------------
# Akaike weights from a vector of AICc or similar values
#' @export
akaike_weights <- function(v){
  v <- v - min(v)
  d <- exp(-0.5 * v)
  return(d/sum(d))
}






# --------------------------------



#===============================================================================
#                                 TESTING                                      #
#===============================================================================






#===============================================================================
#                                 DEPRECATED                                   #
#===============================================================================

get_fixed_terms_old <- function(m) {
  m <- restore_family(m)
  pnames <- m$family$names
  form_list <- list()
  for (i in seq_along(pnames)) {
    form_list[[pnames[i]]] <- formula(m$formula, rhs = i)
  }

  # Identify re for each formula
  clean_terms_list <- list()
  for (param in names(form_list)) {

    # Extract the term labels
    term_obj <- terms(form_list[[param]])
    all_terms <- attr(term_obj, "term.labels")

    # If not just an intercept ~1
    if (length(all_terms) > 0) {

      # - bs\\s*=\\s*['\"]re['\"] : Catches bs='re' or bs="re" with any spacing
      # - ^random\\( | ^re\\(     : Catches legacy gamlss random effect wrappers
      is_random_effect <- grepl("bs\\s*=\\s*['\"]re['\"]", all_terms) |
        grepl("^random\\(", all_terms) |
        grepl("^re\\(", all_terms)

      # Keep only the terms that are not re
      clean_terms_list[[param]] <- c("(Intercept)", all_terms[!is_random_effect])

    } else {
      # return an empty character vector if intercept
      clean_terms_list[[param]] <- "(Intercept)"
    }
  }

  return(clean_terms_list)
}
