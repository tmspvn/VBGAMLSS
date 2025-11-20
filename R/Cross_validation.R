#############################  cross validation  ###############################








#' @export
vbgamlss.cv <- function(image,
                        mask,
                        g.formula,
                        train.data,
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
                        save_states=T, resume=T,
                        ...) {

  cat(paste0("Starting"), fill=T)

  if (! is.null(segmentation)){
    if (segmentation == 'NULL') {segmentation=NULL}}
  if (missing(image)) { stop("image is missing") }
  if (missing(mask)) { stop("mask is missing") }
  if (missing(fold.var)) { stop("variable to use in the stratified CV is missing") }
  if (missing(g.formula)) { stop("formula is missing")}
  check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}


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
    train.data$folds <- stratCVfolds(fold.var, k = n_folds)
    if (save_states){saveRDS(train.data$folds, fold.file)}
  }

  # DEPRECATED
  # subsample?
  # if (! is.null(subsample.factor)){
  #   subsample_ <- get_subsample_indices(sum(antsImageRead(mask)>0),
  #                                      subsample.type,
  #                                      subsample.factor)
  # } else {subsample_ <- NULL}


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

  # NaN check for params missfits
  if (all(is.nan(nfit$y))) {
    out <- list(TGD = NA, predictError = NA, resid = nfit$mu*NA)
    return(out)
  }
  # Na check for params missfits
  if (all(is.na(nfit$y))) {
    out <- list(TGD = NA, predictError = NA, resid = nfit$mu*NA)
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
      devi <- call(dfun, x = y1, mu = nfit$mu, bd = bd,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfit$mu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu)
    }
  }
  else if (lpar == 2) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
                   bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$nmu, sigma = nfit$sigma,
                   bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma)
    }
  }
  else if (lpar == 3) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu, bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu, bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu)
    }
  }
  else {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu, tau = nfit$tau, bd = bd,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu, tau = nfit$tau, bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu, tau = nfit$tau, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma,
                   nu = nfit$nu, tau = nfit$tau)
    }
  }

  Vresid <- qNO(eval(ures))
  dev <- -2 * sum(eval(devi))

  # Recompose subsetted NaNs
  Vresid_na <- yisnan * NA
  Vresid_na[! yisnan] <- Vresid

  # output
  out <- list(TGD = dev, predictError = dev/length(nfit$mu),
              resid = Vresid_na)
  return(out)

}







# Summarize fold Global Deviance statistics
statGD <- function(GDs, k.penalty=NULL, deg.fre=1, return_all_GD=F) {
  missfits <- sum(is.na(GDs))
  nsub <- length(GDs[[1]]$resid)
  nvxl <- length(GDs)
  TGDs <- numeric(nvxl) * NA
  predictErrors <- numeric(nvxl) * NA
  resids <- matrix(data=NA, nrow=nsub, ncol=nvxl)
  resids_4math <- c()

  for (i in seq_len(nvxl)) {
    if (! is.na(GDs)[i]) {
      TGDs[i] <- GDs[[i]]$TGD
      predictErrors[i] <- GDs[[i]]$predictError
      resids[, i] <- GDs[[i]]$resid
      resids_4math <- c(resids_4math, GDs[[i]]$resid)
       } else {
      TGDs[i] <- NA
      predictErrors[i] <- NA
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
  resids_4math[is.infinite(resids_4math)] <- NA
  RS = describe_stats(resids_4math, '')

  to_return = list(GD=GD, GDpen=GDp, predErr=PE, Resids=RS,
                   missFits=missfits, missFitsPerc=missfits/length(GDs))

  if (return_all_GD) {to_return <- append(to_return, list(allGD=TGDs,
                                                          allPredErr=predictErrors,
                                                          allResid=resids))}
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








# Akaike weights from a vector of AICc or similar values
#' @export
akaike_weights <- function(v){
  # https://link.springer.com/content/pdf/10.3758/BF03206482.pdf
  # Dinga et al 2021
  v <- v - min(v)
  d <- exp(-0.5 * v)
  return(d/sum(d))
}






