#############################  cross validation  ###############################








vbgamlss.cv <- function(image, mask, g.formula, train.data, fold.var, g.family = NO,
                      segmentation = NULL, num_cores = NULL, chunk_max_mb = 64,
                      n_folds = 10, k.penalty=NULL, verbose=F,
                      return_all_GD=T, subsample.factor=NULL,
                      subsample.type=c('regular', 'random'),
                      ...) {
  if (missing(image)) { stop("image is missing") }
  if (missing(mask)) { stop("mask is missing") }
  if (missing(fold.var)) { stop("variable to use in the stratified CV is missing") }
  if (missing(g.formula)) { stop("formula is missing")}
  check_formula_LHS(g.formula)
  if (missing(train.data)) { stop("subjData is missing")}
  # Force character columns to factors
  train.data <- as.data.frame(unclass(train.data),stringsAsFactors=TRUE)

  # Create stratified folds
  train.data$folds <- stratCVfolds(fold.var, k = n_folds)
  # subsample?
  if (! is.null(subsample.factor)){
    subsample_ <- get_subsample_indices(sum(antsImageRead(mask)>0),
                                       subsample.type,
                                       subsample.factor)
  } else {subsample_ <- NULL}

  # Prepare storage for results
  cvresults <- c()

  # Perform stratified cross-validation
  for (fold in 1:n_folds) {
    cat(paste0("\n\nProcessing fold ", fold, " of ", n_folds, "\n"), fill=T)

    # Split data into training and validation sets
    training_fold = train.data$folds != fold
    test_indices <- which(train.data$folds == fold)
    test_fold_data <- train.data[test_indices, ]
    cat('-Fitting fold', fill=T)
    # Fit the model on the training fold
    model <- quite(vbgamlss_chunks(image,
                                   mask,
                                   as.formula(g.formula),
                                   train.data,
                                   g.family,
                                   segmentation,
                                   num_cores,
                                   chunk_max_mb=128,
                                   afold=training_fold, # pass fold indexes
                                   subsample=subsample_),
                   skip=verbose)
    # predict new fold
    cat('-Estimating GD', fill=T)
    # options(future.globals.maxSize=2000*1024^2) # 2000mb limit may be needed
    GDs <- predictGD(model, newdata = test_fold_data, verbose=verbose,
                     segmentation=segmentation,
                     mask=mask,
                     afold=test_indices,
                     subsample=subsample_)
    # get statistics for validation Global Deviance
    cat('-Summarizing statistics', fill=T)
    stats <- statGD(GDs,
                  k.penalty,
                  deg.fre=model[[1]]$df,
                  return_all_GD=return_all_GD)
    # Store the model and validation set
    cvresults[[fold]] <- stats
    # clean
    rm(model)
    rm(GDs)
    gc()
  }
  return(cvresults)
}


predictGD <- function (object, newdata = NULL, verbose=F,
                       segmentation=NULL, mask=NULL, afold=NULL, subsample=NULL, ...) {
  if (is.null(newdata)){stop("newdata is not set")}

  ## predict GD ##
  fname <- object$family[[1]]
  fname <- as.character(object[[1]]$family)
  familyobj <- gamlss2:::complete_family(get(fname))
  quite(cat('Predicting fold parameters', fill=T), skip=verbose)
  nfitted <- quite(predict.vbgamlss(object, newdata = newdata, ptype='parameter',
                                    segmentation=segmentation,
                                    mask=mask, afold=afold,
                                    subsample=subsample), skip=verbose)
  nsub <- length(nfitted[[1]]['mu'])
  # add response to nfitted (from gamlss2.predict)
  quite(cat('Predicting fold response', fill=T), skip=verbose)
  resp <- quite(predict.vbgamlss(object, newdata = newdata, ptype='response',
                                 segmentation=segmentation,
                                 mask=mask, afold=afold,
                                 subsample=subsample), skip=verbose)
  for(i in 1:length(object)){
    nfitted[[i]]['y'] <- unlist(resp[[i]][1:nsub])
  }
  rm(resp)
  # test GD
  quite(cat('Predicting test fold GD ', fill=T), skip=verbose)
  plan(cluster)
  GDs <- foreach(i=1:length(object)) %dofuture% {
    vxlGD <- testGD(nfitted[[1]], familyobj)
    vxlGD$vxl <- i
    vxlGD
  }
  class(GDs) <- "vbgamlss.predictions.GD"
  GDs
}


testGD <- function(nfitted, familyobj){
  ## Internal of predictGD ##

  dfun <- paste("d", familyobj$family, sep = "")
  pfun <- paste("p", familyobj$family, sep = "")
  lpar <- length(familyobj$names)
  if (is.null(nfitted$y))
    stop("the response variables is missing in the newdata")
  if (familyobj$family %in% .gamlss.bi.list) {
    if (NCOL(nfitted$y) == 1) {
      y1 <- nfitted$y
      bd <- nfitted$bd
    }
    else {
      bd <- nfitted$y[, 1] + nfitted$y[, 2]
      y1 <- nfitted$y[, 1]
    }
  }
  else {
    y1 <- nfitted$y
  }
  if (lpar == 1) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, bd = bd,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu)
    }
  }
  else if (lpar == 2) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$nmu, sigma = nfitted$sigma,
                   bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma)
    }
  }
  else if (lpar == 3) {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu)
    }
  }
  else {
    if (familyobj$family %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau, bd = bd,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau, bd = bd)
    }
    else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau)
    }
  }
  Vresid <- qNO(eval(ures))
  dev <- -2 * sum(eval(devi))
  out <- list()
  out <- list(TGD = dev, predictError = dev/length(nfitted$mu),
              resid = Vresid)
  return(out)
}


statGD <- function(GDs, k.penalty=NULL, deg.fre=1, return_all_GD=F) {
  nsub <- length(GDs[[1]]$resid)
  nvxl <- length(GDs)
  TGDs <- numeric(nvxl) * NA
  predictErrors <- numeric(nvxl) * NA
  resids <- matrix(data=NA, nrow=nsub, ncol=nvxl)

  for (i in seq_len(nvxl)) {
    TGDs[i] <- GDs[[i]]$TGD
    predictErrors[i] <- GDs[[i]]$predictError
    resids[, i] <- GDs[[i]]$resid
  }

  # penalize TGDs
  if (! is.null(k.penalty)){
    TGDs_pen <- TGDs - (k.penalty * deg.fre)
  }

  # statistics
  GD = describe_stats(TGDs, '')
  if (! is.null(k.penalty)){
    GDp = describe_stats(TGDs_pen, '')
  }
  PE = describe_stats(predictErrors, '')
  RS = describe_stats(resids, '')

  to_return = list(GD=GD, GDpen=GDp, predErr=PE, Resids=RS)

  if (return_all_GD) {to_return <- append(to_return, list(allGD=TGDs))}
  return(to_return)
}


describe_stats <- function(x, vname) {
  out <- c()
  for (fn in c('mean', 'sd', 'quantile', 'min', 'max')) {
    eval(parse(text=paste0(fn, vname,'<- ', fn,'(x)')))
    eval(parse(text=paste0('out <- c(out,', fn, vname, ' = ', fn, vname, ')')))
  }
  return(out)
}


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


getCVGD <- function(cvresults, term='mean') {
  # mean, sd,
  # quantile.0%, quantile.25%, quantile.50%, quantile.75%, quantile.100%,
  # min, max
  CVGD = 0
  for (fold in cvresults) {CVGD <- CVGD + fold$GD[[term]]}
  return(CVGD)
}


getCVGD.pen <- function(cvresults, term='quantile.50%') {
  CVGD = 0
  for (fold in cvresults) {CVGD <- CVGD + fold$GDpen[[term]]}
  return(CVGD)
}


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


akaike_weights <- function(v){
  # https://link.springer.com/content/pdf/10.3758/BF03206482.pdf
  # Dinga et al 2021
  v <- v - min(v)
  d <- exp(-0.5 * v)
  return(d/sum(d))
}






