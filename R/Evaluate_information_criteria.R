#############################  Model Evaluation  ###############################

#' @export
vbgamlss.evaluate <- function(imageframe,
                              g.formula,
                              data,
                              g.family = NO,
                              segmentation = NULL,
                              segmentation_target = NULL,
                              num_cores = NULL,
                              chunk_max_mb = 64,
                              force_constraints = c(1e-8, +Inf),
                              verbose = FALSE,
                              debug = TRUE,
                              return_all_metrics = TRUE,
                              ...) {

  cat("Fitting model on full dataset...\n")

  # 1. Fit the model once
  model <- quite(
    vbgamlss(imageframe = imageframe,
             g.formula = g.formula,
             train.data = data,
             g.family = g.family,
             segmentation = segmentation,
             segmentation_target = segmentation_target,
             force_constraints = force_constraints,
             num_cores = num_cores,
             chunk_max_mb = chunk_max_mb,
             debug = debug,
             ...),
    skip = verbose)

  # 2. Predict on the same data
  cat("Predicting and computing metrics...\n")
  GDs <- predict_metrics(model,
                         test_imageframe = imageframe,
                         newdata = data,
                         verbose = verbose,
                         segmentation = segmentation,
                         segmentation_target = segmentation_target)

  # 3. Summarize statistics
  cat("Summarizing brain-wide statistics...\n")
  all_dfs <- unlist(pbmcapply::pbmclapply(model, function(m_ser) {
                    # Deserializing each voxel fit
                    m_obj <- qs2::qs_read(m_ser)
                    return(m_obj$df)
                  }, mc.cores = num_cores))

  stats <- statGD_EIC(GDs,
                  deg.fre = all_dfs,
                  n_obs = nrow(data),
                  return_all = return_all_metrics)

  rm(model, GDs)
  gc()

  return(stats)
}

# --------------------------------
# Streamlined Prediction (No CV/Disk States)
predict_metrics <- function(object, test_imageframe, newdata, verbose, segmentation, segmentation_target) {

  familyobj <- restore_family(object[[1]])$family

  # Predict Parameters
  cat("\033[34m")
  nfitted <- quite(
    predict.vbgamlss(object, newdata = newdata, ptype = 'parameter',
                     segmentation = segmentation, segmentation_target = segmentation_target),
    skip = verbose)

  # Predict Response
  resp <- quite(
    predict.vbgamlss(object, newdata = newdata, ptype = 'response',
                     segmentation = segmentation, segmentation_target = segmentation_target),
    skip = verbose)
  cat("\033[0m")

  # Merge observed y and yhat
  nsub <- dim(newdata)[1]
  for(i in seq_along(nfitted)) {
    if (!"try-error" %in% class(nfitted[[i]]) && !all(is.na(nfitted[[i]]))) {
      nfitted[[i]][['y']] <- test_imageframe[, i]
      nfitted[[i]][['yhat']] <- as.numeric(unlist(resp[[i]])[1:nsub])
    } else {
      nfitted[[i]] <- NA
    }
  }

  # Compute voxel-wise Global Deviance and other metrics
  plan(strategy = "future::cluster", workers = availableCores())
  GDs <- foreach(i = seq_along(nfitted)) %dofuture% {
    if (!is.na(nfitted)[i]) {
      tryCatch({ testGD(nfitted[[i]], familyobj) }, error = function(e) NA)
    } else {
      NA
    }
  }

  return(GDs)
}

# --------------------------------
# Voxel-wise metrics (Unchanged core math, returns deviance, MAE, LL, CLL)
testGD <- function(nfit, familyobj) {
  require(gamlss.dist, quietly = TRUE)

  dfun <- paste("d", familyobj$family, sep = "")
  pfun <- paste("p", familyobj$family, sep = "")
  lpar <- length(familyobj$names)

  if (all(is.nan(nfit$y)) || all(is.na(nfit$y))) {
    return(list(TGD = NA, MAE = NA, LL = NA, CLL = NA))
  }

  yisnan <- is.nan(nfit$y)
  for (mtd in c('y', 'yhat', 'mu', 'sigma', 'nu', 'tau')[1:(lpar+1)]) {
    eval(parse(text=paste0('nfit$',mtd,' <- nfit$',mtd,'[! yisnan]')))
  }

  y1 <- nfit$y

  if (lpar == 1) {
    devi <- call(dfun, x = y1, mu = nfit$mu, log = TRUE)
    ures <- call(pfun, q = y1, mu = nfit$mu)
  } else if (lpar == 2) {
    devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, log = TRUE)
    ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma)
  } else if (lpar == 3) {
    devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, log = TRUE)
    ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu)
  } else {
    devi <- call(dfun, x = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, tau = nfit$tau, log = TRUE)
    ures <- call(pfun, q = y1, mu = nfit$mu, sigma = nfit$sigma, nu = nfit$nu, tau = nfit$tau)
  }

  dev <- -2 * sum(eval(devi))
  vxl_mae <- mean(abs(y1 - nfit$yhat), na.rm = TRUE)

  ll_obs_raw <- -eval(devi)
  cdf_vals <- eval(ures)

  vxl_ll <- mean(ll_obs_raw, na.rm = TRUE)

  ll_obs_censored <- ll_obs_raw
  ll_obs_censored[cdf_vals < 0.01 | cdf_vals > 0.99] <- -log(0.02)
  vxl_cll <- mean(ll_obs_censored, na.rm = TRUE)

  list(TGD = dev, MAE = vxl_mae, LL = vxl_ll, CLL = vxl_cll)
}

# --------------------------------
# Brain-wide Summary (Now explicitly calculating AIC & BIC)
statGD_EIC <- function(GDs, deg.fre, n_obs, return_all = FALSE) {
  missfits <- sum(is.na(GDs))
  nvxl <- length(GDs)

  TGDs <- numeric(nvxl) * NA
  MAEs <- numeric(nvxl) * NA
  LLs  <- numeric(nvxl) * NA
  CLLs <- numeric(nvxl) * NA

  for (i in seq_len(nvxl)) {
    if (!is.na(GDs)[i]) {
      TGDs[i] <- GDs[[i]]$TGD
      MAEs[i] <- GDs[[i]]$MAE
      LLs[i]  <- GDs[[i]]$LL
      CLLs[i] <- GDs[[i]]$CLL
    }
  }

  # Calculate AIC and BIC using Deviance (TGD) + Penalties
  AICs <- TGDs + (2 * deg.fre)
  BICs <- TGDs + (log(n_obs) * deg.fre)

  to_return <- list(
    GD  = describe_stats(TGDs, ''),
    MAE = describe_stats(MAEs, ''),
    LL  = describe_stats(LLs, ''),
    CLL = describe_stats(CLLs, ''),
    AIC = describe_stats(AICs, ''),
    BIC = describe_stats(BICs, ''),
    missFitsPerc = missfits / nvxl
  )

  if (return_all) {
    to_return <- c(to_return, list(
      allGD = TGDs, allMAE = MAEs, allLL = LLs, allCLL = CLLs, allAIC = AICs, allBIC = BICs
    ))
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
