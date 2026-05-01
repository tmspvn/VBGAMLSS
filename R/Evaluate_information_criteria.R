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
                              logdir = getwd(),
                              force_predict_recompute = FALSE,
                              ...) {

  # Ensure qs2 is loaded for fast I/O
  require(qs2, quietly = TRUE)

  if (is.null(num_cores)) {num_cores <- future::availableCores()}

  # Determine cache directory
  current_cache_dir <- get_cache_folders(logdir)

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
             cachedir = current_cache_dir,
             ...),
    skip = verbose)

  # 2. Predict on the same data (WITH CACHING via qs2)
  # Changed file extension to .qs
  prediction_cache_file <- file.path(current_cache_dir, "GDs_prediction_cache.qs")

  if (file.exists(prediction_cache_file) && !force_predict_recompute) {
    cat("Loading precomputed prediction metrics from cache (qs2)...\n")
    GDs <- qs_read(prediction_cache_file)
  } else {
    cat("Predicting and computing metrics...\n")
    GDs <- predict_metrics_EIC(model,
                           test_imageframe = imageframe,
                           newdata = data,
                           verbose = verbose,
                           segmentation = segmentation,
                           segmentation_target = segmentation_target)

    # Save the computed metrics to cache for future runs using qs2
    cat("Saving prediction metrics to cache (qs2)...\n")
    qs_save(GDs, file = prediction_cache_file)
  }

  # 3. Summarize statistics
  cat("Summarizing brain-wide statistics...\n")
  stats <- statGD_EIC(GDs,
                      n_obs = nrow(data),
                      return_all = return_all_metrics)

  rm(model, GDs)
  gc()

  return(stats)
}



# --------------------------------
# Streamlined Prediction (No CV/Disk States)
predict_metrics_EIC <- function(object, test_imageframe, newdata, verbose, segmentation, segmentation_target) {

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
      tryCatch({ testGD_EIC(nfitted[[i]], familyobj) }, error = function(e) NA)
    } else {
      NA
    }
  }

  return(GDs)
}

# --------------------------------
# Voxel-wise metrics (Unchanged core math, returns deviance, MAE, LL, CLL)
testGD_EIC <- function(nfit, familyobj) {
  require(gamlss.dist, quietly = TRUE)

  dfun <- paste("d", familyobj$family, sep = "")
  pfun <- paste("p", familyobj$family, sep = "")
  lpar <- length(familyobj$names)

  if (all(is.nan(nfit$y)) || all(is.na(nfit$y))) {
    return(list(TGD = NA, MAE = NA, LL = NA, CLL = NA, df = NA))
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

  df_val <- if(is.null(nfit$df)) NA else nfit$df

  list(TGD = dev, MAE = vxl_mae, LL = vxl_ll, CLL = vxl_cll, df = df_val)
}


# --------------------------------
# Brain-wide Summary (Now explicitly calculating AIC & BIC)
statGD_EIC <- function(GDs, n_obs, return_all = FALSE) {
  # Note: is.na on a list might not work as intended for missing elements.
  # We count missfits by checking if the element is NOT a list.
  missfits <- sum(!vapply(GDs, is.list, logical(1)))
  nvxl <- length(GDs)

  TGDs <- numeric(nvxl) * NA
  MAEs <- numeric(nvxl) * NA
  LLs  <- numeric(nvxl) * NA
  CLLs <- numeric(nvxl) * NA
  dfs  <- numeric(nvxl) * NA

  for (i in seq_len(nvxl)) {
    if (is.list(GDs[[i]])) {
      TGDs[i] <- if (is.null(GDs[[i]]$TGD)) NA else GDs[[i]]$TGD
      MAEs[i] <- if (is.null(GDs[[i]]$MAE)) NA else GDs[[i]]$MAE
      LLs[i]  <- if (is.null(GDs[[i]]$LL)) NA else GDs[[i]]$LL
      CLLs[i] <- if (is.null(GDs[[i]]$CLL)) NA else GDs[[i]]$CLL
      dfs[i]  <- if (is.null(GDs[[i]]$df)) NA else GDs[[i]]$df
    }
  }

  # Calculate AIC and BIC using Deviance (TGD) + Penalties
  AICs <- TGDs + (2 * dfs)
  BICs <- TGDs + (log(n_obs) * dfs)

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


# --------------------------------
get_cache_folders <- function(parent_dir) {
  if (!dir.exists(parent_dir)) {
    stop("The specified parent directory does not exist.")
  }

  all_dirs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  cache_dirs <- all_dirs[grepl("^\\.vbgamlss\\.cache", basename(all_dirs))]

  if (length(cache_dirs) == 0) {
    return(parent_dir)
  } else {
    # FIX: Return only the first element to prevent passing a vector of paths
    return(cache_dirs[1])
  }
}












