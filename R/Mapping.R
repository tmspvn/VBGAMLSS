################################  Mapping  #####################################






#' Write coefficient maps to NIfTI files
#'
#' Convert per-voxel model coefficients from a fitted \code{vbgamlss} object
#' into coefficient images, one image per coefficient, and save them to disk.
#'
#' @param fittedobj A fitted object of class \code{"vbgamlss"} produced by voxel-wise fitting. Each element must contain a \code{$coefficients} vector with identical names and length across voxels.
#' @param mask An ANTs image (or path) defining the analysis mask that encodes image geometry for mapping matrices back to images.
#' @param filename Character prefix for output files. The function appends \code{"_par-<PAR>_coef-(<TERM>).nii.gz"} for each coefficient.
#' @param return_files Logical, if \code{TRUE} return the vector of file paths, otherwise return \code{invisible(NULL)}. Default \code{FALSE}.
#' @return If \code{return_files = TRUE}, a character vector of output paths. Files are written in NIfTI format.
#' @export
map_model_coefficients <- function(fittedobj, mask, filename, return_files=FALSE){
  if (class(fittedobj) != "vbgamlss") { stop("fittedobj must be of class vbgamlss.")}
  message('Warning, specific coefficients from special model terms cannot be map (e.g. pb(), s())')
  nvox <- length(fittedobj)
  first_mod_coefs <- unlist(fittedobj[[1]]$coefficients)
  name_coefs <- names(first_mod_coefs)
  ncoefs <- length(first_mod_coefs)
  coefs_mat <- matrix(nrow = ncoefs, ncol = nvox)
  for (i in 1:nvox) {
    coefs_mat[,i] <- unlist(fittedobj[[i]]$coefficients)
  }
  # convert mat to maps
  coef_maps_images <- matrixToImages(coefs_mat, antsImageRead(mask, 3))
  # save files
  fnames <- c()
  for (i in 1:length(coef_maps_images)) {
    # parse coef name
    parts <- strsplit(name_coefs[i], "\\.")[[1]]
    if (parts[2] != '(Intercept)') {parts[2] = paste0('(', parts[2], ')')}
    fname <- paste0(filename, '_par-', toupper(parts[1]), '_coef-', parts[2],'.nii.gz')
    # save
    antsImageWrite(coef_maps_images[[i]], fname)
    fnames[i] <- fname
  }
  if (return_files) {return(fnames)}
}




#' Write prediction maps to NIfTI files
#'
#' Convert per-voxel predicted parameter (μ,σ,ν,τ) values from a \code{vbgamlss.predictions} object into images and save them to disk.
#'
#' @param obj A \code{"vbgamlss.predictions"} object. Each voxel entry should  contain \code{$family}, optional \code{$vxl}, and one element per modelparameter (for example \code{mu}, \code{sigma}), each being a numeric vector over subjects with consistent ordering across voxels.
#' @param mask An ANTs image (or path) providing geometry, used to map matrices back to images.
#' @param filename Character prefix for output files. The function appends \code{"_subj-<ID>_fam-<FAMILY>_par-<PARAM>.nii.gz"}.
#' @param index Optional integer vector of subject indices to map. If \code{NULL} all subjects are processed.
#' @param return_files Logical, if \code{TRUE} return the vector of file paths written for the last parameter processed. Default \code{FALSE}.
#' @export
map_model_predictions <- function(obj, mask, filename, index=NULL,
                                  return_files=FALSE){
  if (class(obj) != "vbgamlss.predictions") { stop("obj must be of class vbgamlss.predictions")}
  if (! is.null(index)) {cat(paste0('Mapping predictions index: ', list(index)), fill=T)}

  # prepare useful info
  nvox <- length(obj)
  first_pred <- obj[[1]]
  family <- first_pred$family
  name_param <- names(first_pred)[! names(first_pred) %in% c('family', 'vxl')]
  nparam <- length(name_param)
  mask_img <- antsImageRead(mask, 3)

  # save a subset?
  if (is.null(index)) {
    nsubj <- length(first_pred[[name_param[1]]])
    subj = 1:nsubj # all
  } else {
    nsubj = length(subj)
    subj = index # subset
  }

  # process
  for (pname in name_param) {
    param_mat <- matrix(nrow=nsubj, ncol=nvox)
    for (ic in 1:nvox) {
      # voxel   #parameter  #subjects
      param_mat[, ic] <- obj[[ic]][[pname]][subj]
    }
    # convert mat to maps & save per subj
    param_maps_images <- matrixToImages(param_mat, mask_img)
    # save files
    fnames <- c()
    for (ip in 1:length(param_maps_images)) {
      fname <- paste0(filename,
                      '_subj-', subj[ip],
                      '_fam-', family,
                      '_par-' , toupper(pname),
                      '.nii.gz')
      antsImageWrite(param_maps_images[[ip]], fname)
      fnames[ip] <- fname
    }
  }
  if (return_files) {return(fnames)}
}




#' Write voxel-wise z-score maps to NIfTI files
#'
#' Convert per-voxel z-scores stored in a \code{"vbgamlss.zscores"} object into images and save one NIfTI file per subject.
#'
#' @param zscores A \code{"vbgamlss.zscores"} object. Each list element corresponds to a voxel and contains a numeric vector of z-scores over subjects (consistent ordering across voxels).
#' @param mask ANTs image (or path) providing geometry for mapping matrices back to images.
#' @param filename Character prefix for output files. The function appends \code{"_subj-<ID>.zscore.nii.gz"}.
#' @param index Optional integer vector of subject indices to export. If \code{NULL} (default), all subjects are exported.
#' @param return_files Logical, if \code{TRUE} return the vector of written file paths. Default \code{FALSE}.
#' @export
map_zscores <- function(zscores, mask, filename, index=NULL,
                        return_files=FALSE){
  if (class(zscores) != "vbgamlss.zscores") { stop("zscores must be of class vbgamlss.zscores")}
  if (! is.null(index)) {cat(paste0('Mapping predictions index: ', list(index)), fill=T)}

  # prepare usefull info
  nvox <- length(zscores)
  first_vxl <- zscores[[1]]
  # save a subset?
  if (is.null(index)) {
    nsubj <- length(first_vxl)
    subj = 1:nsubj # all
  } else {
    subj = index # subset
    nsubj = length(subj)
  }
  # process
  z_mat <- matrix(nrow=nsubj, ncol=nvox)
  for (ic in 1:nvox) {
    # voxel #subjects
    z_mat[, ic] <- zscores[[ic]][subj]
  }
  # convert mat to maps & save per subj
  z_maps_images <- matrixToImages(z_mat, antsImageRead(mask, 3))
  # save files
  fnames <- c()
  for (ip in 1:length(z_maps_images)) {
    fname <- paste0(filename,
                    '_subj-', index[ip],
                    '.zscore.nii.gz')
    antsImageWrite(z_maps_images[[ip]], fname)
    fnames[ip] <- fname
  }
  if (return_files) {return(fnames)}
}




#' Compute voxel-wise Z-scores from VBGAMLSS coefficients maps
#'
#' Given a subject image and predicted parameter maps (in nii.gz) for a GAMLSS family, compute zscore at each voxel inside a mask and return a z-score image.
#'
#' @param yimage A NIFTI image (or path readable by \code{antsImageRead}) containing the observed data to score.
#' @param mask An NIFTI image (or path) defining voxels to evaluate (nonzero entries are used).
#' @param family Character scalar naming a \pkg{gamlss.dist} family with an available CDF function \code{p<family>} (for example \code{"NO"}, \code{"SHASH"}). Default is \code{"SHASH"}.
#' @param mu_hat,sigma_hat,nu_hat,tau_hat NIFTI images (or paths) with the fitted parameter maps corresponding to the chosen family. Supply only the parameters required by the family.
#' @param num_cores Integer number of parallel workers. Defaults to \code{parallelly::availableCores()}.
#' @return An z-score image with voxel-wise z-scores.
#' @export
map_zscores_from_map <- function(yimage,
                                  mask,
                                  family='SHASH',
                                  mu_hat,
                                  sigma_hat=NULL,
                                  nu_hat=NULL,
                                  tau_hat=NULL,
                                  num_cores=NULL){
  if (is.null(num_cores)) {num_cores <- availableCores()}

  # make df
  frame = do.call(rbind,
                  list(y = images2matrix( yimage, mask),
                       mu = images2matrix( mu_hat, mask),
                       sigma = images2matrix( sigma_hat, mask),
                       nu = images2matrix( nu_hat, mask),
                       tau = images2matrix( tau_hat, mask)
                  ))
  frame <- as.data.frame(t(frame))

  # parallel function
  do.zscore <- function(i) {
    # get number of params
    lpar <- sum(names(frame) %in% c("mu", "sigma", "nu", "tau"))
    yval <- frame[i, 'y']
    pred <- frame[i, ]

    qfun <- paste("p", family, sep="")
    if (lpar == 1) {
      newcall <- call(qfun, yval, mu = pred$mu)
    } else if (lpar == 2) {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma)
    } else if (lpar == 3) {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu)
    } else {
      newcall <- call(qfun, yval, mu = pred$mu, sigma = pred$sigma, nu = pred$nu, tau = pred$tau)
    }
    cdf <- eval(newcall)
    rqres <- qnorm(cdf)
    return(rqres)
  }

  # compute z-scores
  subzs <- pbmclapply(1:dim(frame)[1],
                      do.zscore,
                      mc.cores=num_cores)
  zmap <- matrixToImages(matrix(unlist(subzs), nrow = 1), antsImageRead(mask,3))
  return(zmap[[1]])
}



