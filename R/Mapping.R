################################  Mapping  #####################################







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

#' @export
map_model_predictions <- function(obj, mask, filename, index=NULL,
                                  return_files=FALSE){
  if (class(obj) != "vbgamlss.predictions") { stop("obj must be of class vbgamlss.predictions")}
  if (! is.null(index)) {cat(paste0('Mapping predictions index: ', list(index)), fill=T)}
  # prepare usefull info
  nvox <- length(obj)
  first_pred <- obj[[1]]
  family <- first_pred$family
  name_param <- names(first_pred)[! names(first_pred) %in% c('family', 'vxl')]
  nparam <- length(name_param)
  # save a subset?
  if (is.null(index)) {
    nsubj <- length(first_pred[[name_param[1]]])
    subj = 1:nsubj # all
  } else {
    subj = index # subset
    nsubj = length(subj)
  }
  # process
  for (pname in name_param) {
    param_mat <- matrix(nrow=nsubj, ncol=nvox)
    for (ic in 1:nvox) {
      # voxel   #parameter  #subjects
      param_mat[, ic] <- obj[[ic]][[pname]][subj]
    }
    # convert mat to maps & save per subj
    param_maps_images <- matrixToImages(param_mat, antsImageRead(mask, 3))
    # save files
    fnames <- c()
    for (ip in 1:length(param_maps_images)) {
      fname <- paste0(filename,
                      '_subj-', index[ip],
                      '_fam-', family,
                      '_par-' , toupper(pname),
                      '.nii.gz')
      antsImageWrite(param_maps_images[[ip]], fname)
      fnames[ip] <- fname
    }
  }
  if (return_files) {return(fnames)}
}

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
