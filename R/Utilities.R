################################  Utilities  ###################################








#' Estimate number of chunks based on object size
#'
#' @param image_list list containing image paths.
#' @param mask mask image paths.
#' @return Returns matrix of dimension (numImages, numVoxelsInMask)
#' @export
images2matrix <- function(image_list, mask) {

  if (missing(image_list)) { stop("image_list is missing")}
  if (missing(mask)) { stop("mask is missing")}
  # load mask
  mask_ <- antsImageRead(mask, 3)
  mask_arr <- as.array(mask_) > 0
  numVoxels <- length(which(mask_arr))
  # check if list of images of image
  if (class(image_list) == "list"){
    cat('Converting from list of images paths', fill=T)
    numImages <- length(image_list)
    dataMatrix <- matrix(nrow = numImages, ncol = numVoxels)
    for (i in 1:length(image_list))
    #{dataMatrix[i, ] <- as.numeric(image_list[[i]], mask = mask_arr)}
    {dataMatrix[i, ] <- as.numeric(antsImageRead(image_list[[i]], 3))[mask_arr]}

  } else {
    if (class(image_list) != "character") {
      stop("image_list is not a list of nifties or a path to a 4D image")}
    cat('Converting from 4D image path', fill=T)
    large_image <- antsImageRead(image_list, 4)
    numImages <- dim(large_image)[4]
    dataMatrix <- matrix(nrow = numImages, ncol = numVoxels)
    for (i in 1:numImages)
    {dataMatrix[i, ] <- as.numeric(as.antsImage(drop(large_image[,,,i])), mask = mask_arr)}

  }
  return(as.data.frame(dataMatrix))
}



#' Check input data
#'
#' @param image path of the image to load and check. RDS, nifti or dataframe
#' @param mask  path to mask.
#' @return image as a dataframe.
#' @export
load_input_image <- function(image, mask=NULL){
  if (missing(image)) { stop("image is missing")}
  # check image input format
  if (is.character(image)){
    # check if file exist
    if (!file.exists(image)) {stop("Image file does not exist.")}

    # try load rds?
    is_rds <- tryCatch({
              voxeldata <- readRDS(image)
              cat("loading RDS file, mask is ignored.")
              TRUE
            }, error = function(e){F})

    # if not rds try load nifti?
    if (!is_rds){
      if (is.null(mask)) { stop("mask must be specified if image is a path to .nii")}
      tryCatch({
        voxeldata <- images2matrix(image, mask)
      }, error = function(e){
        cat('Error: ', e$call, '\n')
        cat('Error: ', e$message, '\n')
        stop("Image input must be a *path* to a RDS, nifti *or* a data frame object.")
      })
    }

    # if not a character, try load data frame
  } else {
    if (is.data.frame(image)){
      voxeldata <- image
      warning("Data.frame passed, mask is ignored.")
      } else{
      stop("Image input must be a *path* to a RDS, nifti *or* a data frame object.")
    }
  }
  return(voxeldata)
}


########### NOT EXPORTED ###########

estimate_nchunks <- function(object, from_files=F, chunk_max_Mb=256) {
  # compute chunk size of max 256 mb per job
  if (! from_files){
    memory_size_mb <- object.size(object) / (1024^2)
    Nchunks <- ceiling(memory_size_mb / chunk_max_Mb)
  } else {
    memory_size_mb <- file.info(object)$size / (1024^2)
    Nchunks <- ceiling(memory_size_mb / chunk_max_Mb)
  }
  return(as.numeric(Nchunks))
}



get_subsample_indices <- function(len, sampling_type, factor) {
if (!is.numeric(factor) && factor>0 && factor<=1) {
      stop("Error: factor (or subsample) must be numeric between (0, 1]")
    }

  if (!sampling_type %in% c("random", "regular")) {
    stop("Sampling type must be either 'random' or 'regular'.")
  }

  # Calculate the number of elements to sample
  n <- ceiling(len * factor)

  cat(paste0('Subsampling by a factor of ', factor, ' (', n,'/', len,' vxl)',
             ', with strategy: ', sampling_type, '\n\n'), fill=T)

  if (sampling_type == "regular") {
    # Calculate the step size for regular sampling
    step_size <- floor(len / n)
    # Generate a sequence of indices to select from the vector
    indices <- seq(1, len, by = step_size)
    # Ensure we only select `n` indices
    indices <- indices[1:n]
  } else if (sampling_type == "random") {
    # Randomly sample `n` indices from the vector
    indices <- sample(1:len, n)
  }

  return(indices)
}





combine_formulas_gamlss2 <- function(mu_list, sigma_list,
                                     nu_list, tau_list,
                                     families=NULL) {
  # warning
  dowarn <- (length(mu_list) * length(sigma_list) * length(nu_list) * length(tau_list) * length(families))>200
  if (dowarn){warning('Number of formula to test in >200!')}
  # check if there's a ~1 else force it
  add_intercept_if_missing <- function(lst) {
    if ( ! '1' %in% lst) {
        lst[length(lst)[1]+1]  <- '1'
    }
    return(lst)
  }
  # add intercept?
  sigma_list <- add_intercept_if_missing(sigma_list)
  nu_list <- add_intercept_if_missing(nu_list)
  tau_list <- add_intercept_if_missing(tau_list)
  # combine
  combined_formulas <- c()
  i <- 0
  for (formula1 in mu_list) {
    for (formula2 in sigma_list) {
      for (formula3 in nu_list) {
        for (formula4 in tau_list) {
          combined_formula <- paste(formula1, " | ", formula2, " | ",
                                    formula3, " | ", formula4)
          combined_formulas <- c(combined_formulas, combined_formula)
          i <- i+1
      }
     }
    }
  }
  # map families if requested
  combined_formulas_with_fam <- c()
  if (!is.null(families)){
    for (fam in families){
        fam_form <- paste0(fam, ' :: ', combined_formulas)
        combined_formulas_with_fam <- c(combined_formulas_with_fam, fam_form)
    }
    combined_formulas <- combined_formulas_with_fam
  }
  return(combined_formulas)
}


quite <- function(x, skip=F) {
  if (! skip){
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  } else {force(x)}
}


rand_names <- function(num_names, l=10){
  return(replicate(num_names, paste0(sample(c(letters, LETTERS), l, replace = TRUE), collapse = "")))
}


check_formula_LHS <- function(formula) {
  formula <- as.formula(formula)
  terms_formula <- terms(formula)
  lhs <- attr(terms_formula, "variables")[[2]]

  if (deparse(lhs) != "Y") {
    stop("Error: The left term of the formula must be 'Y'.")
  }
}



TRY <- function(expr, logfile=NULL, save.env.and.stop=F){
  res <-tryCatch(expr={expr},
                  error = function(e) {
                   if (!is.null(logfile)){
                     cat(paste('ERROR:\n\n', e$message), "\n\n",
                         file = logfile, append = TRUE)
                     if (save.env.and.stop) {
                       save(list = ls(all.names = TRUE),
                            file = paste0(logfile, '.ERROR.local.enviroment'))
                       stop('Stopping on first error,
                            *save.env.and.stop* is set to TRUE')
                         }
                     }
                     g <- NA # missfit
                   },
                  warning = function(w) {
                    if (!is.null(logfile)){
                      # Save the error message to a file
                      cat(paste('WARN:\n\n', w$message), "\n\n",
                          file = logfile, append = TRUE)
                    }
                    expr
                  }
  )
  return(res)}




############################  R Package fucntions  ##############################

.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("ANTsR", quietly = TRUE)) {
    missing_functions <- paste(c("antsImageRead", "antsImageWrite"), collapse = ", ")
    packageStartupMessage(
      paste0(
        "*** WARNING: The 'ANTsR' package is not installed. ***\n",
        "  The following functions will not work: ", missing_functions, "\n",
        "  Install it with 'devtools::install_github(\"ANTsX/ANTsR\")' if you need to work with voxel based data."
      ),
      appendLF = TRUE # Ensure a newline after the message
    )
  }
}











