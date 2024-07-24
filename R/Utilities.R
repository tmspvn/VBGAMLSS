################################  Utilities  ###################################








#' Estimate number of chunks based on object size
#'
#' @param object R object for chunk size estimation.
#' @return Estimated number of chunks with max 256MB per job.
#' @examples estimate_nchunks(data)
#' @export
estimate_nchunks <- function(object, from_files=F, chunk_max_Mb=256) {
  # compute chunk size of max 256 mb per job
  if (! from_files){
    memory_size_mb <- object.size(object) / (1024^2)
    Nchunks <- ceiling(memory_size_mb / chunk_max_Mb)
  } else {
    memory_size_mb <- file.info(object)$size / (1024^2)
    Nchunks <- ceiling(memory_size_mb / chunk_max_Mb)
  }
  return(Nchunks)
}



#' @export
images2matrix <- function(image_list, mask) {
  # imageList is a list containing images.  Mask is a mask image Returns matrix of
  # dimension (numImages, numVoxelsInMask)
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
    for (i in 1:length(imageList))
    {dataMatrix[i, ] <- as.numeric(imageList[[i]], mask = mask_arr)}

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




#' @export
matrix2image <- function(matrix, mask){
  if (missing(matix)) { stop("matrix is missing")}
  if (missing(mask)) { stop("mask is missing")}
  return()
}



#' @export
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




#' @export
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

rand_names <- function(num_names){
  return(replicate(num_names, paste0(sample(c(letters, LETTERS), 10, replace = TRUE), collapse = "")))
}

check_formula_LHS <- function(formula) {
  formula <- as.formula(formula)
  terms_formula <- terms(formula)
  lhs <- attr(terms_formula, "variables")[[2]]

  if (deparse(lhs) != "Y") {
    stop("Error: The left term of the formula must be 'Y'.")
  }
}
