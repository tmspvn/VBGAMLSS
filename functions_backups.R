#' Run a Generalized Additive Model on all voxels of a NIfTI image within a mask. 
#'
#' This function is able to run a Generalized Additive Model (GAM) using the mgcv package. 
#' The analysis will run in all voxels in in the mask and will return the model fit for each voxel.
#' 
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will call mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. Must be a binary mask
#' @param fourdOut To be passed to mergeNifti, This is the path and file name without the suffix to save the fourd file. Default (NULL) means script won't write out 4D image.
#' @param formula Must be a formula passed to gam()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param ncores Number of cores to use
#' @param ... Additional arguments passed to gam()
#' 
#' @return List of models fitted to each voxel over the masked images passed to function.
#' @keywords internal
#' @export
#' @examples

library(foreach)
library(doSNOW)
library(gamlss)
library(gamlss2)
library(pbmcapply)
library(voxel)

idivix <- function(n, chunkSize) {
  # from: https://stackoverflow.com/a/37838711
  # function to chunck but keep order
  i <- 1
  # Create an iterator object using the 'idiv' function with the specified chunk size
  it <- idiv(n, chunkSize=chunkSize)
  # Define a function named 'nextEl' to retrieve the next element from the iterator
  nextEl <- function() {
    # Retrieve the next chunk size from the iterator
    m <- nextElem(it)  # may throw 'StopIterator'
    # Create a list containing the current index 'i' and the retrieved chunk size 'm'
    value <- list(i=i, m=m)
    # Increment the index 'i' by the chunk size 'm'
    i <<- i + m
    # Return the list containing the index and chunk size
    value
  }
  # Create an object 'obj' containing the 'nextEl' function
  obj <- list(nextElem=nextEl)
  # Assign the class attributes 'abstractiter' and 'iter' to the object 'obj'
  class(obj) <- c('abstractiter', 'iter')
  obj
}


vbgamlss <- function(image, mask, multi.tissue=NULL, 
                     mu.formula, train.data, 
                     ncores=1, chunks.size=75,
                     fmodel,
                     ...) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(mu.formula)) { stop("formula is missing")}
  if (missing(train.data)) { stop("subjData is missing")}
  
  if (class(mu.formula) != "character") { stop("mu.formula class must be character")}
  
  if (class(image) == "character") {image <- oro.nifti::readNIfTI(fname=image)} 
  
  if (class(mask) == "character") {mask <- oro.nifti::readNIfTI(fname=mask)}
  
  if (class(multi.tissue) == "character") {multi.tissue <- oro.nifti::readNIfTI(fname=multi.tissue)} 
  
  # 4d to 2D and set a global variable: subjects (4th dim) x voxels (columns)
  voxeldata <<- ts2matrix(image, mask)
  if (!is.null(multi.tissue)){multi.tissue <- ts2matrix(multi.tissue, mask)}
  rm(image)
  rm(mask)
  gc()
  
  # Coerce update on left hand mu.formula for parallel fitting
  mu.fo <- update.formula(mu.formula, 'Y ~ .')
  
  print('Foreach sync chunking')
  iterations=ncol(voxeldata)
  #cl <- makeCluster(1, outfile="") # for debug
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  pb <- progressBar(max = iterations) #/chunks.size) # tick per chunk
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # Foreach call
  ctrl = gamlss2_control(light = TRUE, trace=FALSE)
  # no chucking
  models <- foreach(vxlcol = voxeldata,
                    .packages = 'gamlss2',
                    .combine='c',
                    .options.snow = opts) %dofuture% {
                      # fit specific voxel
                      vxl_train_data <- train.data
                      vxl_train_data$Y <- vxlcol
                      # if multi tissue add
                      if (!is.null(multi.tissue)){
                        vxl_train_data$tissue <- multi.tissue[,i]}
                      # gamlss
                      g <- gamlss2::gamlss2(formula=mu.fo,
                                            data=vxl_train_data,
                                            control=ctrl,
                                            ...)
                      g$vxl<-names(vxlcol)
                      wrap <- list(g)
                      wrap
                    }
  
  # Stop parallel processing
  gc()
  close(pb)
  stopCluster(cl)
  return(models)
  
}


















