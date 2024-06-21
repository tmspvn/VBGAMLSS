library(caret)
library(foreach)
library(future)
library(splitTools)

gamlss.cv <- function(image, mask, g.formula, train.data, g.family = NO,
                      segmentation = NULL, num_cores = NULL, chunk_max_mb = 64,
                      n_folds = 5, ...) {
  if (missing(image)) { stop("image is missing") }
  if (missing(mask)) { stop("mask is missing") }
  if (missing(g.formula)) { stop("formula is missing") }
  if (missing(train.data)) { stop("subjData is missing") }
  
  if (class(g.formula) != "character") { stop("g.formula class must be character") }
  if (is.null(num_cores)) { num_cores <- availableCores() }
  
  # Create stratified folds
  train.data$folds <- stratCVfolds(train.data$Y, k = n_folds, list = FALSE)
  
  # Prepare storage for results
  results <- list()
  
  # Perform stratified cross-validation
  for (fold in 1:n_folds) {
    cat(paste0("Processing fold ", fold, " of ", n_folds, "\n"))
    
    # Split data into training and validation sets
    training_fold = train.data$folds == fold
    test_indices <- which(train.data$folds == fold)
    test_fold_data <- train.data[test_indices, ]
    
    # Fit the model on the training fold
    model <- vbgamlss_chunks(image, 
                             mask, 
                             g.formula, 
                             train.data, 
                             g.family=NO,
                             segmentation, 
                             num_cores,
                             chunk_max_mb=128,
                             afold=training_fold)
    # predict new fold
    predictions <- predict.vbgamlss(model, newdata=test_fold_data, num_cores)
    
    
    # Store the model and validation set
    results[[fold]] <- list(
      model = model,
      test_data = test_fold_data
    )
  }
  
  return(results)
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
