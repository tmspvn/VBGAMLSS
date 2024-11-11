# Fit GAMLSS models voxel-wise for normative modelling

### requirements
```
ANTsR,
doFuture,
gamlss,
gamlss2,
itertools,
pbmcapply,
progressr,
tibble
```
#### install gamlss2 & ANTsR
```
install.packages("gamlss2", repos = c("https://gamlss-dev.R-universe.dev", "https://cloud.R-project.org"))
devtools::install_github('ANTsX/ANTsR')
```

# Example 
### package is work in progress, please source individual scripts:
```
set.seed(1)
img <- '~/subjects.nii.gz'  # 90x90x90x258
msk <- '~/mask.nii.gz'  # 90x90x90
nsubj <- 258
covs <- data.frame(x = 1:258, x1 = as.factor(rbinom(258, 1, 0.5)), x2 = rnorm(nsubj))  # 258x3
covs_patients <- data.frame(x = rnorm(nsubj), x1 = rnorm(nsubj)*rnorm(nsubj), x2 = rnorm(nsubj))  # 258x3

# Pre-mask the voxel image and covert to 2d dataframe
imageframe <- images2matrix(img, msk)

# Fit models voxel/vertex-wise
models <- vbgamlss(imageframe, # Data.Frame subject x voxels/vertices
                   g.formula = Y ~ pb(x) + x1 | x1,
                   g.family = NO,
                   num_cores = 20,
                   train.data = covs, # 258x3
                   debug = T)

# Save models
save_model(models, '~/vbgamlss.model/fitted_model')

# Load models
models_loaded <- load_model('~/vbgamlss.model/fitted_model.vbgamlss')

# Predict models response given new data
predictions <- predict.vbgamlss(models_loaded, newdata = covs_patients, ptype='response')

# Compute Z-scores given new data
zscores <- zscore.vbgamlss(predictions, patients_imageframe)

```

## Bugs
```vbgamlss.model_selection``` & ```vbgamlss.cv``` are in development and the current implementation is not working.
















