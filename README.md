Fit GAMLSS models voxel-wise for normative modelling

# requirements
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
# install gamlss2
```
install.packages("gamlss2", repos = c("https://gamlss-dev.R-universe.dev", "https://cloud.R-project.org"))
```

Example, package is work in progress, please source individual scripts:
```
set.seed(1)
img <- '~/subjects.nii.gz'  # 90x90x90x258
msk <- '~/mask.nii.gz'  # 90x90x90
nsubj <- 258
covs <- data.frame(x = 1:258, x1 = as.factor(rbinom(258, 1, 0.5)), x2 = rnorm(nsubj))  # 258x3
covs_patients <- data.frame(x = rnorm(nsubj), x1 = rnorm(nsubj)*rnorm(nsubj), x2 = rnorm(nsubj))  # 258x3
imageframe <- images2matrix(img, msk)
models <- vbgamlss(imageframe, # Data.Frame subject x voxels/vertices
                   g.formula = Y ~ pb(x) + x1 | x1,
                   g.family = NO,
                   num_cores = 20,
                   train.data = covs, # 258x3
                   debug = T)
save_model(models, '~/vbgamlss.model/fitted_model')
models_load <- load_model('~/vbgamlss.model/fitted_model.vbgamlss')
predictions <- predict.vbgamlss(models_load, newdata = covs_patients, ptype='response')
zscores <- predict.vbgamlss(predictions, patients_image, mask)
```



To Be Fixed:
predict.vbgamlss, predict.vbgamlss still uses mask within their code. It needs to be fixed





















