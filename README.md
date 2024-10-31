Fit GAMLSS models voxel-wise for normative modelling




Example:
```
set.seed(1)
img <- '~/subjects.nii.gz'  # 90x90x90x258
msk <- '~/mask.nii.gz'  # 90x90x90
nsubj <- 258
covs <- data.frame(x = 1:258, x1 = as.factor(rbinom(258, 1, 0.5)), x2 = rnorm(nsubj))  # 258x3
covs2 <- data.frame(x = rnorm(nsubj), x1 = rnorm(nsubj)*rnorm(nsubj), x2 = rnorm(nsubj))  # 258x3
models <- vbgamlss_chunks(image = img,
                   mask = msk,
                   g.formula = Y ~ pb(x) + x1 | x1,
                   g.family = NO,
                   num_cores = 20,
                   train.data = covs, # 258x3
                   debug = T)
save_model(models, '~/vbgamlss.model/fitted_model')
models_load <- load_model('~/vbgamlss.model/fitted_model.vbgamlss')
predictions <- predict.vbgamlss(models_load, newdata = covs2, ptype='response')
zscores <- zscore.vbgamlss(predictions, patients_image, mask)
```
























