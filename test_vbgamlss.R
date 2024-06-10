source('~/PycharmProjects/NormativeModellingPsychosis/gamlss/VBGAMLSS/VBGAMLSS.R')

### test ### https://cran.r-project.org/web/packages/voxel/voxel.pdf
image <- '/home/localadmin/PycharmProjects/NormativeModellingPsychosis/gamlss/tests/funneled.nii.gz'
mask <- '/home/localadmin/PycharmProjects/NormativeModellingPsychosis/gamlss/tests/vol.0_mask.nii.gz'
set.seed(1)
nsubj <- 258
covs <- data.frame(x = rnorm(nsubj), x1 = rnorm(nsubj)*runif(nsubj), x2 = rnorm(nsubj))
covs2 <- data.frame(x = rnorm(nsubj), x1 = rnorm(nsubj)*rnorm(nsubj), x2 = rnorm(nsubj))
fm1 <- "~ x + x1 | x + x2 "

models <- vbgamlss_chunks(image=image, 
                   mask=mask, 
                   g.formula=fm1, 
                   g.family=NO,
                   train.data=covs)
save_model(models, '/home/localadmin/PycharmProjects/NormativeModellingPsychosis/gamlss/tests/vbgamlss.model/vbgamlss')
models_load <- load_model_chunks('/home/localadmin/PycharmProjects/NormativeModellingPsychosis/gamlss/tests/vbgamlss.model/vbgamlss')
predictions <- predict.vbgamlss(models_load, newdata = covs2)
zscores <- zscore.vbgamlss(predictions, image, mask)

# not saving properly because of the name()
map_model_coefficients(models, mask, filename='/home/localadmin/PycharmProjects/NormativeModellingPsychosis/gamlss/tests/pmap')


# test antsr
#img <- antsImageRead(image)
#msk <- antsImageRead(mask)
#mat<-imageListToMatrix(img, msk)



# data("abdom", package = "gamlss.data")
# f <- y ~ x | x
# g1 <- gamlss(y ~ x, ~x, data = abdom, family = BCT)
# g2 <- gamlss2(f, data = abdom, family = BCT, light=T)
# g3 <- g2
# family_ <- g3$family
# g3$family = 'BCT'
# g3$control=NULL
# 
# print(object.size(g1) / (1024 * 1024))
# print(object.size(g2) / (1024 * 1024))
# print(object.size(g3))
# obj = g2
# for (attr_name in names(obj)) {
#   attr <- obj[[attr_name]]
#   size_mb <- object.size(attr) / (1024 * 1024)
#   cat(paste("Size of", attr_name, ":", round(size_mb, 2), "MB\n"))
# }
# 
# 
# g3$family = family_
# g3$family = gamlss2:::complete_family(BCT)
# predict(g3)



df <- data.frame(read.csv(file = '~/PycharmProjects/NormativeModellingPsychosis/datasets/normod_df.csv'))
df$sex <- as.factor(df$sex)
df$scanner <- as.factor(df$scanner)
df$group <- as.factor(df$group)

f <- 

## Estimate model.
b <- gamlss2(volume ~ pb(ageY) | sex | group, 
             data = df, 
             family = SHASH)

## Model summary.
summary(b)
























