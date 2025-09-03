# Fit GAMLSS Models Voxel or Vertex-Wise for Normative Modelling

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)](https://cran.r-project.org/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-green.svg)](LICENSE)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](#)
[![Issues](https://img.shields.io/github/issues/tmspvn/VBGAMLSS)](https://github.com/tmspvn/VBGAMLSS/issues)

**VBGAMLSS** fits **Generalized Additive Models for Location, Scale and Shape** voxel-wise or vertex-wise, designed for normative modelling in neuroimaging.

---

## 📦 Requirements

```
ANTsR
doFuture
gamlss
gamlss2
itertools
pbmcapply
progressr
tibble
```

---

## 🔧 Installation

### Install `gamlss2` and `ANTsR`

```r
install.packages(
  "gamlss2", 
  repos = c("https://gamlss-dev.R-universe.dev", "https://cloud.R-project.org")
)
devtools::install_github("ANTsX/ANTsR")
```

### Install `VBGAMLSS`

```r
devtools::install_github("tmspvn/VBGAMLSS", dependencies = TRUE)
```

---

## 🚀 Quick Start

```r
library(VBGAMLSS)

# Paths to example data
img <- "~/subjects.nii.gz"  # 90x90x90x258
msk <- "~/mask.nii.gz"      # 90x90x90
nsubj <- 258

# Covariates
covs <- data.frame(
  x  = 1:nsubj,
  x1 = as.factor(rbinom(nsubj, 1, 0.5)),
)

covs_patients <- data.frame(
  x  = rnorm(nsubj),
  x1 = rnorm(nsubj) * rnorm(nsubj),
)

# Convert to 2D subject × voxel/vertex
imageframe <- images2matrix(img, msk)

# Fit voxel-wise GAMLSS models
models <- vbgamlss(
  imageframe,
  g.formula = Y ~ pb(x) + x1 | x1,
  g.family  = NO,
  num_cores = 20,
  train.data = covs,
  debug = TRUE
)

# Save / load
save_model(models, "~/vbgamlss.model/fitted_model")
models_loaded <- load_model("~/vbgamlss.model/fitted_model.vbgamlss")

# Predict & compute Z-scores
predictions <- predict.vbgamlss(models_loaded, newdata = covs_patients)
zscores <- zscore.vbgamlss(predictions, patients_imageframe)
```

---

## ⚠ Known Issues

* `vbgamlss.model_selection` and `vbgamlss.cv` are **under development** and not fully tested.
* Segmentation handling is experimental.

---

