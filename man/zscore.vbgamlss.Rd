% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Support.R
\name{zscore.vbgamlss}
\alias{zscore.vbgamlss}
\title{Compute Z-scores for vbgamlss predictions given Y voxel data and image mask.}
\usage{
zscore.vbgamlss(predictions, yimage, ymask, num_cores = NULL)
}
\arguments{
\item{predictions}{Predictions from vbgamlss model.}

\item{yimage}{R object or filename containing image data.}

\item{ymask}{R object or filename containing mask data.}

\item{num_cores}{Number of CPU cores to use for parallel processing. Defaults to all available cores.}
}
\value{
A structure containing Z-scores.
}
\description{
Compute Z-scores for vbgamlss predictions given Y voxel data and image mask.
}
