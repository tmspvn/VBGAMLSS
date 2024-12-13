% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{estimate_nchunks}
\alias{estimate_nchunks}
\title{Estimate number of chunks based on object size}
\usage{
estimate_nchunks(object, from_files = F, chunk_max_Mb = 256)
}
\arguments{
\item{object}{R object for chunk size estimation.}
}
\value{
Estimated number of chunks with max 256MB per job.
}
\description{
Estimate number of chunks based on object size
}
\examples{
estimate_nchunks(data)
}
