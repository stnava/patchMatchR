% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchMatch.R
\name{RANSACAlt}
\alias{RANSACAlt}
\title{Alternative random sample consensus (\)}
\usage{
RANSACAlt(
  fixedPoints,
  movingPoints,
  transformType = "Affine",
  nToTrim = 2,
  minProportionPoints = 0.5,
  nCVGroups = 0,
  lambda = 1e-06,
  domainImage = NULL,
  numberOfFittingLevels = 4,
  meshSize = 1,
  dataWeights = NULL,
  verbose = FALSE
)
}
\arguments{
\item{fixedPoints}{fixed points matrix}

\item{movingPoints}{moving points matrix}

\item{transformType}{Affine, Rigid and Similarity currently supported}

\item{nToTrim}{the number of points to throw away at each iteration; if this
is a two-vector then we will sample values between these; can be used to
accelerate the search with a higher number earlier.}

\item{minProportionPoints}{the minimum proportion of points to return}

\item{nCVGroups}{number of cross-validation groups to determine error}

\item{lambda}{ridge penalty in zero to one}

\item{domainImage}{image defining the domain for deformation maps.}

\item{numberOfFittingLevels}{integer specifying the number of fitting levels.}

\item{meshSize}{vector defining the mesh size at the initial fitting level.}

\item{dataWeights}{vector defining the individual weighting of the corresponding
scattered data value.  Default = NULL meaning all values are weighted the same.}

\item{verbose}{boolean}
}
\value{
output list contains best fitted model, inliers, outliers
}
\description{
Alternative random sample consensus (\)
}
