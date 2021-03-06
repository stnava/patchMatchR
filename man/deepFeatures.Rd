% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchMatch.R
\name{deepFeatures}
\alias{deepFeatures}
\title{extract deep features from 2D or 3D image}
\usage{
deepFeatures(
  x,
  mask,
  patchSize = 64,
  featureSubset,
  block_name = "block2_conv2",
  vggmodel,
  subtractor = 127.5,
  patchVarEx = 0.95,
  meanCenter = FALSE
)
}
\arguments{
\item{x}{input input image}

\item{mask}{defines the object of interest in the fixedImage}

\item{patchSize}{vector or scalar defining patch dimensions}

\item{featureSubset}{a vector that selects a subset of features}

\item{block_name}{name of vgg feature block, either block2_conv2 or integer.
use the former for smaller patch sizes.  Or try ripmmarc.}

\item{vggmodel}{prebuilt feature model}

\item{subtractor}{value to subtract when scaling image intensity; should be
chosen to match training paradigm eg 127.5 for vgg and 0.5 for resnet like.}

\item{patchVarEx}{patch variance explained for ripmmarc}

\item{meanCenter}{boolean mean center the patch for ripmmarc}
}
\value{
feature array, patches and patch coordinates
}
\description{
High-level function for extracting features based on a pretrained network.
}
\examples{

library(ANTsR)
img <- ri( 1 ) \%>\% iMath( "Normalize" )
mask = randomMask( getMask( img ), 20 )
features = deepFeatures( img, mask, patchSize = 32 )

}
\author{
Avants BB
}
