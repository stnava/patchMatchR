% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchMatch.R
\name{MSE}
\alias{MSE}
\title{Mean square error of a single image or between two images.}
\usage{
MSE(x, y = NULL)
}
\arguments{
\item{x}{input image.}

\item{y}{input image.}
}
\value{
the mean squared error
}
\description{
Mean square error of a single image or between two images.
}
\examples{

library( ANTsR )

r16 <- antsImageRead( getANTsRData( 'r16' ) )
r85 <- antsImageRead( getANTsRData( 'r85' ) )
mseValue <- MSE( r16, r85 )

}
\author{
Avants BB (from redr)
}
