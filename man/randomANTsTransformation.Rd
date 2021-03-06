% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchMatch.R
\name{randomANTsTransformation}
\alias{randomANTsTransformation}
\title{random transformation matrix}
\usage{
randomANTsTransformation(
  loctx,
  transformType = c("Rigid", "Affine", "ScaleShear"),
  sdAffine,
  idparams,
  fixParams,
  seeder,
  idim
)
}
\arguments{
\item{loctx}{initial affine transformation to modify}

\item{transformType}{one of Rigid, Affine and ScaleShear}

\item{sdAffine}{standard deviation parameter e.g. 0.15}

\item{idparams}{identity parameters}

\item{fixParams}{fixed parameters for ANTs or ITK transformation}

\item{seeder}{random seed}

\item{idim}{spatial dimension of transformation}
}
\value{
the input transformation matrix with parameters generated randomly
}
\description{
generates a random transformation matrix in ants style.  takes an initial
transformation, its parameters and a seed that allows repeatability.
}
\examples{

}
\author{
Avants BB
}
