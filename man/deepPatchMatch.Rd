% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchMatch.R
\name{deepPatchMatch}
\alias{deepPatchMatch}
\title{patch match two images with deep features}
\usage{
deepPatchMatch(movingImage, fixedImage, movingImageMask, fixedImageMask,
  movingPatchSize = 32, fixedPatchSize = 32, knn = 1)
}
\arguments{
\item{movingImage}{input image from which we extract patches that are
transformed to the space of the fixed image}

\item{fixedImage}{input image that provides the fixed reference domain.}

\item{movingImageMask}{defines the object of interest in the movingImage}

\item{fixedImageMask}{defines the object of interest in the fixedImage}

\item{movingPatchSize}{integer greater than or equal to 32.}

\item{fixedPatchSize}{integer greater than or equal to 32.}

\item{knn}{k-nearest neighbors ( should be 1, for now )}
}
\value{
correspondence data
}
\description{
High-level function for deep patch matching that makes many assumptions and
therefore minimizes the number of parameters the user needs to choose.
}
\examples{

library( keras )
library( ANTsR )
nP1 = 5
nP2 = 20
psz = 32
img <- ri( 1 ) \%>\% iMath( "Normalize" )
img2 <- ri( 2 ) \%>\% iMath( "Normalize" )
mask = randomMask( getMask( img ), nP1 )
mask2 = randomMask( getMask( img2 ), nP2 )
match = deepPatchMatch( img2, img, mask, mask2 )
\dontrun{
library( ANTsR )
img <- ri( 1 ) \%>\% iMath( "Normalize" ) \%>\% resampleImage( c( 2, 2 ) )
nP1 = 10
nP2 = 40
psz = 32
mask = randomMask( getMask( img ), nP1 )
features = deepFeatures( img, mask, patchSize = psz )
img2 <- ri( 5 ) \%>\% iMath( "Normalize" ) \%>\% resampleImage( c( 2, 2 ) )
txStretch = createAntsrTransform( "AffineTransform", dim=2 )
params = getAntsrTransformParameters( txStretch )
params[1] = 0.8
setAntsrTransformParameters(txStretch, params)
cos45 = cos(pi*45/180)
sin45 = sin(pi*45/180)
txRotate <- createAntsrTransform( precision="float", type="AffineTransform", dim=2 )
setAntsrTransformParameters(txRotate, c(cos45,-sin45,sin45,cos45,0,0) )
setAntsrTransformFixedParameters(txRotate, c(128,128))
rotateFirst = composeAntsrTransforms(list(txStretch, txRotate))
# img2 = applyAntsrTransform(rotateFirst, img2, img2)
mask2 = randomMask( getMask( img2 ), nP2 )
match = deepPatchMatch( img2, img, mask2, mask, 64, 64 )

for ( k in 1:nrow( match$matches ) ) {
  if ( ! is.na( match$matches[k,1] ) ) {
    layout( matrix(1:2,nrow=1) )
    plot( as.antsImage( match$ffeatures$patches[k,,] ) )
    plot( as.antsImage( match$mfeatures$patches[match$matches[k,1],,] ) )
    print( k )
    print( match$ffeatures$patchCoords[k,] )
    print( match$mfeatures$patchCoords[match$matches[k,1],] )
    Sys.sleep(1)
    }
}
}

}
\author{
Avants BB
}