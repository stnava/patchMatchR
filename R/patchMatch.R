#' patch match two images
#'
#' High-level function for patch matching that makes many assumptions and
#' therefore minimizes the number of parameters the user needs to choose.
#' This prioritizes usability at the cost of optimality.
#'
#' @param movingImage input image from which we extract patches that are
#' transformed to the space of the fixed image
#' @param fixedImage input image that provides the fixed reference domain.
#' @param fixedImageMask defines the object of interest in the fixedImage
#' @param finalTransform defaults to "Rigid" but can be any antsRegistration
#' @param fixedPatchRadius integer greater than zero.
#' @param visualize boolean, will plot to screen
#' @param verbose boolean, will print to screen
#' @return data frame of corresponding points
#' @author Avants BB
#' @examples
#'
#' library(ANTsR)
#' img <- ri( 1 ) %>% iMath( "Normalize" )
#' img2 <- ri( 2 ) %>% iMath( "Normalize" )
#' mask = randomMask( getMask( img ), 2 )
#' match = patchMatch( img2, img, mask, fixedPatchRadius = 3 )
#'
#' @export patchMatch
#' @importFrom graphics plot rasterImage rect plot.new text layout
#' @importFrom ANTsR getCenterOfMass
#' @importFrom ANTsRCore antsRegistration antsApplyTransforms applyAntsrTransformToImage antsApplyTransformsToPoints antsGetSpacing applyAntsrTransformToImage createAntsrTransform  cropIndices getNeighborhoodInMask iMath antsTransformPhysicalPointToIndex antsImageClone antsImageMutualInformation
#' @importFrom ANTsRCore antsGetDirection antsGetOrigin resampleImage labelStats
patchMatch <- function(
  movingImage,
  fixedImage,
  fixedImageMask,
  finalTransform = 'Rigid',
  fixedPatchRadius = 31,
  visualize = FALSE,
  verbose = FALSE ) {

  if ( visualize ) layout( matrix( 1:3, nrow=1 ))
  # we use this to extract spatial information, not an actual matrix
  intmat0 = getNeighborhoodInMask( fixedImage, fixedImageMask,
    rep( 0, fixedImage@dimension ),
    boundary.condition = "image", spatial.info = TRUE,
    physical.coordinates = TRUE, get.gradient = FALSE)

  intmat0ind = getNeighborhoodInMask( fixedImage, fixedImageMask,
    rep( 0, fixedImage@dimension ),
    boundary.condition = "image", spatial.info = TRUE,
    physical.coordinates = FALSE, get.gradient = FALSE)

  i1reg = antsRegistration(  fixedImage, movingImage, 'SyN', verbose=F,
        gradStep = 0.1, regIterations = c(40, 20, 0 ), synSampling=2,
        totalSigma=3, flowSigma=6, synMetric='CC' )

  i1r = i1reg$warpedmovout
  mapPts = antsApplyTransformsToPoints( fixedImage@dimension,
      intmat0$indices, transformlist = i1reg$fwdtransforms  )

  txtype = "Euler2DTransform"
  if ( fixedImage@dimension == 3 )
    txtype = "Euler3DTransform"

  off = rep( fixedPatchRadius, fixedImage@dimension )
  scl = antsGetSpacing( movingImage )/antsGetSpacing( fixedImage )
  searchOff = max( round( scl ) )
  off2 = round( off / scl ) - 1
  if ( verbose ) {
    print( paste( "Search Offset:", searchOff ) )
    print( paste( "Scale Difference:", paste0(scl, collapse='x' ) ) )
    }
  nna = rep( NA, nrow( intmat0ind$indices ) )
  outdf = data.frame(
    indicesFixed = intmat0ind$indices,
    spatialFixed = intmat0$indices,
    indicesMoving = intmat0ind$indices,
    spatialMoving = intmat0$indices,
    PSNR=nna, SSIM=nna, MI=nna )
  spatminds = grep( "spatialMoving", colnames( outdf ) )
  inidminds = grep( "indicesMoving", colnames( outdf ) )
  didviz = 0
  for ( k in (1:nrow( intmat0$indices )) ) {
    locind = intmat0ind$indices[k,]
    indlo = locind - off
    indhi = locind + off + 1
    idim = dim( fixedImage )
    if ( all( indlo > 0 ) & all( indhi <= idim )  ) {
      i0patch = cropIndices( fixedImage, indlo, indhi )
      if ( verbose & didviz == 0 ) {
        print( paste( "fixPatch",paste0( dim( i0patch ), collapse='x' ) ) )
        didviz = 1
        }
      mapInd = antsTransformPhysicalPointToIndex(
        movingImage,as.numeric(mapPts[k,]))
      indlo = round( mapInd ) - off2*searchOff
      indhi = round( mapInd ) + off2*searchOff + 1
      if ( all( indlo > 0 ) & all( indhi <= dim(movingImage) )  ) {
        i1patch = cropIndices( movingImage, indlo, indhi )
      if ( verbose & didviz == 1 ) {
        print( paste( "movPatch",paste0( dim( i1patch ), collapse='x' ) ) )
        didviz = 2
        }
      centerOfMassTemplate <- getCenterOfMass( i0patch*0+1 )
      centerOfMassImage <- getCenterOfMass( i1patch * 0 + 1 )
      xfrm <- createAntsrTransform( type = txtype,
        center = centerOfMassTemplate,
        translation = centerOfMassImage - centerOfMassTemplate )
      i1rpatch = applyAntsrTransformToImage( xfrm, i1patch, i0patch )
      fix = iMath( i0patch, "Normalize" )
      mov = iMath( i1patch, "Normalize" )
      trans=antsRegistration(
        fix,
        mov, 'Translation',initialTransform=xfrm,
        affMetric = 'Mattes',
        affSampling = 20, affIterations = c( 50, 20, 20  ) )
      rigid=antsRegistration(
        fix,
        mov, 'Rigid',
        initialTransform=trans$fwdtransforms,
        affMetric = 'Mattes',
        affSampling = 20, affIterations = c(  50, 20, 20 ) )
      i1rpatchBreg=antsRegistration(
          fix,
          mov, finalTransform,
          initialTransform=rigid$fwdtransforms,
          affMetric = 'Mattes',
          affSampling = 32, affIterations = c( 11, 5 ) )
      i1rpatchB = antsApplyTransforms( fix, mov, i1rpatchBreg$fwdtransforms,
        interpolator = "bSpline" )
      mypt = antsApplyTransformsToPoints( fixedImage@dimension,
            matrix(centerOfMassTemplate, ncol=fixedImage@dimension),
            transformlist = i1rpatchBreg$fwdtransforms  )
      outdf[k, spatminds ] = mypt
      mapInd = antsTransformPhysicalPointToIndex( movingImage, as.numeric(mypt))
      outdf[k, inidminds ] = mapInd
      indlo = round( mapInd ) - off2
      indhi = round( mapInd ) + off2 + 1
      idim = dim( movingImage )
      if ( all( indlo > 0 ) & all( indhi <= idim )  ) {
        i1patchB = cropIndices( movingImage, indlo, indhi )
        if ( verbose & didviz == 2 ) {
          print( paste( "movPatchFinal",paste0( dim( i1patchB ), collapse='x' ) ) )
          didviz = 3
          }
        centerOfMassImage <- getCenterOfMass( i1patchB * 0 + 1 )
        xfrm <- createAntsrTransform( type = txtype,
          center = centerOfMassTemplate,
          translation = centerOfMassImage - centerOfMassTemplate )
          #  i1rpatchB = applyAntsrTransformToImage( xfrm, i1patchB, i0patch )
          mypsn = PSNR( i0patch, i1rpatch  )
          myssm = SSIM( i0patch, i1rpatch  )
          mymi = antsImageMutualInformation( i0patch, i1rpatch )
          mypsnB = PSNR( i0patch, i1rpatchB  )
          myssmB = SSIM( i0patch, i1rpatchB  )
          mymiB = antsImageMutualInformation( i0patch, i1rpatchB )
          outdf[ k, "PSNR" ] = mypsnB
          outdf[ k, "SSIM" ] = myssmB
          outdf[ k, "MI" ] = mymiB
          if ( visualize ) {
            plot(i0patch*10000,doCropping=F)
            plot(i1rpatchB*10000,doCropping=F)
            plot(i1patch*10000,doCropping=F)
            print( paste( k, mymi, mymiB, myssm , myssmB ) )
            Sys.sleep( 3 )
            }
          }
        }}
      }
  return( outdf )
}

#' Mean square error of a single image or between two images.
#'
#' @param x input image.
#' @param y input image.
#'
#' @return the mean squared error
#' @author Avants BB (from redr)
#' @examples
#'
#' library( ANTsR )
#'
#' r16 <- antsImageRead( getANTsRData( 'r16' ) )
#' r85 <- antsImageRead( getANTsRData( 'r85' ) )
#' mseValue <- MSE( r16, r85 )
#'
#' @export
MSE <- function( x, y = NULL )
{
  if( is.null( y ) )
    {
    return( mean( x^2 ) )
    } else {
    return( mean ( ( x - y )^2 ) )
    }
}

#' Peak signal-to-noise ratio between two images.
#'
#' @param x input image.
#' @param y input image.
#'
#' @return the peak signal-to-noise ratio
#' @author Avants BB
#' @examples
#'
#' library( ANTsR )
#'
#' r16 <- antsImageRead( getANTsRData( 'r16' ) )
#' r85 <- antsImageRead( getANTsRData( 'r85' ) )
#' psnrValue <- PSNR( r16, r85 )
#'
#' @export
PSNR <- function( x, y )
{
  return( 20 * log10( max( x ) ) - 10 * log10( MSE( x, y ) ) )
}

#' Structural similarity index (SSI) between two images.
#'
#' Implementation of the SSI quantity for two images proposed in
#'
#' Z. Wang, A.C. Bovik, H.R. Sheikh, E.P. Simoncelli. "Image quality
#' assessment: from error visibility to structural similarity". IEEE TIP.
#' 13 (4): 600–612.
#'
#' @param x input image.
#' @param y input image.
#' @param K vector of length 2 which contain SSI parameters meant to stabilize
#' the formula in case of weak denominators.
#'
#' @return the structural similarity index
#' @author Avants BB
#' @examples
#'
#' library( ANTsR )
#'
#' r16 <- antsImageRead( getANTsRData( 'r16' ) )
#' r85 <- antsImageRead( getANTsRData( 'r85' ) )
#' ssimValue <- SSIM( r16, r85 )
#'
#' @export
SSIM <- function( x, y, K = c( 0.01, 0.03 ) )
{
  globalMax <- max( max( x ), max( y ) )
  globalMin <- abs( min( min( x ), min( y ) ) )
  L <- globalMax - globalMin

  C1 <- ( K[1] * L )^2
  C2 <- ( K[2] * L )^2
  C3 <- C2 / 2

  mu_x <- mean( x )
  mu_y <- mean( y )

  mu_x_sq <- mu_x * mu_x
  mu_y_sq <- mu_y * mu_y
  mu_xy <- mu_x * mu_y

  sigma_x_sq <- mean( x * x ) - mu_x_sq
  sigma_y_sq <- mean( y * y ) - mu_y_sq
  sigma_xy <- mean( x * y ) - mu_xy

  numerator <- ( 2 * mu_xy + C1 ) * ( 2 * sigma_xy + C2 )
  denominator <- ( mu_x_sq + mu_y_sq + C1 ) * ( sigma_x_sq + sigma_y_sq + C2 )

  SSI <- numerator / denominator

  return( SSI )
}


#' extract matched patches between two images
#'
#' provides the matched patches given output of \code{patchMatch}
#'
#' @param movingImage input image from which we extract patches that are
#' transformed to the space of the fixed image
#' @param fixedImage input image that provides the fixed reference domain.
#' @param patchMatchOutput the data frame output from \code{patchMatch}.
#' @param fixedPatchRadius integer greater than zero.
#' @param verbose boolean, will print to screen.
#' @return lists of corresponding patches
#' @author Avants BB
#' @examples
#'
#' library(ANTsR)
#' img <- ri( 1 ) %>% iMath( "Normalize" )
#' img2 <- ri( 2 ) %>% iMath( "Normalize" )
#' mask = randomMask( getMask( img ), 2 )
#' match = patchMatch( img2, img, mask, fixedPatchRadius = 3 )
#' myMatches = matchedPatches( img2, img, match, fixedPatchRadius = 3 )
#'
#' @export
matchedPatches <- function(
  movingImage,
  fixedImage,
  patchMatchOutput,
  fixedPatchRadius = 31,
  verbose = FALSE )
{
  off = rep( fixedPatchRadius, fixedImage@dimension )
  scl = antsGetSpacing( movingImage )/antsGetSpacing( fixedImage )
  searchOff = max( round( scl ) )
  off2 = round( off / scl ) - 1
  cnms = colnames( patchMatchOutput )
  spatfinds = grep( "spatialFixed", cnms )
  inidfinds = grep( "indicesFixed", cnms )
  spatminds = grep( "spatialMoving", cnms )
  inidminds = grep( "indicesMoving", cnms )

  if ( verbose ) {
    print( paste( "Search Offset:", searchOff ) )
    print( paste( "Scale Difference:", paste0(scl, collapse='x' ) ) )
    }

  fixPatchList = list()
  movPatchList = list()
  didviz = 0
  for ( k in (1:nrow( patchMatchOutput ) ) ) {
    if ( ! is.na( patchMatchOutput$PSNR[k]  ) )
      {
      locind = as.numeric( patchMatchOutput[k,inidfinds] )
      indlo = locind - off
      indhi = locind + off + 1
      idim = dim( fixedImage )
      if ( all( indlo > 0 ) & all( indhi <= idim )  ) {
        i0patch = cropIndices( fixedImage, indlo, indhi )
        if ( verbose & didviz == 0 ) {
          print( paste( "fixPatch",paste0( dim( i0patch ), collapse='x' ) ) )
          didviz = 1
          }
        fixPatchList[[k]] = i0patch
        } else fixPatchList[[k]] = NA
      mapInd = as.numeric( patchMatchOutput[k,inidminds] )
      indlo = round( mapInd ) - off2
      indhi = round( mapInd ) + off2 + 1
      idim = dim( movingImage )
      if ( all( indlo > 0 ) & all( indhi <= idim )  ) {
        i1patch = cropIndices( movingImage, indlo, indhi )
        if ( verbose & didviz == 1 ) {
          print( paste( "movPatch",paste0( dim( i1patch ), collapse='x' ) ) )
          didviz = 2
          }
        movPatchList[[k]] = i1patch
      } else movPatchList[[k]] = NA
    }
  }
  return( list( fixPatchList = fixPatchList, movPatchList = movPatchList ) )
}
