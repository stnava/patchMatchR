#' make a coordinate image
#'
#' This function will return n-dimensionality coordinate images that encode
#' the spatial ITK coordinates of the image domain.
#'
#' @param mask defining wherein we create coordinates
#' @param physicalCoordinates boolean select index or phyiscal coordinates
#' @return list of coordinate images
#' @author Avants BB
#' @examples
#'
#' library(ANTsR)
#' img <- ri( 1 )
#' coords = coordinateImages( img * 0 + 1, TRUE )
#'
#' @export
coordinateImages <- function( mask, physicalCoordinates = TRUE ) {
  temp = getNeighborhoodInMask( mask, mask,
    rep( 0, mask@dimension ),
    boundary.condition = "image", spatial.info = TRUE,
    physical.coordinates = physicalCoordinates, get.gradient = FALSE)
  ilist = list()
  for ( i in 1:mask@dimension )
    ilist[[i]] = makeImage( mask, temp$indices[,i] )
  ilist
}



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
#' @param finalTransform defaults to "Rigid" but can be any \code{antsRegistration}.
#' @param fixedPatchRadius integer greater than zero.
#' @param initialMap this is the output of a call to \code{antsRegistration}.
#' if it is not present, a quick SyN map will be computed.
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
#' @importFrom stats predict rnorm sd
#' @importFrom magrittr %>%
#' @importFrom graphics plot rasterImage rect plot.new text layout
#' @importFrom ANTsR getCenterOfMass sparseDistanceMatrixXY makePointsImage
#' @importFrom ANTsRCore antsRegistration antsApplyTransforms applyAntsrTransformToImage antsApplyTransformsToPoints antsGetSpacing applyAntsrTransformToImage createAntsrTransform  cropIndices getNeighborhoodInMask iMath antsTransformPhysicalPointToIndex readAntsrTransform antsImageClone antsImageMutualInformation
#' @importFrom ANTsRCore antsGetDirection antsGetOrigin resampleImage labelStats antsTransformIndexToPhysicalPoint applyAntsrTransformToPoint antsSetSpacing getAntsrTransformFixedParameters invertAntsrTransform randomMask labelClusters getAntsrTransformParameters getCentroids getMask readAntsrTransform ri
patchMatch <- function(
  movingImage,
  fixedImage,
  fixedImageMask,
  finalTransform = 'Rigid',
  fixedPatchRadius = 31,
  initialMap,
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
  if ( missing( initialMap ) ) {
    if ( verbose ) print("Begin SyNning")
    initialMap = antsRegistration(
      fixedImage, movingImage, 'SyN', verbose=FALSE,
        gradStep = 0.1, regIterations = c(40, 20, 0 ), synSampling=2,
        totalSigma=3, flowSigma=6, synMetric='CC' )
    if ( verbose ) print( "Confess!" )
    }

  i1r = initialMap$warpedmovout
  if ( verbose ) print("map points")
  mapPts = antsApplyTransformsToPoints( fixedImage@dimension,
      intmat0$indices, transformlist = initialMap$fwdtransforms  )

  txtype = "Euler2DTransform"
  if ( fixedImage@dimension == 3 )
    txtype = "AffineTransform"

  if ( verbose ) print( paste( "txtype =", txtype ) )

  off = rep( fixedPatchRadius, fixedImage@dimension )
  scl = antsGetSpacing( movingImage ) / antsGetSpacing( fixedImage )
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
    MSE = nna, PSNR=nna, SSIM=nna, MI=nna )
  spatminds = grep( "spatialMoving", colnames( outdf ) )
  inidminds = grep( "indicesMoving", colnames( outdf ) )
  didviz = 0
  if ( verbose ) print("Begin point-wise optimization")
  for ( k in (1:nrow( intmat0$indices )) ) {
    locind = as.numeric( intmat0ind$indices[k,] )
    indlo = locind - off
    indhi = locind + off + 1
    idim = dim( fixedImage )
    if ( all( indlo > 0 ) & all( indhi <= idim )  ) {
      i0patch = cropIndices( fixedImage, indlo, indhi )
      if ( verbose & didviz == 0 ) {
        print( paste( "fixPatch",paste0( dim( i0patch ), collapse='x' ) ) )
        didviz = 1
        }
      mapInd = as.numeric( antsTransformPhysicalPointToIndex(
        movingImage,as.numeric(mapPts[k,])) )
      outdf[k, spatminds ] = as.numeric(mapPts[k,])
      outdf[k, inidminds ] = mapInd
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
        mov, 'Translation', initialTransform=xfrm,
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
      mypt2 = antsApplyTransformsToPoints( fixedImage@dimension,
            matrix(centerOfMassTemplate, ncol=fixedImage@dimension),
            transformlist = i1rpatchBreg$fwdtransforms  )
      mapInd2 = antsTransformPhysicalPointToIndex( movingImage,
        as.numeric( mypt2 ) )
      indlo = round( mapInd ) - off2
      indhi = round( mapInd ) + off2 + 1
      idim = dim( movingImage )
      if ( all( indlo > 0 ) & all( indhi <= idim )  ) {
        i1patchB = cropIndices( movingImage, indlo, indhi )
        if ( verbose & didviz == 2 ) {
          print( paste( "movPatchFinal",
            paste0( dim( i1patchB ), collapse='x' ) ) )
          didviz = 3
          }
        mymi = antsImageMutualInformation( i0patch, i1rpatch, sampling.percentage=0.8, nBins=16 )
        mymiB = antsImageMutualInformation( i0patch, i1rpatchB, sampling.percentage=0.8, nBins=16 )
        if ( is.na( mymiB ) ) mymiB = Inf
	if ( is.na( mymi ) ) mymi = Inf
	if ( mymiB < mymi ) {
          outdf[ k, spatminds ] = mypt2
          outdf[ k, inidminds ] = mapInd2
          outdf[ k, "MI" ] = mymiB
          outdf[ k, "MSE" ] = MSE( i0patch, i1rpatchB  )
          outdf[ k, "PSNR" ] = PSNR( i0patch, i1rpatchB  )
          outdf[ k, "SSIM" ] = SSIM( i0patch, i1rpatchB  )
          } else {
          outdf[ k, "MI" ] = mymi
          outdf[ k, "MSE" ] = MSE( i0patch, i1rpatch  )
          outdf[ k, "PSNR" ] = PSNR( i0patch, i1rpatch  )
          outdf[ k, "SSIM" ] = SSIM( i0patch, i1rpatch  )
          }
        if ( visualize ) {
            plot(i0patch*10000,doCropping=F)
            plot(i1rpatchB*10000,doCropping=F)
            plot(i1patch*10000,doCropping=F)
            print( paste( k, mymi, mymiB,
              outdf[ k, "SSIM" ], outdf[ k, "PSNR" ] ) )
            Sys.sleep( verbose )
            }
          }
        }}
      if ( verbose & ( ( k %% 100 ) == 0 ) ) {
        cat( paste0( k,'..' ) )
        }
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
  x = x / max( x )
  if( is.null( y ) )
    {
    return( mean( x^2 ) )
    } else {
    y = y / max( y )
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
  x = x / max( x )
  y = y / max( y )
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
  x = x / max( x )
  y = y / max( y )
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

#' convert coordinates of deep feature patches to image-based physical coordinates
#'
#' High-level function for handling the converion between voxel and corner-based
#' indices of patches to center and physical-based landmark coordinates.
#'
#' @param patchCoords input patch coordinates (usually from \code{deepFeatures})
#' @param img reference image defining spatial domain
#' @param patchSize vector or scalar defining patch dimensions
#' @return spatial patch coordinates
#' @author Avants BB
#' @examples
#'
#' library(ANTsR)
#' img <- ri( 1 ) %>% iMath( "Normalize" )
#' mask = randomMask( getMask( img ), 20 )
#' features = deepFeatures( img, mask, patchSize = 32 )
#' coords = convertPatchCoordsToSpatialPoints( features$patchCoords, img )
#'
#' @export convertPatchCoordsToSpatialPoints
convertPatchCoordsToSpatialPoints<-function( patchCoords, img, patchSize = 32 ) {
  idim = img@dimension
  patchSizeDivBy2 = patchSize / 2
  locpts = matrix( patchCoords+patchSizeDivBy2, ncol = idim )
  locpts = antsTransformIndexToPhysicalPoint( img, locpts )
  #  ptrad = ptscl * sqrt( sum( antsGetSpacing( img )^2 ) )
  return( locpts )
}

#' patch match two images with deep features
#'
#' High-level function for deep patch matching that makes many assumptions and
#' therefore minimizes the number of parameters the user needs to choose.
#'
#' @param movingImage input image from which we extract patches that are
#' transformed to the space of the fixed image
#' @param fixedImage input image that provides the fixed reference domain.
#' @param movingImageMask defines the object of interest in the movingImage
#' @param fixedImageMask defines the object of interest in the fixedImage
#' @param movingPatchSize integer greater than or equal to 32.
#' @param fixedPatchSize integer greater than or equal to 32.
#' @param knn k-nearest neighbors ( should be >= 1  )
#' @param knnSpatial k-nearest neighbors for spatial localization (optional).
#' this will constrain the search to more proximal locations.  will perform
#' better if the images are in the same physical space. currently, the units
#' for the spatial distance is in voxels.  may add physical space option later.
#' FIXME - allow a transformation to be passed to this step s.t. moving points
#' can be transformed to fixed space before distance assessment.
#' @param featureSubset a vector that selects a subset of features
#' @param block_name name of vgg feature block, either block2_conv2 or integer.
#' use the former for smaller patch sizes.
#' @param switchMatchDirection boolean
#' @param kPackage name of package to use for knn
#' @param vggmodel prebuilt feature model
#' @param subtractor value to subtract when scaling image intensity; should be
#' chosen to match training paradigm eg 127.5 for vgg and 0.5 for resnet like.
#' @param verbose boolean
#' @return correspondence data
#' @author Avants BB
#' @examples
#'
#' library( keras )
#' library( ANTsR )
#' nP1 = 5
#' nP2 = 20
#' psz = 32
#' img <- ri( 1 ) %>% iMath( "Normalize" )
#' img2 <- ri( 2 ) %>% iMath( "Normalize" )
#' mask = randomMask( getMask( img ), nP1 )
#' mask2 = randomMask( getMask( img2 ), nP2 )
#' match = deepPatchMatch( img2, img, mask, mask2 )
#' \dontrun{
#' library( ANTsR )
#' img <- ri( 1 ) %>% iMath( "Normalize" ) %>% resampleImage( c( 2, 2 ) )
#' nP1 = 10
#' nP2 = 40
#' psz = 32
#' mask = randomMask( getMask( img ), nP1 )
#' features = deepFeatures( img, mask, patchSize = psz )
#' img2 <- ri( 5 ) %>% iMath( "Normalize" ) %>% resampleImage( c( 2, 2 ) )
#' txStretch = createAntsrTransform( "AffineTransform", dim=2 )
#' params = getAntsrTransformParameters( txStretch )
#' params[1] = 0.8
#' setAntsrTransformParameters(txStretch, params)
#' cos45 = cos(pi*45/180)
#' sin45 = sin(pi*45/180)
#' txRotate <- createAntsrTransform( precision="float", type="AffineTransform", dim=2 )
#' setAntsrTransformParameters(txRotate, c(cos45,-sin45,sin45,cos45,0,0) )
#' setAntsrTransformFixedParameters(txRotate, c(128,128))
#' rotateFirst = composeAntsrTransforms(list(txStretch, txRotate))
#' # img2 = applyAntsrTransform(rotateFirst, img2, img2)
#' mask2 = randomMask( getMask( img2 ), nP2 )
#' match = deepPatchMatch( img2, img, mask2, mask, 64, 64 )
#'
#' for ( k in 1:nrow( match$matches ) ) {
#'   if ( ! is.na( match$matches[k,1] ) ) {
#'     layout( matrix(1:2,nrow=1) )
#'     plot( as.antsImage( match$ffeatures$patches[k,,] ) )
#'     plot( as.antsImage( match$mfeatures$patches[match$matches[k,1],,] ) )
#'     print( k )
#'     print( match$ffeatures$patchCoords[k,] )
#'     print( match$mfeatures$patchCoords[match$matches[k,1],] )
#'     Sys.sleep(1)
#'     }
#' }
#' }
#'
#' @export deepPatchMatch
#' @importFrom ANTsRNet extractImagePatches extractImagePatchCoordinates createVggModel3D createVggModel2D
#' @importFrom ANTsRNet createFullyConvolutionalVggModel2D createFullyConvolutionalVggModel3D
#' @importFrom qlcMatrix colMin
#' @importFrom abind abind
#' @importFrom keras application_vgg19 keras_model get_layer get_weights set_weights
deepPatchMatch <- function(
  movingImage,
  fixedImage,
  movingImageMask,
  fixedImageMask,
  movingPatchSize = 32,
  fixedPatchSize = 32,
  knn = 1,
  knnSpatial = 0,
  featureSubset,
  block_name = 'block2_conv2',
  switchMatchDirection = FALSE,
  kPackage = 'FNN',
  vggmodel,
  subtractor = 127.5,
  verbose = FALSE )
{
  if ( verbose ) print( Sys.time() )
  if ( ! missing( vggmodel ) ) {
    if ( verbose ) print("DF1")
    ffeatures = deepFeatures( fixedImage, fixedImageMask,
      patchSize = fixedPatchSize, block_name = block_name, vggmodel=vggmodel,
      subtractor = subtractor  )
    if ( verbose ) print("DF2")
    if ( verbose ) print( Sys.time() )
    mfeatures = deepFeatures( movingImage, movingImageMask,
      patchSize = movingPatchSize, block_name = block_name, vggmodel=vggmodel,
      subtractor = subtractor    )
  }
  if (  missing( vggmodel ) )  {
    if ( verbose ) print("DF1")
    ffeatures = deepFeatures( fixedImage, fixedImageMask,
      patchSize = fixedPatchSize, block_name = block_name,
      subtractor = subtractor   )
    if ( verbose ) print("DF2x")
    if ( verbose ) print( Sys.time() )
    mfeatures = deepFeatures( movingImage, movingImageMask,
      patchSize = movingPatchSize, block_name = block_name, vggmodel=ffeatures$featureModel,
      subtractor = subtractor   )
    if ( verbose ) print("DF2x-end")
    if ( verbose ) print( Sys.time() )
  }
  if ( ! missing( featureSubset ) ) {
    if ( verbose ) print("DF1-subset")
    ffeatures = deepFeatures( fixedImage, fixedImageMask, patchSize = fixedPatchSize,
      featureSubset = featureSubset, block_name = block_name,
      subtractor = subtractor  )
    if ( verbose ) print("DF2-subset")
    mfeatures = deepFeatures( movingImage, movingImageMask, patchSize = movingPatchSize,
      featureSubset = featureSubset, block_name = block_name, vggmodel=ffeatures$featureModel,
      subtractor = subtractor   )
  }
  if ( knnSpatial > 0 ) {
    if ( verbose ) print("spatial-distance-begin")
    if ( verbose ) print( Sys.time() )
    fdistmat <- antsTransformIndexToPhysicalPoint(fixedImage,ffeatures$patchCoords)
    mdistmat <- antsTransformIndexToPhysicalPoint(movingImage,mfeatures$patchCoords)
#    fdistmat <- ffeatures$patchCoords
#    mdistmat <- mfeatures$patchCoords
    # FIXME - add jitter to prevent zero distances
    fspc = sqrt( sum( antsGetSpacing(fixedImage )))
    mspc = sqrt( sum( antsGetSpacing(movingImage )))
    jitterF = matrix( rnorm(length(fdistmat),0,1e-4*fspc),
      ncol = fixedImage@dimension )
    jitterM = matrix( rnorm(length(mdistmat),0,1e-4*mspc),
      ncol = fixedImage@dimension )
    if ( !switchMatchDirection ) spatialDistMat = sparseDistanceMatrixXY(
      t(mdistmat+jitterM), t(fdistmat+jitterF),
      k = knnSpatial, kmetric='euclidean', kPackage=kPackage)
    if ( switchMatchDirection ) spatialDistMat = sparseDistanceMatrixXY(
      t(fdistmat+jitterF), t(mdistmat+jitterM),
      k = knnSpatial, kmetric='euclidean', kPackage=kPackage)
    if ( verbose ) print("spatial-distance-end")
    if ( verbose ) print( Sys.time() )
    # this will constrain the search
    spatialDistMat[ spatialDistMat > 0] = 1
  }
  if ( verbose ) print("sdxy-begin")
  if ( verbose ) print( Sys.time() )
  matches = matrix( nrow = nrow( ffeatures$patches  ), ncol = 1 )
  costs = matrix( nrow = nrow( ffeatures$patches  ), ncol = 1 )
  matchesKNN = matrix( nrow = nrow( ffeatures$patches  ), ncol = knn )
  costsKNN = matrix( nrow = nrow( ffeatures$patches  ), ncol = knn )
  if ( switchMatchDirection ) {
    mydist = sparseDistanceMatrixXY(
      t(ffeatures$features), t(mfeatures$features), k = knn,
      kmetric='euclidean', kPackage = kPackage )
    if ( knnSpatial > 0 ) mydist = mydist * spatialDistMat
    best1s = qlcMatrix::colMin( mydist, which = TRUE  )
    for ( k in 1:ncol(best1s$which) ) {
      ww = which( best1s$which[,k] )
      if ( length( ww ) > 0 ) {
        matches[ k, ] = ww[1]
        costs[k, ] = as.numeric( best1s$min[k] )
      }
      ww = which( mydist[,k] > 0 )
      if ( length( ww ) > 0 ) {
        availind = 1:length( ww )
        matchesKNN[ k, availind ] = ww
        costsKNN[ k, availind ] = mydist[ matchesKNN[ k, availind ], k ]
      }
    }
  } else {
    mydist = sparseDistanceMatrixXY(
      t(mfeatures$features), t(ffeatures$features), k = knn,
      kmetric='euclidean', kPackage = kPackage )
    if ( verbose ) print( Sys.time() )
    if ( knnSpatial > 0 ) mydist = mydist * spatialDistMat
    best1s = qlcMatrix::rowMin( mydist, which = TRUE  )
    for ( k in 1:nrow(best1s$which) ) {
      ww = which( best1s$which[k,] )
      if ( length( ww ) > 0 ) {
        matches[ k, ] = ww[1]
        costs[k, ] = as.numeric( best1s$min[k] )
      } # length
      ww = which( mydist[k,] > 0 )
      if ( length( ww ) > 0 ) {
        availind = 1:length( ww )
        matchesKNN[ k, availind ] = ww
        costsKNN[k, availind ] = mydist[k, matchesKNN[ k, availind ] ]
        }
    } # row
  } # else
  if ( verbose ) print("sdxy-fin")
  if ( verbose ) print( Sys.time() )
  if ( knn > 1 ) {
    # if knn > 0 => make probabilities
    # (exp( -1.0 * (vv-min(vv))^2/ (mean(vv)*9) ))
    }
  return(
    list(
      distanceMatrix = mydist,
      ffeatures = ffeatures,
      mfeatures = mfeatures,
      matches = matches,
      costs = costs,
      matchesK = matchesKNN,
      costsK = costsKNN )
   )
}




#' locally constrained deep feature patch matching
#'
#' High-level function for deep patch matching constrained to be local.
#'
#' @param movingImage input image from which we extract patches that are
#' transformed to the space of the fixed image
#' @param fixedImage input image that provides the fixed reference domain.
#' @param movingImageMask defines the object of interest in the movingImage
#' @param fixedImageMask defines the object of interest in the fixedImage
#' @param patchSize integer greater than or equal to 32.
#' @param knn k-nearest neighbors ( should be >= 1  )
#' @param localSearchRadius radius value passed to \code{makePointsImage}
#' @param nSamples number of local samples (optional), can speed things up at
#' the cost of some accuracy
#' @param block_name name of vgg feature block, either block2_conv2 or integer.
#' use the former for smaller patch sizes.
#' @param kPackage name of package to use for knn
#' @param verbose boolean
#' @return correspondence data
#' @author Avants BB
#' @examples
#'
#' library( keras )
#' library( ANTsR )
#' nP1 = 5
#' nP2 = 20
#' psz = 32
#' img <- ri( 1 ) %>% iMath( "Normalize" )
#' img2 <- ri( 2 ) %>% iMath( "Normalize" )
#' mask = randomMask( getMask( img ), nP1 )
#' mask2 = randomMask( getMask( img2 ), nP2 )
#' match = deepLocalPatchMatch( img2, img, mask, mask2 )
#'
#' @export deepLocalPatchMatch
deepLocalPatchMatch <- function(
  movingImage,
  fixedImage,
  movingImageMask,
  fixedImageMask,
  patchSize = 32,
  knn = 1,
  localSearchRadius = 5,
  nSamples,
  block_name = 'block2_conv2',
  kPackage = 'FNN',
  verbose = FALSE )
{
  ffeatures = deepFeatures( fixedImage, fixedImageMask,
      patchSize = patchSize, block_name = block_name  )
  featureModel = ffeatures$featureModel
  # create a point mask for local regions
  createPointMask <- function( ptPhys, referenceMask, radius, nsamples ) {
    if ( ! missing( nsamples ) ) {
      randomMask( makePointsImage( matrix( ptPhys, nrow = 1 ), referenceMask,
        radius=radius ),  nsamples )
      } else makePointsImage( matrix( ptPhys, nrow = 1 ), referenceMask,
        radius=radius )
  }

  if ( verbose ) print("local-match-begin")
  matches = matrix( nrow = nrow( ffeatures$patches  ), ncol = 1 )
  costs = matrix( nrow = nrow( ffeatures$patches  ), ncol = 1 )
  matchedCoords = fixedCoords = ffeatures$patchCoords
  fixedCoords[ ] = matchedCoords[ ] = NA
  for ( i in 1:nrow( ffeatures$patches  ) ) {
    cat(paste(i,'...'))
    off = 0
    fpt = ffeatures$patchCoords[i, ] + off
    ptPhys = antsTransformIndexToPhysicalPoint( fixedImageMask, fpt )
    fixedCoords[i,] = ptPhys
    if ( missing( nSamples ) )
      localMask = createPointMask( ptPhys, movingImageMask, localSearchRadius )
    else localMask = createPointMask( ptPhys, movingImageMask, localSearchRadius,
      nSamples )
    mfeatures = deepFeatures( movingImage, localMask,
      patchSize = patchSize, block_name = block_name, vggmodel = featureModel)
    mydist = sparseDistanceMatrixXY(
      t(mfeatures$features), matrix( t(ffeatures$features[i,]),ncol=1),
      k = knn, kmetric='euclidean')
    best1s = qlcMatrix::rowMin( mydist, which = TRUE  )
    ww = which( best1s$which[1,] )
    if ( length( ww ) > 0 ) {
      matches[ i, ] = ww[1]
      costs[i, ] = as.numeric( best1s$min[1] )
      mpt = mfeatures$patchCoords[ww[1],] + off
      matchedCoords[i, ] = antsTransformIndexToPhysicalPoint( fixedImageMask, mpt )
      } # length( ww )
    } # row
  if ( verbose ) print("local-match-fin")
  return(
    list(
      fixedPoints = fixedCoords,
      movingPoints = matchedCoords,
      matches = matches,
      costs = costs )
   )
}


#' extract deep features from 2D or 3D image
#'
#' High-level function for extracting features based on a pretrained network.
#'
#' @param x input input image
#' @param mask defines the object of interest in the fixedImage
#' @param patchSize vector or scalar defining patch dimensions
#' @param featureSubset a vector that selects a subset of features
#' @param block_name name of vgg feature block, either block2_conv2 or integer.
#' use the former for smaller patch sizes.
#' @param vggmodel prebuilt feature model
#' @param subtractor value to subtract when scaling image intensity; should be
#' chosen to match training paradigm eg 127.5 for vgg and 0.5 for resnet like.
#' @return feature array, patches and patch coordinates
#' @author Avants BB
#' @examples
#'
#' library(ANTsR)
#' img <- ri( 1 ) %>% iMath( "Normalize" )
#' mask = randomMask( getMask( img ), 20 )
#' features = deepFeatures( img, mask, patchSize = 32 )
#'
#' @export deepFeatures
deepFeatures <- function( x, mask, patchSize = 64,
  featureSubset, block_name = 'block2_conv2', vggmodel,
  subtractor = 127.5 ) {
  idim = x@dimension
  if ( length( patchSize ) == 1 ) patchSize = rep( patchSize, idim )
  vggp = patchSize
  if ( missing( vggmodel ) ) {
    if ( idim == 2 ) {
      if ( block_name == 'block2_conv2' ) {
        vggmodel = createFullyConvolutionalVggModel2D(
          list( NULL, NULL, 3 ),
          layers = c( 1, 2, 3 ), lowestResolution = 64,
          convolutionKernelSize = c(3, 3), poolSize = c(2, 2),
          strides = c(2, 2), dropoutRate = 0, style = 19 )
        vggmodel <- keras_model( inputs = vggmodel$input,
          outputs = get_layer(vggmodel, index = 6 )$output)
      } else {
      lays = c(1, 2, 3, 4, 4 )
      if ( block_name <= 6 ) lays = c( 1, 2, 3 )
      vggmodel = createFullyConvolutionalVggModel2D(
        list( NULL, NULL, 3 ),
        layers = lays, lowestResolution = 64,
        convolutionKernelSize = c(3, 3), poolSize = c(2, 2),
        strides = c(2, 2), dropoutRate = 0, style = 19 )
      vggmodel <- keras_model( inputs = vggmodel$input,
        outputs = get_layer(vggmodel, index = block_name )$output)
      }
      vgg19 = application_vgg19(
          include_top = FALSE, weights = "imagenet",
          input_shape = list( 32, 32, 3 ),
          classes = 1000)
      if ( block_name == 'block2_conv2' )
        vggmodelRaw <- keras_model( inputs = vgg19$input,
              outputs = get_layer(vgg19, block_name )$output)
      if ( is.numeric( block_name ) ) {
        vggmodelRaw <- keras_model( inputs = vgg19$input,
              outputs = get_layer(vgg19, index = block_name + 1 )$output)
      }
      set_weights( vggmodel, get_weights( vggmodelRaw ) )
    }
    if ( idim == 3 ) {
      vgg19 = application_vgg19(
        include_top = FALSE, weights = "imagenet",
        input_shape = c( 32, 32, 3 ),
        classes = 1000)
      if ( block_name == 'block2_conv2' )
        vggmodel2D <- keras_model( inputs = vgg19$input,
          outputs = get_layer(vgg19, block_name )$output)
      if ( is.numeric( block_name ) )
        vggmodel2D <- keras_model( inputs = vgg19$input,
          outputs = get_layer(vgg19, index = block_name + 1 )$output)
      ######################################################################################
      nchan = 1
      if ( block_name == 'block2_conv2' ) {
        vggmodel = createFullyConvolutionalVggModel3D(
          list( NULL, NULL, NULL, 1 ),
          layers = c( 1, 2, 3 ), lowestResolution = 64,
          convolutionKernelSize = c(3, 3, 3), poolSize = c(2, 2, 2),
          strides = c(2, 2, 2),  dropoutRate = 0, style = 19 )
        vggmodel <- keras_model( inputs = vggmodel$input,
          outputs = get_layer(vggmodel, index = 6 )$output)
      } else {
        lays = c(1, 2, 3, 4, 4 )
        if ( block_name <= 6 ) lays = c( 1, 2, 3 )
        vggmodel = createFullyConvolutionalVggModel3D(
          list( NULL, NULL, NULL, 1 ),
          lays, lowestResolution = 64,
          convolutionKernelSize = c(3, 3, 3), poolSize = c(2, 2, 2),
          strides = c(2, 2, 2),  dropoutRate = 0, style = 19 )
        vggmodel <- keras_model( inputs = vggmodel$input,
          outputs = get_layer(vggmodel, index = block_name )$output)
      }
      vgg3Dweights = get_weights( vggmodel )
      vgg2Dweights = get_weights( vggmodel2D )
      for ( j in 1:nchan )
        vgg3Dweights[[1]][ , , , j , ] = vgg2Dweights[[1]] / nchan
      for ( k in seq( from=2, length( vgg3Dweights ), by=2 ) )
        vgg3Dweights[[k]] = vgg2Dweights[[k]]
      for ( k in seq( from=3, length( vgg3Dweights )-1, by=2 ) )
        for ( j in 1:idim )
          vgg3Dweights[[k]][ , , j, , ] = vgg2Dweights[[k]] / idim
      set_weights( vggmodel, vgg3Dweights )
    }
  } # exists vggmodel
  x = iMath( x, "Normalize" ) * abs(subtractor*2) - subtractor
  patches0 = extractImagePatches( x, patchSize, maskImage = mask,
    maxNumberOfPatches=sum(mask), returnAsArray = T, randomSeed = 1, randomize=FALSE   )
  patchCoords = extractImagePatchCoordinates( x, patchSize, maskImage = mask,
    maxNumberOfPatches=sum(mask), physicalCoordinates = FALSE, cornerCoordinates=TRUE, randomSeed = 1, randomize=FALSE   )
  patches = patches0
  nChannels = unlist(dim(vggmodel$inputs[[1]]))
  if ( idim == 2 & nChannels == 3 ) {
    for( k in 2:3 )
      patches = abind( patches, patches0, along = idim+2)
    } else {
      patches = array( patches, dim = c( dim( patches ), 1  ) )
    }
  features = predict( vggmodel, patches )
  vecdim = prod( dim( features )[-1]  )
  if ( ! missing( featureSubset ) )
    if ( any( featureSubset > vecdim ) )
      featureSubset = featureSubset[ featureSubset <= vecdim ]
  featuresMat = as.matrix( array( features,  dim = c( nrow( features), vecdim ) ) )
  if ( ! missing( featureSubset ) )
    featuresMat = featuresMat[,featureSubset]
  return(
    list(
      features=featuresMat,
      featuresRaw=features,
      patches=patches0,
      patchCoords = patchCoords + floor( patchSize / 2 - 1 ),
      featureModel = vggmodel ) )
}

#' Fit transform to points
#'
#' This function will use either the Kabsch algorithm or a least squares fitting
#' algorithm to match the pairs of points that the user provides.  An antsr
#' transform is returned.
#'
#' @param movingPoints moving points matrix
#' @param fixedPoints fixed points matrix
#' @param transformType Affine, Rigid and Similarity currently supported
#' @param lambda ridge penalty
#' @param domainImage image defining the domain for deformation maps.
#' @param numberOfFittingLevels integer specifying the number of fitting levels.
#' @param meshSize vector defining the mesh size at the initial fitting level.
#' @param dataWeights vector defining the individual weighting of the corresponding
#' scattered data value.  Default = NULL meaning all values are weighted the same.
#' Currently, this is only relevant to BSpline fits but FIXME will be generalized.
#' @return antsTransform that maps the moving image to the fixed image space.
#' the inverse transform maps the moving points to the fixed space.  Associated
#' error is also returned.
#' @export
fitTransformToPairedPoints <-function(
  movingPoints,
  fixedPoints,
  transformType = "Affine",
  lambda = 1e-4,
  domainImage, numberOfFittingLevels=4, meshSize=1, dataWeights
 ) {
    if ( ! any( transformType %in% c( "Rigid", "Affine", "Similarity", "BSpline" ) ) )
      stop("Transform not supported")
######################
# x: fixedLandmarks
# y: movingLandmarks
# (A,t,c) : affine transform, A:3*3, t: 3*1 c: 3*1 (c is the center of all points in x)
# y-c = A*(x-c) + t;
# steps:
# 1. c = average of points of x
# 2. let y1 = y-c; x1 = x - c; x11 = [x1; 1 ... 1] # extend x11
# 3. minimize (y1-A1*x11)^2, A1 is a 3*4 matrix
# 4. A = A1(1:3, 1:3), t = A1(1:3, 4);
# step 3:
#   A11 = (y1*x11')*(x11*x11')^(-1)

  # https://github.com/ANTsX/ANTs/blob/3f3cd4b775036345a28898ca9fe5a56f04ed4973/Examples/ANTSUseLandmarkImagesToGetAffineTransform.cxx#L84-L180
  generateData = FALSE
  if ( generateData ) {
    img = ri( 1 )
    antsSetSpacing( img, rep( 0.8 , 2 ) )
    img2 = ri( 5 )
    aff = antsRegistration( img, img2, "Affine" )
    trueTx = readAntsrTransform( aff$fwdtransforms )
    trueTx = invertAntsrTransform( trueTx )
    fixedParams = getAntsrTransformFixedParameters( trueTx )
    txParams = getAntsrTransformParameters( trueTx )
    # find some fixed and moving points
    msk = getMask( img )
    rmsk = randomMask( msk, 555 ) %>% labelClusters( 1 ) %>% iMath("GD",2)
    rmskTx = antsApplyTransforms( img2, rmsk, transformlist = aff$fwdtransforms,
      whichtoinvert = c(TRUE), interpolator = 'nearestNeighbor' )
    fixedPoints = getCentroids( rmsk )[,1:2]
    movingPoints = getCentroids( rmskTx )[,1:2]
    movingPoints2 = applyAntsrTransformToPoint( trueTx, fixedPoints )
    }
  n = nrow( fixedPoints )
  idim = ncol( fixedPoints )
  if ( transformType %in% c( "Rigid", "Affine", "Similarity" ) ) {
    x = fixedPoints
    y = movingPoints

    # 1. c = average of points of x
    # 2. let y1 = y-c; x1 = x - c; x11 = [x1; 1 ... 1] # extend x11
    centerX = colMeans( x )
    centerY = colMeans( y )
    for ( i in 1:nrow( x ) ) {
      x[i,] = x[i,] - centerX
      y[i,] = y[i,] - centerY
      }
    x11 = cbind( x, rep( 1, nrow(x)))
    # 3. minimize (y1-A1*x11)^2, A1 is a 3*4 matrix
    temp = qr.solve( x11, y )
    A = t( temp[1:idim, 1:idim ] )
    trans = temp[idim+1,] + centerY - centerX
    if ( transformType %in% c("Rigid", "Similarity" ) ) {
      # http://web.stanford.edu/class/cs273/refs/umeyama.pdf
  #    Kabsch Algorithm.
      covmat = ( t( y ) %*% x )
      x_svd <- svd( covmat  + diag(idim) * lambda)
      myd = det( x_svd$u %*% t( x_svd$v ) )
      signadj = diag( idim )
      if ( myd > 0 ) A = x_svd$u %*% t( x_svd$v ) else {
        signadj = diag( c( rep( 1, idim - 1 ), -1 ) )
        A = ( x_svd$u %*% signadj ) %*% t( x_svd$v )
      }
      scaling = 1
      if ( transformType == "Similarity" ) {
        scaling =  sqrt( mean( rowSums( y^2 )/n )  ) /
                   sqrt( mean( rowSums( x^2 )/n )  )
        }
      A = A %*% ( diag( idim ) * scaling )
    }
    aff = createAntsrTransform( matrix = A, translation = trans, dimension = idim,
      center = centerX  )
    if ( generateData ) {
      aff = invertAntsrTransform( aff )
      movingPointsTx = applyAntsrTransformToPoint( aff, movingPoints )
      movingPointsTx2 = applyAntsrTransformToPoint( trueTx, movingPoints )
      print( paste( myd,
        norm( movingPointsTx - fixedPoints, "F" ),
        norm( movingPointsTx2 - fixedPoints, "F" ) ) )
      } else {
        err = norm( movingPoints - applyAntsrTransformToPoint( aff, fixedPoints ), "F" )
      }
    return( list( transform = aff, error = err/n ) )
  }
  if ( transformType == "BSpline" ) {
    if ( length( meshSize ) == domainImage@dimension ) mymeshsize = meshSize
    if ( length( meshSize ) == 1 ) mymeshsize = rep( meshSize, domainImage@dimension )
    txfpts = fixedPoints
    txmpts = movingPoints
    mydir = antsGetDirection( )
    scatteredData = movingPoints - fixedPoints
    bsplineImage <- fitBsplineObjectToScatteredData( scatteredData, fixedPoints,
           parametricDomainOrigin = antsGetOrigin(domainImage),
           parametricDomainSpacing = antsGetSpacing(domainImage),
           parametricDomainSize = dim( domainImage ),
           numberOfFittingLevels = numberOfFittingLevels, meshSize = mymeshsize )
    tx = antsrTransformFromDisplacementField( bsplineImage )
    err = norm( movingPoints - applyAntsrTransformToPoint( tx, fixedPoints ), "F" )
    return( list( transform = tx, error = err/n, displacement=bsplineImage ) )
  } else stop( "Invalid transformation type." )
}




#' Fit transform to points with tensorflow
#'
#' This function will use either the Kabsch algorithm or a least squares fitting
#' algorithm to match the pairs of points that the user provides.  A tensorflow
#' tensor is returned.
#'
#' @param movingPoints moving points matrix
#' @param fixedPoints fixed points matrix
#' @param numberOfPoints per sample in batch
#' @param dimensionality of the point set
#' @param transformType Rigid, Similarity or Affine currently supported
#' @param batch_size the batch size
#' @param preventReflection boolean
#' @return tensorflow tensor object.
#' @export
fitTransformToPairedPointsTF <-function(
  movingPoints,
  fixedPoints,
  numberOfPoints,
  dimensionality,
  transformType = c( "Rigid", "Similarity", "Affine"  ),
  batch_size = 1,
  preventReflection = TRUE ) {
  transformType = match.arg( transformType )
  if ( class( fixedPoints )[1] != "tensorflow.tensor" ) {
    xyz0 = array( fixedPoints, dim = c( batch_size, numberOfPoints, dimensionality ) )
    xyz1 = array( movingPoints, dim = c( batch_size, numberOfPoints, dimensionality ) )
  } else {
    xyz0 = tf$cast( fixedPoints, "float64" )
    xyz1 = tf$cast(movingPoints, "float64" )
  }
  cen0 = tf$reduce_mean(xyz0, 1L, keepdims=TRUE)
  cen1 = tf$reduce_mean(xyz1, 1L, keepdims=TRUE)
  xtf = xyz0 - cen0
  ytf = xyz1 - cen1
  if ( transformType == 'Affine' ) {
    myones = tf$expand_dims( tf$expand_dims(
      tf$cast( tf$ones( numberOfPoints ), 'float64' ), axis = c( 0L ) ), axis=2L )
    myonesB = tf$expand_dims( tf$expand_dims(
      tf$cast( tf$ones( numberOfPoints ), 'float64' ), axis = c( 0L ) ), axis=2L )
    if ( batch_size > 1 )
      for ( k in 1:(batch_size-1) ) myones = tf$concat( list(myones, myonesB ), axis=0L )
    xtfu = tf$concat( list( xtf, myones ), axis=2L )
    return( tf$linalg$matrix_transpose(
      tf$linalg$lstsq( xtfu, ytf )[,1:dimensionality,1:dimensionality] ) )
    }
  cov = tf$matmul( ytf, xtf,  transpose_b = F, transpose_a=T )
  mysvd = tf$linalg$svd(cov, full_matrices=T ) # returns _, u, v
  u = mysvd[[2]]
  v = mysvd[[3]]
  d = tf$linalg$det(tf$matmul(v, u, transpose_b=F))
  newu=list()
  if ( preventReflection | transformType == "Similarity" ) {
    tmp = tf$cast( numberOfPoints,'float64' )
    for ( k in 1:batch_size ) {
      newu[[k]] = u[[k-1]]
      if ( as.numeric( d[[ k - 1 ]] ) < 0 ) {
        signadj = diag( c( rep( 1, dimensionality - 1 ), -1 ) )
        newu[[k]] = tf$cast( tf$matmul( u[[k-1]], signadj ), "float64" )
        }
      if (  transformType == "Similarity" ) {
        scaling =  tf$math$sqrt( tf$reduce_mean( tf$reduce_sum( ytf[[k-1]]^2 )/tmp ) ) /
                   tf$math$sqrt( tf$reduce_mean( tf$reduce_sum( xtf[[k-1]]^2 )/tmp ) )
        scladj = diag( dimensionality ) * as.numeric( scaling )
        newu[[k]] = tf$matmul(newu[[k]], tf$cast(scladj,'float64'), transpose_b=FALSE)
        }
      }
    u = newu
  }
  return( tf$matmul(u, v, transpose_b=TRUE) ) # the rotation matrix
}


#' Random sample consensus (RANSAC)
#'
#' @param movingPoints moving points matrix
#' @param fixedPoints fixed points matrix
#' @param transformType Affine, Rigid and Similarity currently supported
#' @param minNtoFit the minimum number of data values required to fit the model.
#' @param maxIterations the maximum number of iterations allowed in the algorithm
#' @param errorThreshold a threshold value for determining when a test data point fits a model.
#' this parameter is set based on the standard deviation in the random subset model.
#' that is, a point fits the model error distribution if it is within the bracket
#' of values between mean error plus or minus sd error times errorThreshold.
#' @param goodProportion the fraction of close data values required to assert that a model
#' fits well to data.  that is, if equal to 0.5, then one would need 50 points to
#' assert that a model fit is good if the whole dataset contains 100 points.
#' @param lambda ridge penalty
#' @param verbose boolean
#'
#' @return output list contains best fitted model, inliers, outliers
#'
#' @export
RANSAC <- function(
  fixedPoints,
  movingPoints,
  transformType = "Affine",
  minNtoFit = 16,
  maxIterations = 20,
  errorThreshold = 1,
  goodProportion = 0.5,
  lambda = 1e-4,
  verbose = FALSE ) {
  # 1 Select a random subset of the original data. Call this subset the hypothetical inliers.
  # 2 A model is fitted to the set of hypothetical inliers.
  # 3 All other data are then tested against the fitted model.
  #     Those points that fit the estimated model well, according to some
  #     model-specific loss function, are considered as part of the consensus set.
  # 4 The estimated model is reasonably good if sufficiently many points have been classified as part of the consensus set.
  # 5 Afterwards, the model may be improved by reestimating it using all members of the consensus set.
  nMax = nrow( fixedPoints )
  d = round( nMax * goodProportion )
  nInliers = 0
  nIterations = 0
  bestModel = NULL
  bestErr = Inf
#  while ( nInliers < d & nIterations < maxIterations ) {
  while ( nIterations < maxIterations ) {
    nIterations = nIterations + 1
    randSubset = sample( 1:nMax, minNtoFit )         # step 1
    modelFit = fitTransformToPairedPoints(   # step 2
      movingPoints[ randSubset, ],
      fixedPoints[ randSubset, ],
      transformType = transformType, lambda = lambda )
    mapComplement = applyAntsrTransformToPoint( modelFit$transform, fixedPoints )
    err = sqrt( rowMeans( ( movingPoints - mapComplement )^2 ) )
    mn1 = mean( err[ randSubset ] )
    sd1 = sd( err[ randSubset ] )
    meanInValBracket = mn1 + sd1 * c( -1, 1 ) * errorThreshold
    inliers = which( err > meanInValBracket[1] & err < meanInValBracket[2] )
    if ( modelFit$err < bestErr & length( inliers ) > minNtoFit  ) {
      bestErr = modelFit$err
      bestModel = modelFit
      fullInliers = inliers
      nInliers = length( fullInliers )
      }
    if ( verbose )
      print( paste( "It:", nIterations, "nIn:", length( inliers ),
        mean( err[ randSubset ] ), "v",  mean( err[ -randSubset ] ) ) )
    }
  # fit with full set
  if ( nInliers > 0 ) {
    finalFit = fitTransformToPairedPoints(   # step 5
      movingPoints[ fullInliers, ],
      fixedPoints[ fullInliers, ],
      transformType = transformType, lambda = lambda )
    } else {
      finalFit = NULL
      fullInliers = NULL
    }
  return(
    list(
      finalModel=finalFit,
      bestModel=bestModel,
      inliers = fullInliers ) )
}



#' Alternative random sample consensus (RANSACAlt)
#'
#' @param movingPoints moving points matrix
#' @param fixedPoints fixed points matrix
#' @param transformType Affine, Rigid and Similarity currently supported
#' @param nToTrim the number of points to throw away at each iteration
#' @param minProportionPoints the minimum proportion of points to return
#' @param nCVGroups number of cross-validation groups to determine error
#' @param lambda ridge penalty
#' @param domainImage image defining the domain for deformation maps.
#' @param numberOfFittingLevels integer specifying the number of fitting levels.
#' @param meshSize vector defining the mesh size at the initial fitting level.
#' @param dataWeights vector defining the individual weighting of the corresponding
#' scattered data value.  Default = NULL meaning all values are weighted the same.
#' @param verbose boolean
#'
#' @return output list contains best fitted model, inliers, outliers
#'
#' @export
RANSACAlt <- function(
  fixedPoints,
  movingPoints,
  transformType = "Affine",
  nToTrim = 2,
  minProportionPoints = 0.5,
  nCVGroups = 0,
  lambda = 1e-4,
  domainImage=NULL, numberOfFittingLevels=4, meshSize=1, dataWeights=NULL,
  verbose = FALSE ) {

#   1 Initialize the set S   with all points
#   2 Fit a line through the points in the set S
#   3 Calculate distance between line and each sample in S
#   4 Select the N    samples with the highest distance to the line, remove them from S
#   5 Repeat step 2-4 until an error criterion (e.g. sum of all squared distances between the line and points in S
#    S is below a certain threshold) is reached.

  myFP = fixedPoints
  myMP = movingPoints
  nMax = nrow( myFP )
  minn = round( minProportionPoints * nMax )
  its = 0
  rejectFixedPoints = NULL
  rejectMovingPoints = NULL
  bestErr = Inf
  while ( nMax >= minn  ) {
    modelFit = fitTransformToPairedPoints(   # step 2
      myMP,
      myFP,
      transformType = transformType, lambda = lambda,
      domainImage,
      numberOfFittingLevels,
      meshSize,
      dataWeights )
    mapComplement = applyAntsrTransformToPoint( modelFit$transform,
      myFP)
    err = sqrt( rowMeans( ( myMP - mapComplement )^2 ) )
    nToSelect = nMax - nToTrim
    inliers = sort( order( err )[1:nToSelect] )
    outliers = (1:nMax)[ -inliers ]
    if ( is.null( rejectFixedPoints ) ) {
      rejectFixedPoints=myFP[outliers,]
      rejectMovingPoints=myFP[outliers,]
    } else {
      rejectFixedPoints=rbind( rejectFixedPoints, myFP[outliers,] )
      rejectMovingPoints=rbind( rejectMovingPoints, myMP[outliers,] )
    }

    myFP = myFP[inliers,]
    myMP = myMP[inliers,]
    nMax = nrow( myFP )
    nCVGroups = 8
    if ( nCVGroups > 0 ) useCV = TRUE
    cvErr = rep( NA, nCVGroups )
    if ( useCV ) {
      cvgroups = sample(c(1:nCVGroups),nrow( myFP ),replace=T)
      group = 1
      for ( group in 1:nCVGroups) {
        if ( sum( cvgroups==group ) > 1 ) {
          tempfit = fitTransformToPairedPoints(
            myMP[cvgroups!=group,],
            myFP[cvgroups!=group,],
            transformType = transformType, lambda = lambda,
            domainImage,
            numberOfFittingLevels,
            meshSize,
            dataWeights )
          tempcv = applyAntsrTransformToPoint( tempfit$transform,
            myFP[cvgroups==group,])
          cvErr[group] = sqrt( rowMeans( ( myMP[cvgroups==group,] - tempcv )^2 ) )
          } else cvErr[group] = Inf
        }
      }
    if ( useCV ) meanErr = mean( cvErr, na.rm=TRUE ) else meanErr = mean( err )
    cvErr = mean( cvErr )
    isBest = FALSE
    if ( meanErr < bestErr ) {
      bestErr = meanErr
      bestFP = myFP
      bestMP = myMP
      bestRejectFixedPoints = rejectFixedPoints
      bestRejectMovingPoints = rejectMovingPoints
      isBest = TRUE
      }
    if ( verbose ) {
      print( paste( "It:", its,
        "mean(err):",mean(err),
        "nMax:",nMax,
        "cvErr:",cvErr,
        "isBest:",isBest ) )
      }
    its = its + 1
    }
  return(
    list(
      finalModel=modelFit,
      fixedPoints=bestFP,
      movingPoints = bestMP,
      rejectFixedPoints=bestRejectFixedPoints,
      rejectMovingPoints=bestRejectMovingPoints,
      bestErr=bestErr ) )
}





#' Convert match output to landmarks in physical space
#'
#' @param matchObject object, the output of \code{deepPatchMatch}
#' @param referenceImage the fixed image
#' @param movingImage the image that will be matched to the fixed image
#' @param patchSize size of patch features
#' @param whichK which matched point set (e.g. 1 gives the best, 2 second best and so on)
#' @return output list contains fixed and matched points
#'
#' @export
#' @examples
#' \dontrun{
#' library( keras )
#' library( ANTsR )
#' layout( matrix(1:2,nrow=1) )
#' nP1 = 50
#' nP2 = 200
#' psz = 32
#' img <- ri( 1 ) %>% iMath( "Normalize" )  %>% resampleImage( c( 2, 2 ) )
#' img2 <- ri( 2 ) %>% iMath( "Normalize" )  %>% resampleImage( c( 2, 2 ) )
#' mask = randomMask( getMask( img ), nP1 )
#' mask2 = randomMask( getMask( img2 ), nP2 )
#' matchO = deepPatchMatch( img2, img, mask, mask2 )
#' mlm = matchedLandmarks( matchO, img, img2, c(psz,psz) )
#' ct = 0
#' mxct = 18
#' lmImage1 = img * 0
#' lmImage2 = img2 * 0
#' for ( k in 1:nrow( mlm$fixedPoints ) ) {
#'   if ( ct < mxct ) {
#'     pt1i = makePointsImage( matrix(mlm$fixedPoints[k,],ncol=2), img, radius = 2 ) * k
#'     pt2i = makePointsImage( matrix(mlm$movingPoints[k,],ncol=2), img2, radius = 2 ) * k
#'     lmImage1 = lmImage1 + pt1i
#'     lmImage2 = lmImage2 + pt2i
#'     }
#'   ct = ct + 1
#'  }
#'
#' plot( img, lmImage1 )
#' plot( img2, lmImage2 )
#' }
matchedLandmarks <- function(
  matchObject,
  referenceImage,
  movingImage,
  patchSize,
  whichK = 1
) {
fixPoints = matrix( nrow = nrow( matchObject$matches ),
  ncol =  referenceImage@dimension )
matchedPoints = matrix( nrow = nrow( matchObject$matches ),
  ncol =  referenceImage@dimension )
off = 0 #
for ( k in 1:nrow( matchObject$matches ) ) {
  fpt =  matchObject$ffeatures$patchCoords[k,] + off
  pt1phys = antsTransformIndexToPhysicalPoint( referenceImage, fpt )
  fixPoints[k,] = pt1phys
  myselected = order( matchObject$costsK[k,])[ whichK ]
  if ( ! is.na( myselected ) )
  if ( ! is.na( matchObject$matchesK[k,myselected] ) ) {
#    pt1i = makePointsImage( pt1phys, img, radius = 2 ) * k
    fpt2 = matchObject$mfeatures$patchCoords[matchObject$matchesK[k,myselected],] + off
    pt2phys = antsTransformIndexToPhysicalPoint( movingImage, fpt2 )
    matchedPoints[k,] = pt2phys
#    pt2i = makePointsImage( pt2phys, img2, radius = 2 ) * k
#      lmImage1 = lmImage1 + pt1i
#      lmImage2 = lmImage2 + pt2i
#      print( paste( ct, k ) )
#      locimg1 = as.antsImage( matchObject$ffeatures$patches[k,,] )+127.5
#      locimg2 = as.antsImage(matchObject$mfeatures$patches[matchObject$matches[k,1],,])+127.5
    }
  }
nna=!is.na( matchedPoints[,1] )
return( list( fixedPoints=fixPoints[nna,], movingPoints=matchedPoints[nna,] ) )
}





#' Feature distance map
#'
#' @param image1 image 1
#' @param image2 image 2
#' @param jointMask mask covering features of both images
#' @param patchSize patchSize for both images
#' @param ... additional options to pass to deepFeatures e.g. vggmodel, block_name, substractor
#' @return error map
#'
#' @export
featureDistanceMap <- function( image1, image2, jointMask, patchSize=32, ... ) {
  df1=deepFeatures(image1,jointMask,patchSize=patchSize, ... )
  df2=deepFeatures(image2,jointMask,patchSize=patchSize, ... )
  err = rep( NA, nrow( df2$features ) )
  for ( k in 1:nrow(df2$features) ) {
    err[k] = mean(df2$features[k,]-df1$features[k,])
  }
  return( makeImage( jointMask, err ) )
}





#' Deep landmark regression stage
#'
#' @param model input deep model, presumably a unet
#' @param activation the activation function for the regression maps
#' @param theta the theta parameter for thresholded relu
#' @return the augmented model
#'
#' @export
deepLandmarkRegression <- function(
  model,
  activation = c("none","relu","trelu","softmax"),
  theta ) {
  if ( length( model$input_shape) == 5 ) {
    targetDimensionality = 3
    coordConv <- layer_input(  list( NULL, NULL, NULL, targetDimensionality ) )
  }
  if ( length( model$input_shape) == 4 ) {
    targetDimensionality = 2
    coordConv <- layer_input(  list( NULL, NULL, targetDimensionality ) )
  }
  nPoints = tail( unlist( model$output_shape ), 1 )
  # perform soft thresholding to get positive component of unet output
  if ( activation == 'none' ) {
    unet_output = unet$outputs[[1]]
  } else if ( activation == 'softmax') {
    unet_output = unet$outputs[[1]] %>%
      layer_activation_softmax()
  } else {
    if ( is.na( theta ) )
      unet_output = unet$outputs[[1]] %>%
        layer_activation_relu()
    if ( ! is.na( theta ) )
      unet_output = unet$outputs[[1]] %>%
        layer_activation_thresholded_relu( theta = theta )
    }
  targetDimensionality = length( coordConv )
  weightedRegressionList = tf$split( unet_output, as.integer(nPoints),
    axis=as.integer(targetDimensionality+1) )
  K <- keras::backend()
  regressAxes = as.integer(1:(targetDimensionality+1))
  regressAxes2 = as.integer(1:(targetDimensionality))
  regPoints = list()
  for ( k in 1:length(weightedRegressionList) ) {
    # forced averaging function
    weightedRegressionList[[k]] = weightedRegressionList[[k]] /
      ( K$sum( weightedRegressionList[[k]], regressAxes, keepdims = TRUE ) + 1e-19 )
    weightedRegressionList[[k]] =
      layer_multiply( list( coordConv, weightedRegressionList[[k]] ), trainable = FALSE )
    regPoints[[k]] = K$sum( weightedRegressionList[[k]], regressAxes2, keepdims = FALSE )
    }
  catout = layer_concatenate( regPoints, axis=1L ) %>%
    layer_reshape( c( nPoints , targetDimensionality ))
  keras_model( list( unet$inputs[[1]], coordConv), list( unet_output, catout  ) )
  # pp=predict( temp, list( images, coordconvbatch ) ) => gives landmark coordinates
}
