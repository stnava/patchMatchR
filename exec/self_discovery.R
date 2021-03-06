library( keras )
library( ANTsR )
library( patchMatchR )
dev.new(width=6,height=4)
layout( matrix(1:6,nrow=2))
# build a graph reprsentation from a single image
idim = 2
# count the number of neighbors with correlations > thresh
countcors <- function( x, thresh=0.5 ) {
  mycounters = rep( NA, nrow( x  ) )
  for ( k in 1:nrow( x ) )
    mycounters[k] = as.numeric( table( x[k,] > thresh ) )[2]
  mycounters
}
gaussKernel <- function (X = NULL, sigma = NULL) {
  return(exp(-1 * as.matrix(dist(X)^2)/sigma))
}

mysink <- function( x, thresh, iterations = 3 ) {
  diag(x) = 0
  x[ x < thresh ] = 0
  for ( i in 1:iterations ) {
    for ( j in 1:nrow( x ) ) {
      mysum = sum(x[j,],na.rm=T)
      if ( mysum == 0 ) mysum = 1
      x[j,] = x[j,]/mysum
    }
    for ( j in 1:ncol( x ) ) {
      mysum = sum(x[,j],na.rm=T)
      if ( mysum == 0 ) mysum = 1
      x[,j] = x[,j]/mysum
    }
  }
  x
}
set.seed( Sys.time() )
if ( ! exists( "ss" ) ) ss = sample(1:6)[1]
print(ss)
selfPoints = list()
for ( ss in 1:2 ) {
  img2 = ri( ss )
  scaleParam = 4.5
  fullMask2 = getMask( img2 )
#  fullMask2 = iMath(img2*2000,"Canny",scaleParam,8,10)
  npts = 5000
  mask2 = randomMask( fullMask2, npts )
  patchSize = 32
  patchSizeDivBy2 = patchSize/2
  myFeats = deepFeatures( img2, mask2, patchSize = patchSize  )
#  t1=Sys.time()
  mycor = cor( t(myFeats$features ) )
#  t2=Sys.time()
#  mycor2 = gaussKernel( myFeats$features, 2000 )
#  t3=Sys.time()
  myut = upper.tri( mycor )
  useSink = TRUE
  if ( useSink ) mycor = mysink( mycor, 0.0 )
  mycounts = countcors( mycor, thresh = quantile( mycor[myut], 0.9 ) )
  # bestk = sort( mycounts )[ round( npts * 0.01 )]
  bestk = sort( mycounts )[ 500 ] # N-th index
  goodones = which( mycounts <= bestk )
  print(length(goodones))
  uniquePoints = makePointsImage(
    matrix(myFeats$patchCoords[goodones,]+patchSizeDivBy2,ncol=idim), img2, radius = 1 )
  plot( img2, uniquePoints )
  selfPoints[[ ss ]] = uniquePoints
  }
#############################
message("first pass for unique LMs")
message("now pass the self-aware LMs to match images")
mask1 = thresholdImage( selfPoints[[1]], 1, Inf )
mask2 = thresholdImage( selfPoints[[2]], 1, Inf )
img1 = ri( 1 )
img2 = ri( 2 )
myknn = 1
matchO = deepPatchMatch(
  img2, img1, knnSpatial = 50,
  mask2, mask1, block_name = 'block2_conv2',  knn = myknn )
mlm = matchedLandmarks( matchO, img1, img2, rep(patchSize, idim) )
subsam = sample( 1:nrow(mlm$fixedPoints), 10 )
lmImage1 = makePointsImage(
  matrix(mlm$fixedPoints[subsam,],ncol=idim), img1, radius = 2 )
lmImage2 = makePointsImage(
  matrix(mlm$movingPoints[subsam,],ncol=idim), img2, radius = 2 )
layout( matrix(1:2,nrow=1))
plot( img1*222, lmImage1, doCropping=T  )
plot( img2*222, lmImage2, doCropping=T   )

# sdmat = sparseDistanceMatrix( t(myFeats$features), k = 25,  sinkhorn = FALSE )
message("represent as network")
# library( network )
# net = network( mycor )
# plot( net )
