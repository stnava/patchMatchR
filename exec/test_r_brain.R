library( keras )
library( ANTsR )
library( patchMatchR )
set.seed( Sys.time() )
layout( matrix(1:6,nrow=2, byrow=T ) )
idim = 2 # image dimensionality
nP1 = 2
nP2 = 3000
psz = 40
myknn = 10 # how many points to find in target image
img1 <- ri( 2 ) %>% iMath( "Normalize" )
img2 <- ri( 2 ) %>% iMath( "Normalize" )
img3 <- ri( 4 ) %>% iMath( "Normalize" )
img2 = antsRegistration( img1, img2, "Rigid" )$warpedmovout
img3 = antsRegistration( img1, img3, "Rigid" )$warpedmovout
scaleParam = 4.5
fullMask1 = iMath(img1*2000,"Canny",scaleParam,8,10) * getMask( img1 )
mask1 = randomMask( fullMask1, nP1 )
fullMask2 = iMath(img2*2000,"Canny",scaleParam,8,10)
mask2 = randomMask( fullMask2, nP2 )
plot( fullMask2 )
fullMask3 = iMath(img3*2000,"Canny",scaleParam,8,10)
mask3 = randomMask( fullMask3, nP2 )
matchO = deepPatchMatch(
  img2, img1,
  mask2, mask1, block_name = 'block2_conv2',  knn = myknn ) # knnSpatial=50 )
mlm = matchedLandmarks( matchO, img1, img2, rep(psz, idim) )
lmImage1 = makePointsImage(
  matrix(mlm$fixedPoints,ncol=idim), img1, radius = 2 )
lmImage2 = lmImage1 * 0
for ( k in 1:myknn ) {
  mlm2 = matchedLandmarks( matchO, img1, img2, rep(psz, idim), whichK = k )
  lmImage2 = lmImage2 +
    makePointsImage( matrix(mlm2$movingPoints,ncol=idim), img2, radius = 2 )
  }
plot( img1*222, lmImage1, doCropping=T  )
plot( img2*222, lmImage2, doCropping=T   )

#### now do the next example
matchO = deepPatchMatch(
  img3, img1,
  mask3, mask1, block_name = 'block2_conv2',  knn = myknn ) # knnSpatial=50 )
mlm = matchedLandmarks( matchO, img1, img3, rep(psz, idim) )
lmImage1 = makePointsImage(
  matrix(mlm$fixedPoints,ncol=idim), img1, radius = 2 )
lmImage3 = lmImage1 * 0
for ( k in 1:myknn ) {
  mlm2 = matchedLandmarks( matchO, img1, img3, rep(psz, idim), whichK = k )
  lmImage3 = lmImage3 +
    makePointsImage( matrix(mlm2$movingPoints,ncol=idim), img3, radius = 2 )
  }
plot( fullMask3 )
plot( img1*222, lmImage1, doCropping=T  )
plot( img3*222, lmImage3, doCropping=T   )

stop("Look past this point for a network representation")
# build a graph reprsentation from a single image
countcors <- function( x, thresh=0.5 ) {
  mycounters = rep( NA, nrow( x  ) )
  for ( k in 1:nrow( x ) )
    mycounters[k] = as.numeric( table( x[k,] > thresh ) )[2]
  mycounters
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
ss = sample(1:6)[1]
print(ss)
img2 = ri( ss )
scaleParam = 4.5
fullMask2 = getMask( img2 ) # iMath(img2*2000,"Canny",scaleParam,8,10)
npts = 2000
mask2 = randomMask( fullMask2, npts )
myFeats = deepFeatures( img2, mask2, patchSize = 32  )
mycor = cor( t(myFeats$features ) )
mycor = mysink( mycor, 0.2 )
mycounts = countcors( mycor, thresh = quantile(mycor,0.95) )
# bestk = sort( mycounts )[ round( npts * 0.01 )]
bestk = sort( mycounts )[ 50 ]
goodones = which( mycounts <= bestk )
length(goodones)
patchSizeDivBy2 = 16
uniquePoints = makePointsImage(
  matrix(myFeats$patchCoords[goodones,]+patchSizeDivBy2,ncol=idim), img2, radius = 3 )
plot( img2, uniquePoints )
#############################
stop("first pass for unique LMs")
kk = mysink( mycor, 0.2 )
hist( kk[kk>0] )
table( kk > 0.8 )
cth = 0.98
for ( k in 1:10 ) {
  kk = mysink( kk, quantile(kk,cth) )
  hist( kk[kk>0] )
  print(  quantile(kk,cth) )
  }

# sdmat = sparseDistanceMatrix( t(myFeats$features), k = 25,  sinkhorn = FALSE )
library( network )
net = network( sdmat )
plot( net )
