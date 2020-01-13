library( keras )
library( ANTsR )
library( patchMatchR )
set.seed( 1 ) # Sys.time() )
idim = 2 # image dimensionality
nP1 = 2
nP2 = 500
psz = 32
myknn = 12 # how many points to find in target image
img1 <- ri( 1 ) %>% iMath( "Normalize" )
img2 <- ri( 6 ) %>% iMath( "Normalize" )
img2 = antsRegistration( img1, img2, "Rigid" )$warpedmovout
scaleParam = 4.5
fullMask1 = iMath(img1*2000,"Canny",scaleParam,8,10) * getMask( img1 )
mask1 = randomMask( fullMask1, nP1 )
fullMask2 = iMath(img2*2000,"Canny",scaleParam,8,10)
mask2 = randomMask( fullMask2, nP2 )
dev.new(width=6,height=2)
layout( matrix(1:3,nrow=1, byrow=T ))
plot( fullMask2 )
matchO = deepPatchMatch(
  img2, img1,
  mask2, mask1, block_name = 'block2_conv2',  knn = myknn , knnSpatial=25 )
mlm = matchedLandmarks( matchO, img1, img2, rep(psz, idim) )
lmImage1 = makePointsImage(
  matrix(mlm$fixedPoints,ncol=idim), img1, radius = 2 )
lmImage2 = lmImage1 * 0
for ( k in 1:myknn ) {
  mlm2 = matchedLandmarks( matchO, img1, img2, rep(psz, idim), whichK = k )
  tempmat = matrix(mlm2$movingPoints,ncol=idim)
  if ( nrow( tempmat ) > 0 )
    lmImage2 = lmImage2 + makePointsImage( tempmat, img2, radius = 2 )
  }
plot( img1*222, lmImage1, doCropping=T  )
plot( img2*222, lmImage2, doCropping=T   )
