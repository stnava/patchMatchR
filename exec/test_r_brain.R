library( keras )
library( ANTsR )
library( patchMatchR )
layout( matrix(1:3,nrow=1) )
idim = 2 # image dimensionality
nP1 = 2
nP2 = 3000
psz = 40
img1 <- ri( 2 ) %>% iMath( "Normalize" )
img2 <- ri( 4 ) %>% iMath( "Normalize" )
img2 = antsRegistration( img1, img2, "Rigid" )$warpedmovout
scaleParam = 2.5
fullMask1 = iMath(img1*2000,"Canny",scaleParam,8,10)
plot( fullMask1 )
fullMask2 = iMath(img2*2000,"Canny",scaleParam,8,10)
mask1 = randomMask( fullMask1, nP1 )
mask2 = randomMask( fullMask2, nP2 )
myknn = 10
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
