library( keras )
library( ANTsR )
library( patchMatchR )
idim = 2 # image dimensionality
nP1 = 30
nP2 = 300
psz = 32
img1 <- ri( 1 ) %>% iMath( "Normalize" )
img2 <- ri( 2 ) %>% iMath( "Normalize" )
img2 = antsRegistration( img1, img2, "Rigid" )$warpedmovout
scaleParam = 2.0
fullMask1 = iMath(img1*2000,"Canny",scaleParam,9,10)
fullMask2 = iMath(img2*2000,"Canny",scaleParam,9,10)
mask1 = randomMask( fullMask1, nP1 )
mask2 = randomMask( fullMask2, nP2 )
matchO = deepPatchMatch(
  img2, img1,
  mask2, mask1, block_name = 'block2_conv2',  knn = 1 ) # knnSpatial=50 )
mlm = matchedLandmarks( matchO, img1, img2, rep(psz, idim) )
k = 1:nrow( mlm$fixedPoints)
lmImage1 = makePointsImage(
  matrix(mlm$fixedPoints[k,],ncol=idim), img1, radius = 2 )
lmImage2 = makePointsImage(
  matrix(mlm$movingPoints[k,],ncol=idim), img2, radius = 2 )
layout( matrix(1:2,nrow=1) )
plot( img1*222, lmImage1, doCropping=T  )
plot( img2*222, lmImage2, doCropping=T   )
