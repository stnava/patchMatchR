library( keras )
library( ANTsR )
library( patchMatchR )
layout( matrix(1:2,nrow=1) )
idim = 2 # image dimensionality
nP1 = 30
nP2 = 300
psz = 32
img1 <- ri( 1 ) %>% iMath( "Normalize" )
img2 <- ri( 2 ) %>% iMath( "Normalize" )
fn1 = system.file("extdata", "fish0709.png", package = "patchMatchR")
fn2 = system.file("extdata", "fish0724.png", package = "patchMatchR")
img1=antsImageRead( fn1 ) %>% iMath("PadImage",20)
img2=antsImageRead( fn2 )
img2 = antsRegistration( img1, img2, "Rigid" )$warpedmovout %>%
  thresholdImage( 0.33, Inf )
fullMask1 = getMask( img1 )
fullMask2 = getMask( img2 )
bordMask1 = fullMask1 - iMath(fullMask1,"ME",1)
bordMask2 = fullMask2 - iMath(fullMask2,"ME",1)
mask1 = randomMask( bordMask1, nP1 )
mask2 = randomMask( bordMask2, nP2 )
matchO = deepPatchMatch(
#  iMath(bordMask2,"D"), iMath(bordMask1,"D"),
 smoothImage(img2,0.5), smoothImage(img1,0.5),
 mask2, mask1, block_name = 'block2_conv2',  knn = 1, knnSpatial=50 )
mlm = matchedLandmarks( matchO, img1, img2, rep(psz, idim) )
k = 1:nrow( mlm$fixedPoints)
lmImage1 = makePointsImage(
  matrix(mlm$fixedPoints[k,],ncol=idim), img1, radius = 2 )
lmImage2 = makePointsImage(
  matrix(mlm$movingPoints[k,],ncol=idim), img2, radius = 2 )
plot( img1*222, lmImage1, doCropping=F  )
plot( img2*222, lmImage2, doCropping=F   )
