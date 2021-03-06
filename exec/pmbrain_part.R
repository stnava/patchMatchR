Sys.setenv("TF_NUM_INTEROP_THREADS"=28)
Sys.setenv("TF_NUM_INTRAOP_THREADS"=28)
set.seed(Sys.time())
library( tensorflow )
library( keras )
library( ANTsR )
library( ANTsRNet )
library( abind )
library( patchMatchR)
message("If you get a weird error then run it again - h5py load issue")
# these images are in inst extdata
fn1=system.file("extdata", "1000_3.nii.gz", package = "patchMatchR")
fn2=system.file("extdata", "1017_3.nii.gz", package = "patchMatchR")
rdims = rep( 2, 3 )
img <- antsImageRead(fn1) %>% iMath( "Normalize" ) %>% resampleImage( rdims )
img2 <- antsImageRead(fn2) %>% iMath( "Normalize" ) %>% resampleImage( rdims )
# take care of physical space differences here - or could use registration
img2 = resampleImageToTarget( img2, img, interpType='nearestNeighbor')

polarX <- function(X) {
        x_svd <- svd(X)
        P <- x_svd$u %*% diag(x_svd$d) %*% t(x_svd$u)
        Z <- x_svd$u %*% t(x_svd$v)
        if (det(Z) < 0)
            Z = Z * (-1)
        return(list(P = P, Z = Z, Xtilde = P %*% Z))
    }
randAff <- function( loctx,  txtype = "Rigid", sdAffine,
      idparams, fixParams ) {
      idim = 3
      noisemat = stats::rnorm(length(idparams), mean = 0, sd = sdAffine)
      if (txtype == "Translation")
        noisemat[1:(length(idparams) - idim )] = 0
      idparams = idparams + noisemat
      idmat = matrix(idparams[1:(length(idparams) - idim )],
                  ncol = idim )
      idmat = polarX(idmat)
              if (txtype == "Rigid")
                  idmat = idmat$Z
              if (txtype == "Affine")
                  idmat = idmat$Xtilde
              if (txtype == "ScaleShear")
                  idmat = idmat$P
              idparams[1:(length(idparams) - idim )] = as.numeric(idmat)
      setAntsrTransformParameters(loctx, idparams)
      setAntsrTransformFixedParameters( loctx, fixParams )
      return(loctx)
      }
################################################################################
randomRotateImage <- function( image, sdAff=0.1 ) {
  loctx = randAff( loctx, sdAffine=sdAff, txtype = 'Rigid',
    idparams = idparams, fixParams = fixedParams)
  imageR = applyAntsrTransformToImage( loctx, image, image,
      interpolation = "nearestNeighbor" )
  return( imageR )
}
fixedParams = getCenterOfMass( img2 * 0 + 1 )
loctx <- createAntsrTransform(precision = "float",
  type = "AffineTransform", dimension = 3  )
setAntsrTransformFixedParameters(loctx, fixedParams)
idparams = getAntsrTransformParameters( loctx )
setAntsrTransformParameters( loctx, idparams )
setAntsrTransformFixedParameters(loctx, fixedParams)
# img[55:128,1:128,1:144]=0 # half-brain coverage
# img[1:128,1:128,55:144]=0 #
img[ round(dim(img)[1]*0.4):dim(img)[1], 1:dim(img)[2], 1:dim(img)[3]  ] = 0 # half-brain coverage
# img[ 1:128,  1:128, 77:144 ] = 0 # quarter-brain if you do this
# img2 = randomRotateImage( img2, sdAff=0.33 )
layout(matrix(1:2,nrow=1))
plot(img*255,slice=60)
plot(img2*255,slice=50)
fm1=iMath(img*255, "Canny", 3.0, 18, 25 )
fm2=iMath(img2*255, "Canny", 3.0, 18, 25  )
################################################################################
nP1 = 500
nP2 = 2400 # round( sum( fm2 == 1 )*0.5 )
print( paste(nP1,nP2))
mask = randomMask( fm1, nP1 )
mask2 = randomMask( fm2, nP2 )
psz = 20
print(paste("Begin",Sys.time()))
match = deepPatchMatch( img2, img, mask2, mask, psz, psz )
print(paste("Done",Sys.time()))
mlm = matchedLandmarks( match, img, img2, c(psz,psz) )
tx = fitTransformToPairedPoints( mlm$movingPoints, mlm$fixedPoints,"Affine" )
warped = applyAntsrTransformToImage( tx$transform, img2, img )
################################################################################
rns = RANSACAlt( mlm$fixedPoints, mlm$movingPoints, transformType = "Affine",
    minProportionPoints=0.25, nToTrim = 4, nCVGroups=10, verbose = FALSE )
temp = applyAntsrTransform( rns$finalModel$transform,
    antsImageClone(img2), antsImageClone(img)  )
locmi = antsImageMutualInformation(antsImageClone(img),antsImageClone(temp))
print( paste(
    "(org):",antsImageMutualInformation(img,img2),
#    "(aff):",antsImageMutualInformation(img,reg$warpedmovout),
    "(all):",antsImageMutualInformation(img,warped),
    "(rsc):",locmi,
     "(npt):",nrow( rns$fixedPoints ) ) )

layout( matrix(1:2,nrow=1) )
plot(img2*255,antsImageClone(img)*255,alpha=0.5,axis=2,nslices=21,ncolumns=7)
plot(temp*255,antsImageClone(img)*255,alpha=0.5,axis=2,nslices=21,ncolumns=7)

feat1 = deepFeatures( img, getMask(img), patchSize=20 )
feat2 = deepFeatures( img2, getMask(img2), patchSize=20 )
derka
# jointMask = thresholdImage( getMask( img ) + getMask( img2 ), 1, 2 )
# fdm = featureDistanceMap(img, img2, jointMask, patchSize = 12 )


if ( FALSE ) {

  mmm = img * 0 + 1
  mmm2 = img2 * 0 + 1
  mpi1 = makePointsImage( rns$fixedPoints, mmm, radius=3)
  mpi2 = makePointsImage( rns$movingPoints, mmm2, radius=3)
  mpi1r = makePointsImage( rns$rejectFixedPoints, mmm, radius=3)
  mpi2r = makePointsImage( rns$rejectMovingPoints, mmm2, radius=3)

  layout( matrix(1:5,nrow=1) )
  plot(antsImageClone(img)*11,mpi1,slice=60)
  plot(antsImageClone(img2)*11,mpi2,slice=60)
  plot(antsImageClone(img)*11,temp,alpha=0.5,slice=60)
  plot(antsImageClone(img)*11,mpi1r,slice=60)
  plot(antsImageClone(img2)*11,mpi2r,slice=60)

  antsImageWrite( img, '/tmp/a1.nii.gz' )
  antsImageWrite( mpi1, '/tmp/a2.nii.gz' )
  antsImageWrite( mpi2, '/tmp/b2.nii.gz' )
  antsImageWrite( img2, '/tmp/b1.nii.gz' )
  antsImageWrite( temp, '/tmp/w.nii.gz' )
  }
