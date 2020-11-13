Sys.setenv("TF_NUM_INTEROP_THREADS"=16)
Sys.setenv("TF_NUM_INTRAOP_THREADS"=16)
set.seed( Sys.time())
library( ANTsR )
library( patchMatchR )
library( ANTsRNet )
library( keras )
# these images are in inst extdata
fn1=system.file("extdata", "1000_3.nii.gz", package = "patchMatchR")
fn2=system.file("extdata", "1017_3.nii.gz", package = "patchMatchR")
rdims = rep( 2, 3 )
img <- antsImageRead(fn1) %>% iMath( "Normalize" ) %>% resampleImage( rdims )
img2 <- antsImageRead(fn2) %>% iMath( "Normalize" ) %>% resampleImage( rdims )
trans=antsRegistration( img, img2, "Translation")
img2 = antsApplyTransforms( img, img2, trans$fwdtransforms, interpolator='nearestNeighbor')

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
img2 = randomRotateImage(img2,sdAff=0.25)
plot(img,img2,alpha=0.5)
# img2[55:128,1:128,1:144]=0
# img2[1:128,1:128,55:144]=0
reg = antsRegistration( img, img2, "Affine",
  affIterations=c( 1000, 1000, 100, 11 ) )
fm1=iMath(img*255, "Canny", 3.0, 18, 25 )
fm2=iMath(img2*255, "Canny", 3.0, 18, 25  )
fmdl = createResNetModel3D( list(NULL,NULL,NULL,1), lowestResolution=4,
  mode='regression', numberOfClassificationLabels = 64 )
fmdl = keras_model( fmdl$inputs, fmdl$layers[[38]]$output)
nP1 = 1024*2
nP2 = nP1 # round( sum( fm2 == 1 )*0.5 )
print( paste(nP1,nP2))
mask = randomMask( fm1, nP1 )
mask2 = randomMask( fm2, nP2 )
psz = 16
print("Begin")
match = deepPatchMatch( img2, img, mask2, mask, psz, psz,
  vggmodel=fmdl, subtractor=0.5 )
#
print("Done")
mlm = matchedLandmarks( match, img, img2, c(psz,psz) )
tx = fitTransformToPairedPoints( mlm$movingPoints, mlm$fixedPoints,"Affine" )
warped = applyAntsrTransformToImage( tx$transform, img2, img )
bestMI = 0
while( nrow( mlm$fixedPoints ) > 21 ) {
  rns = RANSACAlt( mlm$fixedPoints, mlm$movingPoints, transformType = "Affine",
    minProportionPoints=0.95, nToTrim = 5, verbose = F )
  temp = applyAntsrTransform( rns$finalModel$transform,
    antsImageClone(img2), antsImageClone(img)  )
  print( paste(
    "(org):",antsImageMutualInformation(img,img2),
    "(aff):",antsImageMutualInformation(img,reg$warpedmovout),
    "(all):",antsImageMutualInformation(img,warped),
    "(rsc):",antsImageMutualInformation(img,temp),
     "(npt):",nrow( mlm$fixedPoints ) ) )
  locmi = antsImageMutualInformation(antsImageClone(img),antsImageClone(temp))
  if ( locmi < bestMI ) {
    mmm = antsImageClone(img)*0+1
    bestMI = locmi
    bestPoints = rns
    }
    if ( FALSE ) {
      mpi1 = makePointsImage( rns$fixedPoints, mmm, radius=3)
      mpi2 = makePointsImage( rns$movingPoints, mmm, radius=3)
      mpi1r = makePointsImage( rns$rejectFixedPoints, mmm, radius=3)
      mpi2r = makePointsImage( rns$rejectMovingPoints, mmm, radius=3)
      layout( matrix(1:5,nrow=1) )
      plot(antsImageClone(img)*11,mpi1,slice=60)
      plot(antsImageClone(img2)*11,mpi2,slice=60)
      plot(antsImageClone(img)*11,temp,alpha=0.5,slice=60)
      plot(antsImageClone(img)*11,mpi1r,slice=60)
      plot(antsImageClone(img2)*11,mpi2r,slice=60)
    }
  mlm=rns
}

temp = applyAntsrTransform( bestPoints$finalModel$transform,
  antsImageClone(img2), antsImageClone(img)  )
mpi1 = makePointsImage( bestPoints$fixedPoints, mmm, radius=3)
mpi2 = makePointsImage( bestPoints$movingPoints, mmm, radius=3)
layout( matrix(1:4,nrow=1) )
plot(antsImageClone(img)*11,mpi1,slice=60)
plot(antsImageClone(img2)*11,mpi2,slice=60)
plot(antsImageClone(img)*11,temp*11,alpha=0.5,slice=60)
plot(antsImageClone(img)*11,reg$warpedmovout*11,alpha=0.5,slice=60)
antsImageWrite( img, '/tmp/a1.nii.gz' )
antsImageWrite( mpi1, '/tmp/a2.nii.gz' )
antsImageWrite( mpi2, '/tmp/b2.nii.gz' )
antsImageWrite( img2, '/tmp/b1.nii.gz' )
antsImageWrite( temp, '/tmp/w.nii.gz' )
