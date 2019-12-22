library( keras )
library( ANTsR )
library( patchMatchR )
dev.new(width=12,height=8)
layout( matrix(1:6,nrow=2))
# build a graph reprsentation from a single image
idim = 2
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
if ( ! exists( "ss" ) ) ss = sample(1:6)[1]
print(ss)
for ( ss in 1:6 ) {
  img2 = ri( ss )
  scaleParam = 4.5
  fullMask2 = getMask( img2 ) # iMath(img2*2000,"Canny",scaleParam,8,10)
  fullMask2 = iMath(img2*2000,"Canny",scaleParam,8,10)
  npts = 2000
  mask2 = randomMask( fullMask2, npts )
  patchSize = 32
  patchSizeDivBy2 = patchSize/2
  myFeats = deepFeatures( img2, mask2, patchSize = patchSize  )
  mycor = cor( t(myFeats$features ) )
  mycor = mysink( mycor, 0.2 )
  mycounts = countcors( mycor, thresh = quantile(mycor,0.95) )
  # bestk = sort( mycounts )[ round( npts * 0.01 )]
  bestk = sort( mycounts )[ 50 ]
  goodones = which( mycounts <= bestk )
  length(goodones)
  uniquePoints = makePointsImage(
    matrix(myFeats$patchCoords[goodones,]+patchSizeDivBy2,ncol=idim), img2, radius = 3 )
  plot( img2, uniquePoints )
  }
#############################
stop("first pass for unique LMs")
# sdmat = sparseDistanceMatrix( t(myFeats$features), k = 25,  sinkhorn = FALSE )
library( network )
net = network( mycor )
plot( net )
