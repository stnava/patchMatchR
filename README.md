# patchMatchR
very basic patch matching functionality based on ANTsR -- see [documentation](https://stnava.github.io/patchMatchR/)


[![Build Status](https://travis-ci.org/stnava/patchMatchR.png?branch=master)](https://travis-ci.org/stnava/patchMatchR)

### installation

```
devtools::install_github( "stnava/patchMatchR")
```

or if you have ANTsR etc installed then clone the repo and do:

```
cd patchMatchR
R CMD INSTALL .
```

### getting started

in `R`

2D

```
library( ANTsR )
library( patchMatchR )
img <- ri( 1 )  %>% iMath( "Normalize" )
img2 <- ri( 2 )%>% resampleImage( 2 ) %>% iMath( "Normalize" )
mask = randomMask( getMask( img ), 10 )
match = patchMatch( img2, img, mask, visualize=T, verbose=T )
myMatches = matchedPatches( img2, img, match, verbose = TRUE )
plot( myMatches[[1]][[4]] * 100 )
plot( myMatches[[2]][[4]] * 100 )
```

2D with pre-existing mapping

```
library( ANTsR )
library( patchMatchR )
img <- ri( 1 )  %>% iMath( "Normalize" )
img2 <- ri( 2 )%>% resampleImage( 2 ) %>% iMath( "Normalize" )
ireg = antsRegistration( img, img2, "SyNCC" )
mask = randomMask( getMask( img ), 10 )
match = patchMatch( img2, img, mask, visualize=T, verbose=T,
  initialMap = ireg )
myMatches = matchedPatches( img2, img, match, verbose = TRUE )
plot( myMatches[[1]][[4]] * 100 )
plot( myMatches[[2]][[4]] * 100 )
```



3D

```
library( ANTsR )
library( patchMatchR )
img <- antsImageRead( getANTsRData( "ch2" ) ) %>% resampleImage( 1 ) %>%
  iMath( "Normalize" )
img2 <- resampleImage( img, 2 ) %>% iMath( "Normalize" )
mask = randomMask( getMask( img ), 25 )
myr = 15
match = patchMatch( img2, img, mask, visualize=T, verbose=T,
  fixedPatchRadius = myr )
myMatches = matchedPatches( img2, img, match, verbose = TRUE,
  fixedPatchRadius = myr+1  )
plot( myMatches[[1]][[4]] * 100 )
plot( myMatches[[2]][[4]] * 100 )
```

TODO


- allow a probabilistic representation of matching

- provide a method for generating transforms from the matched landmarks (my RANSAC implementation is weak)

- test on GPU

- compare to SIFT

- pass a transform to deepPatchMatch s.t. moving point distance can be calculated wrt tx

