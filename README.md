# patchMatchR
very basic patch matching functionality based on ANTsR

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

```
library( ANTsR )
library( patchMatchR )
img <- ri( 1 )  %>% iMath( "Normalize" )
img2 <- ri( 2 )%>% resampleImage( 2 ) %>% iMath( "Normalize" )
mask = randomMask( getMask( img ), 10 )
match = patchMatch( img2, img, mask, visualize=T, verbose=T )
myMatches = matchedPatches( img2, img, match, verbose = TRUE )
plot( myMatches[[1]][[4]] )
plot( myMatches[[2]][[4]] )
```


