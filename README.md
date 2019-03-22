# patchMatchR
very basic patch matching functionality based on ANTsR

### installation

```
devtools::install_github( "stnava/patchMatchR")
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
```


