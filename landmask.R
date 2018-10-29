#!/usr/local/bin/R

## define the size of the mask
n <- 361
## define the data filename
f <- "landmask.dat"

# https://gis.stackexchange.com/questions/75033/given-a-lat-and-lon-identify-if-the-point-is-over-land-or-ocean

library(maptools)
data(wrld_simpl)
gpclibPermit()


## Create a SpatialPoints object
lat <- seq( -90,  90, length=n/2)
lon <- seq(-180, 180, length=n  )
points <- expand.grid(lon, lat)  
pts <- SpatialPoints(points, proj4string=CRS(proj4string(wrld_simpl)))

## Find which points fall over land
ii <- !is.na(over(pts, wrld_simpl)$FIPS)

## Check that it worked
#plot(wrld_simpl)
#points(pts, col=1+ii, pch=16)

write.table(cbind(points,matrix(ii,ncol=1)), file=f,row.names=FALSE,col.names=FALSE)
