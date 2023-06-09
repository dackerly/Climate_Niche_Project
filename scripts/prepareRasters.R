# test script for testing functions
rm(list=ls())
library('mvtnorm')
library('terra')
source('scripts/climNiche2.R')

extFromVect <- function(x) {
  x1 <- min(geom(x)[,3])
  x2 <- max(geom(x)[,3])
  y1 <- min(geom(x)[,4])
  y2 <- max(geom(x)[,4])
  return(c(x1,x2,y1,y2))
}

aet <- rast('big_data/CHELSA_climate/CHELSA_annual_AET_1979_2013.tiff')
cwd <- rast('big_data/CHELSA_climate/CHELSA_annual_CWD_1979_2013.tiff')
tmin <- rast('big_data/CHELSA_climate/CHELSA_tmin10_01_1979-2013_V1.2_land.tif')
plot(tmin)
ext(aet)

CA <- vect("data/gis_data/LatLong/CA.shp")
CAext <- ext(extFromVect(CA))

crs(cwd)
crs(CA) <- crs(cwd)

tmn <- crop(tmin,CAext)
tmn <- mask(tmn,CA)
tmn <- tmn/10
plot(tmn)
plot(CA,add=T)

cwd <- crop(cwd,CAext)
cwd <- mask(cwd,CA)
plot(cwd)
plot(CA,add=T)

aet <- crop(aet,CAext)
aet <- mask(aet,CA)
plot(aet)
plot(CA,add=T)

writeRaster(cwd,'data/gis_data/CAcwd.tiff')
writeRaster(aet,'data/gis_data/CAaet.tiff')
writeRaster(tmn,'data/gis_data/CAtmn.tiff')
