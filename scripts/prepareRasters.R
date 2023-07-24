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
tav <- rast('big_data/CHELSA_climate/CHELSA_bio10_01.tiff')
ppt <- rast('big_data/CHELSA_climate/CHELSA_bio10_12.tiff')
plot(tav)
ext(aet)

CA <- vect("data/gis_data/LatLong/CA.shp")
CAext <- ext(extFromVect(CA))

crs(cwd)
crs(CA) <- crs(cwd)

tav <- crop(tav,CAext)
tav <- mask(tav,CA)
tav <- tav/10
plot(tav)
plot(CA,add=T)

cwd <- crop(cwd,CAext)
cwd <- mask(cwd,CA)
plot(cwd)
plot(CA,add=T)

aet <- crop(aet,CAext)
aet <- mask(aet,CA)
plot(aet)
plot(CA,add=T)

ppt <- crop(ppt,CAext)
ppt <- mask(ppt,CA)
plot(ppt)
plot(CA,add=T)

writeRaster(cwd,'data/gis_data/CAcwd.tiff',overwrite=T)
writeRaster(aet,'data/gis_data/CAaet.tiff',overwrite=T)
writeRaster(tav,'data/gis_data/CAtav.tiff',overwrite=T)
writeRaster(ppt,'data/gis_data/CAppt.tiff',overwrite=T)

CAs <- c(tav,cwd,aet,ppt)
CAv <- data.frame(values(CAs))
names(CAv) <- c('tav','cwd','aet','ppt')
sel <- which(complete.cases(CAv))
CAv <- CAv[sel,]
head(CAv)

rsamp <- sample(1:nrow(CAv),1000)
pairs(CAv[rsamp,])
