# test script for testing functions
rm(list=ls())
library('mvtnorm')
library('terra')
source('scripts/climNiche2.R')

aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
CA <- vect("data/gis_data/LatLong/CA.shp")
plot(cwd)
plot(CA,add=T)

rangeExtend <- function(x,extAmount=0.5) {
  rx <- range(x)
  rx[1] <- rx[1]-extAmount
  rx[2] <- rx[2]+extAmount
  return(rx)
}

cd <- read.csv('data/test_data/test_climDat.csv',row.names = 1)
head(cd)
table(cd$scientificName)

op <- par(mfrow=c(1,2))
plot(aet)
points(cd[,c('longitude','latitude')])
plot(cd$cwd,cd$aet)
par(op)

## see climNiche.R for specification of the input format
cn <- climNiche(cd)

cn$climStats
cn$covMatrix
cn$pcResults
cn$pcResults$loadings
head(cn$climData)

### now the second algorithm, with means weighted by the available environment; and the climate data is re-extraxted so this can be run with cd dataframe providing only lat/long
# clip rasters to domain
st <- c(aet,cwd);names(st) <- c('aet','cwd')
st <- crop(st,ext(c(rangeExtend(cd$longitude),rangeExtend(cd$latitude))))

cn <- climNiche(cd,st,onePerCell = T,trunc=0.01)
cn$climStats

# read in a color ramp
ramp <- read.csv('data/colorRamps/jja6ColRamp.csv',as.is=T)
nrow(ramp)
head(ramp)
plot(1:nrow(ramp),1:nrow(ramp),pch=19,col=ramp$Hex)

## now plot again with colors proportional to normalized multivariate probability
pVals <- as.numeric(cut(cn$climData$mnScaledProb,nrow(ramp)))
plot(cwd,xlim=rangeExtend(cd$longitude),ylim=rangeExtend(cd$latitude))
points(cd[,c('longitude','latitude')],pch=19,col=ramp$Hex[pVals])

## NOW add species distribution modeling to calculate climatic optima
library('dismo')
library('raster')
library('rJava')

aet <- raster('data/gis_data/CAaet.tiff')
cwd <- raster('data/gis_data/CAcwd.tiff')
evar <- stack(aet,cwd)
head(cd)

mx <- maxent(evar,cd[,2:3])

system.file('java',package='dismo')
