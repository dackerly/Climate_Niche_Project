# test script for testing functions
rm(list=ls())
library('mvtnorm')
source('scripts/climNiche2.R')

aet <- raster('big_data/CHELSA_climate/CHELSA_annual_AET_1979_2013.tiff')
cwd <- raster('big_data/CHELSA_climate/CHELSA_annual_CWD_1979_2013.tiff')

rangeExtend <- function(x,extAmount=0.5) {
  rx <- range(x)
  rx[1] <- rx[1]-extAmount
  rx[2] <- rx[2]+extAmount
  return(rx)
}

cd <- read.csv('data/test_data/test_climDat.csv',row.names = 1)
head(cd)
op <- par(mfrow=c(1,2))
plot(aet,xlim=rangeExtend(cd$longitude),ylim=rangeExtend(cd$latitude))
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
st <- stack(cwd,aet)
st <- crop(st,extent(c(rangeExtend(cd$longitude),rangeExtend(cd$latitude))))
names(st) <- c('cwd','aet')

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
