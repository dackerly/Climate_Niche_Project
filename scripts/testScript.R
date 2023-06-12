# test script for testing functions
rm(list=ls())
library('mvtnorm')
source('scripts/climNiche.R')

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
plot(aet,xlim=rangeExtend(cd$longitude),ylim=rangeExtend(cd$latitude))
points(cd$long,cd$lat)
plot(cd$cwd,cd$aet)

## see climNiche.R for specification of the input format
cn <- climNicheDirect(cd)

cn$climStats
cn$covMatrix
cn$pcResults
cn$pcResults$loadings
head(cn$climData)
