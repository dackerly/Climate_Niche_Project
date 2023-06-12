# test script for testing functions
rm(list=ls())
library('raster')
source('scripts/climNiche.R')

# use this to extract data for different species
if (FALSE) {
  cch <- readRDS('big_data/CCH_clean_data/California_Species_clean_All_epsg_3310.Rdata')
  bk <- cch[which(cch$scientificName=='Quercus kelloggii'),]
  dim(bk)
  write.csv(bk,'data/test_data/black_oak.csv',row.names = F)
}

bk <- read.csv('data/test_data/black_oak.csv')
dim(bk)
names(bk)
plot(bk$longitude,bk$latitude,asp=1)

aet <- raster('big_data/CHELSA_climate/CHELSA_annual_AET_1979_2013.tiff')
cwd <- raster('big_data/CHELSA_climate/CHELSA_annual_CWD_1979_2013.tiff')
plot(aet)

bk$aet <- extract(aet,bk[,c('longitude','latitude')])
hist(bk$aet)
bk$cwd <- extract(cwd,bk[,c('longitude','latitude')])
hist(bk$cwd)
plot(bk$cwd,bk$aet)

names(bk)
cd <- bk[,c('scientificName','longitude','latitude','aet','cwd')]
cd <- cd[is.finite(cd$latitude) & is.finite(cd$longitude),]

write.csv(cd,'data/test_data/test_climDat.csv')
