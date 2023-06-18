## now run stats for CA flora
rm(list=ls())
library('raster')
library('sp')
source('scripts/climNiche2.R')

rangeExtend <- function(x,extAmount=0.5) {
  rx <- range(x)
  rx[1] <- rx[1]-extAmount
  rx[2] <- rx[2]+extAmount
  return(rx)
}

# read in CA shapefile....
#CA <- sp::

# read in CCH data set (Baldwin et al.), and extract just binomial scientific name, longitude, and latitude
cch <- readRDS('big_data/CCH_clean_data/California_Species_clean_All_epsg_3310.Rdata')
head(cch)
names(cch)
cch <- cch[,c('current_name_binomial','longitude','latitude')]
names(cch)[1] <- 'name'
head(cch)
dim(cch)

# subsample and plot locations
rsamp <- sample(1:nrow(cch),50000)
plot(cch[rsamp,c('longitude','latitude')],asp=1)
(CAplus <- extent(rangeExtend(cch$longitude),rangeExtend(cch$latitude)))

# read in climate data and crop to CA long,lat bounding box +/- 0.5 degrees latitude
aet <- raster('big_data/CHELSA_climate/CHELSA_annual_AET_1979_2013.tiff')
cwd <- raster('big_data/CHELSA_climate/CHELSA_annual_CWD_1979_2013.tiff')
st <- crop(stack(aet,cwd),CAplus)
plot(st[[1]])
points(cch[rsamp,c('longitude','latitude')],cex=0.2)

# pull out unique names
uNames <- sort(as.character(unique(cch$name)))
head(uNames)
length(uNames)

# analyze first species to set up results data.frames
cd <- cch[cch$name==uNames[1],]
dim(cd)
cn <- climNiche(cd,st,onePerCell=T,trunc=0.05)
cn$climStats
aetAll <- data.frame(name=uNames[1],cn$climStats[1,])
cwdAll <- data.frame(name=uNames[1],cn$climStats[2,])

# analyze all the rest
i=231
for (i in 2:length(uNames)) {
  print(i)
  sname <- uNames[i]
  cd <- cch[cch$name==sname,]
  if (nrow(cd)>=20) {
    cn <- climNiche(cd,st,onePerCell=T,trunc=0.05)
    aetAll <- rbind(aetAll,data.frame(name=sname,cn$climStats[1,]))
    cwdAll <- rbind(cwdAll,data.frame(name=sname,cn$climStats[2,]))
  }
}
head(aetAll)
head(cwdAll)

# examine results
plot(aetAll$mean,aetAll$wtd.mean)
abline(0,1)
plot(cwdAll$mean,cwdAll$wtd.mean)
abline(0,1)

aetVals <- getValues(st[[1]])
cwdVals <- getValues(st[[2]])
length(aetVals)
rsamp <- sample(1:length(aetVals),100000)
plot(cwdVals[rsamp],aetVals[rsamp],col='gray',pch=19,cex=0.5)
points(cwdAll$mean,aetAll$mean,pch=19,col='red')
points(cwdAll$wtd.mean,aetAll$wtd.mean,pch=19,col='blue')

