## now run stats for CA flora
rm(list=ls())
library('mvtnorm')
library('terra')
source('scripts/climNiche2.R')
source('scripts/paSDM.R')

rangeExtend <- function(x,extAmount=0.5) {
  rx <- range(x)
  rx[1] <- rx[1]-extAmount
  rx[2] <- rx[2]+extAmount
  return(rx)
}

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
(CAplus <- ext(rangeExtend(cch$longitude),rangeExtend(cch$latitude)))

# read in climate data, already cropped to CA shape
aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
st <- c(aet,cwd)
plot(st[[1]])
points(cch[rsamp,c('longitude','latitude')],cex=0.2)

head(cch)
cch$aet <- extract(aet,cch[,2:3])[,2]
cch$cwd <- extract(cwd,cch[,2:3])[,2]
head(cch)

# pull out unique names
uNames <- sort(as.character(unique(cch$name)))
head(uNames)
length(uNames)

# analyze first species to set up results data.frames
cd <- cch[cch$name==uNames[1],]
dim(cd)
cn <- climNiche(cd,st,onePerCell=T,trunc=0.05)
cn$climStats

glmRes <- pasdm(cd,st)
cd2 <- glmRes[[1]]
nd <- glmRes[[2]]
mpt <- glmRes[[3]]
mat <- glmRes[[4]]
opt <- glmRes[[5]]

cn$climStats$mpt <- as.numeric(cd2[mpt,4:5])
cn$climStats$mat <- as.numeric(cd2[mat,4:5])
cn$climStats$opt <- as.numeric(nd[opt,1:2])
cn$climStats

aetAll <- data.frame(name=uNames[1],cn$climStats[1,])
cwdAll <- data.frame(name=uNames[1],cn$climStats[2,])

# analyze random sample, or all the rest - comment out as needed
rsamp <- 500
ssel <- sample(2:length(uNames),rsamp)
ssel <- 2:length(uNames)

for (i in ssel) {
  print(i)
  sname <- uNames[i]
  cd <- cch[cch$name==sname,]
  
  # filter to minimum number of observations
  if (nrow(cd)>=20) {
    cn <- climNiche(cd,st,onePerCell=T,trunc=0.05)
    glmRes <- pasdm(cd,st)
    cd2 <- glmRes[[1]]
    nd <- glmRes[[2]]
    cn$climStats$mpt <- as.numeric(cd2[glmRes[[3]],4:5])
    cn$climStats$mat <- as.numeric(cd2[glmRes[[4]],4:5])
    cn$climStats$opt <- as.numeric(nd[glmRes[[5]],1:2])
    
    aetAll <- rbind(aetAll,data.frame(name=sname,cn$climStats[1,]))
    cwdAll <- rbind(cwdAll,data.frame(name=sname,cn$climStats[2,]))
  }
}
dim(aetAll)
head(aetAll)
head(cwdAll)

# examine results
plot(aetAll$mean,aetAll$wtd.mean)
plot(aetAll$wtd.mean,aetAll$mpt)
plot(aetAll$wtd.mean,aetAll$mat)
plot(aetAll$mean,aetAll$opt)

abline(0,1)
plot(cwdAll$mean,cwdAll$wtd.mean)
plot(cwdAll$mean,cwdAll$mat)
plot(cwdAll$mean,cwdAll$opt)

abline(0,1)
plot(aetAll$mean,cwdAll$mean)

aetVals <- values(st[[1]])
cwdVals <- values(st[[2]])
length(aetVals)
rsamp <- sample(1:length(aetVals),100000)
plot(aetVals[rsamp],cwdVals[rsamp],col='gray',pch=19,cex=0.5)
points(aetAll$mean,cwdAll$mean,pch=19,col='red')

plot(aetVals[rsamp],cwdVals[rsamp],col='gray',pch=19,cex=0.5)
points(aetAll$wtd.mean,cwdAll$wtd.mean,pch=19,col='blue')

plot(aetVals[rsamp],cwdVals[rsamp],col='gray',pch=19,cex=0.5)
points(aetAll$wtd.mean,cwdAll$mpt,pch=19,col='blue')

plot(aetVals[rsamp],cwdVals[rsamp],col='gray',pch=19,cex=0.5)
points(aetAll$wtd.mean,cwdAll$mat,pch=19,col='blue')

plot(aetVals[rsamp],cwdVals[rsamp],col='gray',pch=19,cex=0.5,xlim=c(0,600),ylim=c(0,4000))
points(aetAll$wtd.mean,cwdAll$opt,pch=19,col='blue')

