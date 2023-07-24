rm(list=ls())
library('terra')
library('dismo')
source('scripts/climNiche3.R')
source('scripts/prepareSpData.R')

# read in climate data using terra functions, and create stack
aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
tmn <- rast('data/gis_data/CAtmn.tiff')
rastS <- c(aet,cwd)
names(rastS) <- c('aet','cwd')
nlr <- nlyr(rastS)

# also create raster class stack for maxent
rasterS <- stack(rastS)

# read in CA plant biodiversity database
cch <- readRDS('big_data/CCH_clean_data/California_Species_clean_All_epsg_3310.Rdata')
cch$current_name_binomial <- as.character(cch$current_name_binomial)
names(cch)

cval.all <- data.frame(values(rastS))
cval <- cval.all[which(complete.cases(cval.all)),]
dim(cval)
rsamp <- sample(1:nrow(cval),20000)
cvch <- terra::convHull(cval)

cch$aet <- extract(aet,cch[,c('longitude','latitude')])[,2]
cch$cwd <- extract(cwd,cch[,c('longitude','latitude')])[,2]
vCols <- c(which(names(cch)=='aet'),which(names(cch)=='cwd'))

sp <- 'Quercus douglasii'
spr <- which(cch$current_name_binomial==sp)
#spr <- spr[-153]
spch <- terra::convHull(cch[spr,c('aet','cwd')])
head(cch[spr,])
str(spch)
length(spr)

plot(cval[rsamp,],pch=19,cex=0.2)
plot(cvch,add=T)
points(cch$aet[spr],cch$cwd[spr],pch=19,cex=0.4,col='blue')
plot(spch,add=T,lwd=2)

sv <- vect(as.matrix(cch[spr,]),type='points')
ch <- convHull(sv)

summary(cval[,2])
