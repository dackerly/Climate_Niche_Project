# test script for testing functions
rm(list=ls())
library('terra')
library('dismo')
source('scripts/climNiche3.R')
source('scripts/prepareSpData.R')

# read in climate data using terra functions, and create stack
aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
rastS <- c(aet,cwd)
names(rastS) <- c('aet','cwd')
nlr <- nlyr(rastS)

# also create raster class stack for maxent
rasterS <- stack(rastS)

# read in CA plant biodiversity database
cch <- readRDS('big_data/CCH_clean_data/California_Species_clean_All_epsg_3310.Rdata')
cch$current_name_binomial <- as.character(cch$current_name_binomial)
names(cch)

selSp <- 'Quercus douglasii'
inputData <- cch[,c('current_name_binomial','longitude','latitude')]

cd <- prepareSpData(id=inputData,spName=selSp,rst=rastS,npa=10000,psa.source='id',gridEnv=T,minFA=rep('rdif',nlr),minGE=rep(1,nlr),maxFA=rep('rdif',nlr),maxGE=rep(1,nlr),verbose=F)

dim(cd)
head(cd)
tail(cd)
table(cd$pa)
plot(cwd~aet,data=cd[cd$pa==0,])

# fit a glm to the presence data, against pseudoabsence
glm.fit <- glm(pa~aet+cwd+aet2+cwd2+aet:cwd,data=cd[cd$pa %in% c(0,1),],family='binomial')
glm.fit

# fit maxent
mx.fit <- dismo::maxent(cd[cd$pa %in% c(0,1),5:6],cd$pa[cd$pa %in% c(0,1)])
# examine maxent results will open page in web browser
#mx.fit
str(mx.fit)

# calculate climate niche statistics, using presence, absence, and (optionally GriddedEnv data), and the model
head(cd)
tail(cd)
dim(cd)
write.csv(cd,'data/test_data/expanded_test_data.csv')

### CALCULATE CLIMATE NICHE STATISTICS
source('scripts/climNiche3.R')
cn <- climNiche3(cd,trunc=0.01,vCols=5:6,model=glm.fit)
cn$climStats
cnd <- cn$climData

p <- which(cnd$pa==1)
a <- which(cnd$pa==0)
pa <- which(cnd$pa %in% c(0,1))
e <- which(cnd$pa==c(-1))

head(cnd)
tail(cnd)
hist(cnd$pVal[p])

# plot multivariate normal probability against fitted model predicted value
plot(pVal~mnScaledProb,data=cnd[p,])

# read in a color ramp
ramp <- read.csv('data/colorRamps/jja6ColRamp.csv',as.is=T)
nrow(ramp)
head(ramp)
plot(1:nrow(ramp),1:nrow(ramp),pch=19,col=ramp$Hex)

## now plot again with colors proportional to normalized multivariate probability
pVals <- as.numeric(cut(cnd$pVal[p],nrow(ramp)))
plot(cwd)
points(cnd[p,c('longitude','latitude')],pch=19,col=ramp$Hex[pVals])

pVals <- as.numeric(cut(cnd$mnScaledProb[p],nrow(ramp)))
plot(cwd)
points(cnd[p,c('longitude','latitude')],pch=19,col=ramp$Hex[pVals])

# show suitability values for all presence and absence points
pVals <- as.numeric(cut(cnd$pVal[pa],nrow(ramp)))
plot(cwd)
points(cnd[pa,c('longitude','latitude')],pch=19,cex=0.25,col=ramp$Hex[pVals])

pVals <- as.numeric(cut(cnd$mnScaledProb[pa],nrow(ramp)))
plot(cwd)
points(cnd[pa,c('longitude','latitude')],pch=19,cex=0.25,col=ramp$Hex[pVals])

# visualize model in full climate space
pVals <- as.numeric(cut(cnd$pVal[e],nrow(ramp)))
plot(cnd$aet[e],cnd$cwd[e],cex=0.25,type='n')
points(cnd$aet[a],cnd$cwd[a],pch=19,cex=0.25)
points(cnd$aet[p],cnd$cwd[p],pch=19,col='blue',cex=0.25)
points(cnd[e,c('aet','cwd')],pch=1,cex=0.25,col=ramp$Hex[pVals])

pVals <- as.numeric(cut(cnd$mnScaledProb[e],nrow(ramp)))
plot(cnd$aet[e],cnd$cwd[e],cex=0.25,type='n')
points(cnd$aet[pa],cnd$cwd[pa],pch=19,cex=0.25)
points(cnd$aet[p],cnd$cwd[p],pch=19,col='blue',cex=0.25)
points(cnd[e,c('aet','cwd')],pch=1,cex=0.25,col=ramp$Hex[pVals])


