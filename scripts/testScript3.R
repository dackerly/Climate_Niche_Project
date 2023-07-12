# test script for testing functions
rm(list=ls())
library('terra')
library('dismo')
source('scripts/climNiche3.R')

# read in climate data using terra functions, and create stack
aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
rastS <- c(aet,cwd)
names(rastS) <- c('aet','cwd')

# also create raster class stack for maxent
rasterS <- stack(rastS)

# read in CA plant biodiversity database
cch <- readRDS('big_data/CCH_clean_data/California_Species_clean_All_epsg_3310.Rdata')
head(cch)

# prepare input data for one species
spname <- 'Quercus kelloggii'
spname <- 'Quercus douglasii'
cd <- cch[which(cch$current_name_binomial==spname),c('current_name_binomial','longitude','latitude')]
if (!'pa' %in% names(cd)) cd$pa <- 1
names(cd)[1] <- 'name'
cd$name <- as.character(cd$name)
cd <- cd[,c('name','pa','longitude','latitude')]
dim(cd)
head(cd)
tail(cd)

# for black oak there's one weird climate outlier - up near Shasta and in a very low cwd spot, throws off convex hull
if (spname=='Quercus kelloggii') cd <- cd[-688,]

## input data.frame should have minimum of 4 columns:
# c1: species name
# c2: presence/absence - all 1s for original obs data
# c3-4: longitude/latitude

# add pseudoabsence values, sampled in this case from CCH to reflectsampling bias
npa <- 10000
rsamp <- sample(1:nrow(cch),npa)

psa <- data.frame(name=rep('PsAb',npa),pa=0,longitude=cch$longitude[rsamp],latitude=cch$latitude[rsamp])
head(psa)

cd <- rbind(cd,psa)

# sample env data across all data points
cd$aet <- extract(aet,cd[,c('longitude','latitude')])[,2]
cd$cwd <- extract(cwd,cd[,c('longitude','latitude')])[,2]
head(cd)
tail(cd)

# look at presence data
op <- par(mfrow=c(1,2))
plot(aet)
points(cd[cd$pa==1,c('longitude','latitude')])
plot(cwd~aet,data=cd[cd$pa==1,])
par(op)


# create gridded env data, if desired to extrapolate optimal conditions beyond environmental background
aetVals <- seq(0,2*max(cd$aet,na.rm=T),length.out=100)
cwdVals <- seq(0,2*max(cd$cwd,na.rm=T),length.out=100)
xx <- expand.grid(aetVals,cwdVals)
nd <- data.frame(name='GriddedEnv',pa=c(-1),longitude=NA,latitude=NA,aet=xx[,1],cwd=xx[,2])
head(nd)
tail(nd)

cd <- rbind(cd,nd)
table(cd$name)
tail(cd)

# add quadratic values
cd$aet2 <- cd$aet^2
cd$cwd2 <- cd$cwd^2

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
cn <- climNiche3(cd,vCols=5:6,trunc=0.01,model=mx.fit)
cn$climStats

head(cn$climData)
tail(cn$climData)
hist(cn$climData$pVal[cn$climData$pa==1])

# plot multivariate normal probability against fitted model predicted value
plot(pVal~mnScaledProb,data=cn$climData[cn$climData$pa==1,])

# read in a color ramp
ramp <- read.csv('data/colorRamps/jja6ColRamp.csv',as.is=T)
nrow(ramp)
head(ramp)
plot(1:nrow(ramp),1:nrow(ramp),pch=19,col=ramp$Hex)

## now plot again with colors proportional to normalized multivariate probability
pVals <- as.numeric(cut(cn$climData$pVal[cn$climData$pa==1],nrow(ramp)))
plot(cwd)
points(cn$climData[cn$climData$pa==1,c('longitude','latitude')],pch=19,col=ramp$Hex[pVals])

pVals <- as.numeric(cut(cn$climData$mnScaledProb[cn$climData$pa==1],nrow(ramp)))
plot(cwd)
points(cn$climData[cn$climData$pa==1,c('longitude','latitude')],pch=19,col=ramp$Hex[pVals])

# show suitability values for all presence and absence points
pVals <- as.numeric(cut(cn$climData$pVal[cn$climData$pa>=0],nrow(ramp)))
plot(cwd)
points(cn$climData[cn$climData$pa>=0,c('longitude','latitude')],pch=19,cex=0.25,col=ramp$Hex[pVals])



