# test script for testing functions
rm(list=ls())
library('terra')
library('dismo')
library('raster')
source('scripts/climNiche3.R')
source('scripts/paSDM.R')

aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')

# prepare raster stack
st <- c(aet,cwd);names(st) <- c('aet','cwd')

rangeExtend <- function(x,extAmount=0.5) {
  rx <- range(x)
  rx[1] <- rx[1]-extAmount
  rx[2] <- rx[2]+extAmount
  return(rx)
}

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

# fit a glm to the presence data, against pseudoabsence
# add quadratic values
cd$aet2 <- cd$aet^2
cd$cwd2 <- cd$cwd^2

# fit logistic
fit <- glm(pa~aet+cwd+aet2+cwd2+aet:cwd,data=cd[cd$pa %in% c(0,1),],family='binomial')
fit

# calculate climate niche statistics, using presence, absence, and (optionally GriddedEnv data), and the model
head(cd)
tail(cd)
dim(cd)
write.csv(cd,'data/test_data/expanded_test_data.csv')

cn <- climNiche3(cd,vCols=5:6,trunc=0.01,model=fit)
cn$climStats

head(cn$climData)
tail(cn$climData)
hist(cn$climData$pVal[cn$climData$pa==1])

# read in a color ramp
ramp <- read.csv('data/colorRamps/jja6ColRamp.csv',as.is=T)
nrow(ramp)
head(ramp)
plot(1:nrow(ramp),1:nrow(ramp),pch=19,col=ramp$Hex)

## now plot again with colors proportional to normalized multivariate probability
pVals <- as.numeric(cut(cn$climData$pVal[cn$climData$pa==1],nrow(ramp)))
plot(cwd)
points(cd[cd$pa==1,c('longitude','latitude')],pch=19,col=ramp$Hex[pVals])

plot(cwd)
points(cd[cd$pa>=0,c('longitude','latitude')],pch=1,col=ramp$Hex[pVals])

## NOW add species distribution modeling to calculate climatic optima
head(cd)
source('scripts/paSDM.R')

glmRes <- pasdm(cd,st)
cd2 <- glmRes[[1]]
nd <- glmRes[[2]]
mpt <- glmRes[[3]]
mat <- glmRes[[4]]
opt <- glmRes[[5]]

# plot(cd2[,4:5])
# points(cd[,4:5],pch=1,col='red')
# plot(pa~pVal,data=cd2)

plot(cd2[,4:5],xlim=c(0,max(nd$aet)),ylim=c(0,max(nd$cwd)))
points(cd[,4:5],pch=19,col='red')
points(cd2[mpt,4:5],pch=19,col='blue',cex=2)
points(nd[opt,1:2],pch=19,col='green',cex=2)

cn$climStats$mpt <- as.numeric(cd2[mpt,4:5])
cn$climStats$mat <- as.numeric(cd2[mat,4:5])
cn$climStats$opt <- as.numeric(nd[opt,1:2])

cn$climStats

# now try maxent
head(cd)
dim(cd)

aetr <- raster('data/gis_data/CAaet.tiff')
cwdr <- raster('data/gis_data/CAcwd.tiff')
str <- stack(aetr,cwdr)

mx <- maxent(str,cd[,2:3])
mx
str(mx)
plot(mx)
head(mx@presence)
head(mx@absence)
