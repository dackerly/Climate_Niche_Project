# niche simulation to test methods
rm(list=ls())
library('terra')
library('dismo')
source('scripts/climNiche3.R')
source('scripts/prepareSpData.R')
source('scripts/simDist.R')

# read in climate data using terra functions, and create stack
aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
rastS <- c(aet,cwd)
names(rastS) <- c('aet','cwd')
nlr <- nlyr(rastS)

cV <- values(rastS)
head(cV)
hasV <- which(complete.cases(cV))
length(hasV)
mnV <- apply(cV[hasV,],2,mean)
sdV <- apply(cV[hasV,],2,sd)

zTran <- function(x,mn,sd) (x-mn)/sd
cVz <- cV
for (i in 1:nlr) cVz[hasV,i] <- zTran(cV[hasV,i],mnV[i],sdV[i])

rVz <- data.frame(min=rep(NA,nlr),max=rep(NA,nlr))
for (i in 1:nlr) rVz[i,] <- range(cVz[hasV,i])
row.names(rVz) <- names(rastS)
rVz
rVar <- data.frame(min=c(0.2,0.5),max=c(0.5,0.8))
rCov <- c(0,0)

i=1
for (i in 1:100) {
  simD <- simDist(rast=rastS,mnR=rVz,varR=rVar,covR=rCov,maxPts=5000)
  simD$means
  plot(simD$ras)
  
  dVal <- values(simD$ras)
  #table(dVal,useNA = 'always')
  spname <- 'sim1'
  ll.vals <- xyFromCell(aet,1:length(values(aet)))
  ll.vals <- ll.vals[which(dVal==1),]
  cd <- data.frame(name=spname,longitude=ll.vals[,1],latitude=ll.vals[,2])
  #head(cd)
  #plot(cd[,2:3],asp=1)
  
  cd <- prepareSpData(id=cd,spName='sim1',rst=rastS,psa.source='ras')
  #head(cd)
  #tail(cd)
  #table(cd$name)
  
  # fit a glm to the presence data, against pseudoabsence
  glm.fit <- glm(pa~aet+cwd+aet2+cwd2+aet:cwd,data=cd[cd$pa %in% c(0,1),],family='binomial')
  
  # fit maxent
  mx.fit <- dismo::maxent(cd[cd$pa %in% c(0,1),5:6],cd$pa[cd$pa %in% c(0,1)])
  
  cn <- climNiche3(cd,vCols=5:6,model=mx.fit)
  cn$climStats$simMns <- simD$means
  #print(cn$climStats)

  if (i==1) {
    aetVals <- cn$climStats[1,]
    cwdVals <- cn$climStats[2,]
  } else {
    aetVals <- rbind(aetVals,cn$climStats[1,])
    cwdVals <- rbind(cwdVals,cn$climStats[2,])
  }
}
aetVals
plot(simMns~mpt,data=aetVals);abline(0,1)
plot(simMns~mat,data=aetVals);abline(0,1)
plot(simMns~opt,data=aetVals);abline(0,1)
plot(simMns~pmn,data=aetVals);abline(0,1)
plot(mpt~pmn,data=aetVals);abline(0,1)
