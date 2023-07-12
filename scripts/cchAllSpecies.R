# script to calculate climate niche stats for some or all species in CCH database
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

allSpecies <- sort(unique(cch$current_name_binomial))
#quercus <- grep('Quercus',allSpecies)
#allSpecies <- allSpecies[quercus]

i=231
length(allSpecies)

for (i in 1:length(allSpecies)) {
  print(i)
  selSp <- allSpecies[i]
  inputData <- cch[,c('current_name_binomial','longitude','latitude')]
  
  cd <- prepareSpData(id=inputData,spName=selSp,rst=rastS,npa=10000,psa.source='id',gridEnv=T,minFA=rep('set',nlr),minGE=rep(0,nlr),maxFA=rep('rdif',nlr),maxGE=rep(2,nlr),verbose=F)
  
  if (length(which(cd$pa==1))>=20) {
    # fit a glm to the presence data, against pseudoabsence
    # glm.fit <- glm(pa~aet+cwd+aet2+cwd2+aet:cwd,data=cd[cd$pa %in% c(0,1),],family='binomial')
    # glm.fit
    
    # fit maxent
    mx.fit <- dismo::maxent(cd[cd$pa %in% c(0,1),5:6],cd$pa[cd$pa %in% c(0,1)])
    # examine maxent results will open page in web browser
    #mx.fit
    
    ### CALCULATE CLIMATE NICHE STATISTICS
    cn <- climNiche3(cd,vCols=5:6,trunc=0.01,model=mx.fit)
    
    if (i==1) {
      aetAll <- data.frame(name=selSp,cn$climStats[1,])
      cwdAll <- data.frame(name=selSp,cn$climStats[2,])
    } else {
      aetAll <- rbind(aetAll,data.frame(name=selSp,cn$climStats[1,]))
      cwdAll <- rbind(cwdAll,data.frame(name=selSp,cn$climStats[2,]))
    }
  }
}
head(aetAll)
head(cwdAll)

pairs(aetAll[,c('mean','wtd.mean','chc','pwt.mean','mpt','mat','opt')])

names(aetAll)
mvar <- 14
plot(aetAll[,mvar],cwdAll[,mvar])

