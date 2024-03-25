# script to calculate climate niche stats for some or all species in CCH database
rm(list=ls())
library('terra')
library('dismo')
source('scripts/climNiche3.R')
source('scripts/prepareSpData.R')

# read in climate data using terra functions, and create stack
aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
tmn <- rast('data/gis_data/CAtmn.tiff')
ppt <- rast('data/gis_data/CAppt.tiff')
tav <- rast('data/gis_data/CAtav.tiff')

rastS <- c(ppt,tmn)
names(rastS) <- c('ppt','tmn')
nlr <- nlyr(rastS)

# also create raster class stack for maxent
rasterS <- stack(rastS)

# Identify data set
dset <- 'CCH_clean_2016'
#dset <- 'TP_non_native'
#dset <- 'PWD_non_native'

if (dset=='CCH_clean_2016') {
  # read in CA plant biodiversity database
  cch <- readRDS('big_data/CCH_clean_data/California_Species_clean_All_epsg_3310.Rdata')
  cch$current_name_binomial <- as.character(cch$current_name_binomial)
  names(cch)
  
  allSpecies <- sort(unique(cch$current_name_binomial))
  length(allSpecies)
  
  #(allSpecies <- allSpecies[1:10])
  #quercus <- grep('Quercus',allSpecies)
  #allSpecies <- allSpecies[quercus]
}

if (dset=='TP_non_native') {
  cch <- read.csv('big_data/CCH_non_native/TP_non_native.csv',as.is=T)
  names(cch)[grep('eFlora',names(cch))] <- 'current_name_binomial'
  names(cch)[grep('decimalLatitude',names(cch))] <- 'latitude'
  names(cch)[grep('decimalLongitude',names(cch))] <- 'longitude'
  
  allSpecies <- sort(unique(cch$current_name_binomial))
  length(allSpecies)
  
  #(allSpecies <- allSpecies[1:1])
}

if (dset=='PWD_non_native') {
  cch <- read.csv('big_data/CCH_non_native/PWD_non_native.csv',as.is=T)
  names(cch)[grep('ACCName',names(cch))] <- 'current_name_binomial'
  names(cch)[grep('decimalLatitude',names(cch))] <- 'latitude'
  names(cch)[grep('decimalLongitude',names(cch))] <- 'longitude'
  
  allSpecies <- sort(unique(cch$current_name_binomial))
  length(allSpecies)
  
  #(allSpecies <- allSpecies[1:1])
}

if (TRUE) { # CHANGE TO TRUE TO RERUN, OTHERWISE READ IN RESULTS BELOW
  i=1
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
        cnres <- data.frame(name=selSp,cn$climStats)
        #aetAll <- data.frame(name=selSp,cn$climStats[1,])
        #cwdAll <- data.frame(name=selSp,cn$climStats[2,])
        #tmnAll <- data.frame(name=selSp,cn$climStats[3,])
      } else {
        cnres <- rbind(cnres,data.frame(name=selSp,cn$climStats))
        #aetAll <- rbind(aetAll,data.frame(name=selSp,cn$climStats[1,]))
        #cwdAll <- rbind(cwdAll,data.frame(name=selSp,cn$climStats[2,]))
        #tmnAll <- rbind(tmnAll,data.frame(name=selSp,cn$climStats[3,]))
      }
    }
  }
  write.csv(cnres,paste('results/',dset,'_results.csv',sep=''))
  # head(aetAll)
  # head(cwdAll)
  # #head(tmnAll)
  # dim(aetAll)
  # write.csv(aetAll,'results/cchAET.csv')
  # write.csv(cwdAll,'results/cchCWD.csv')
  
}

#aetAll <- read.csv('results/cchAET.csv')
#cwdAll <- read.csv('results/cchCWD.csv')
#names(aetAll)
#head(aetAll)
#dim(aetAll)

#pairs(aetAll[,c('pmn','wtm','chc','mpt','mat','opt')])
#plot(pmn~mpt,data=cwdAll,xlim=c(200,1800),ylim=c(200,1800))
# abline(0,1,lwd=2,col='red')
# 
# par(mfrow=c(1,1))
# plot(pmn~mpt,data=aetAll,pch=19,col='darkgray',xlab='AET at maximum suitability',ylab='Niche mean AET across range',cex.lab=1.5)
# abline(0,1,lwd=2)
# 
# names(aetAll)
# mvar <- 14
# plot(aetAll[,mvar],cwdAll[,mvar])
# #plot(aetAll[,mvar],tmnAll[,mvar])
# 
# 
