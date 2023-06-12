# climate niche functions
library('mvtnorm')

# temporary for code development
cd <- read.csv('data/test_data/test_climDat.csv',row.names = 1)

# cd is a point data set of climate data for a species
# col1: species name
# col2-3: latitude,longitude (in any units)
# col4+ climate data (any number of variables)
climNicheDirect <- function(cd) {
  cn <- list()
  
  climVarNames <- names(cd)[4:length(names(cd))]
  nClim <- length(climVarNames)
  cCols <- 1:nClim+3
  
  #remove records with incomplete climate data
  row.names(cd) <- 1:nrow(cd)
  allFinite <- function(x) all(is.finite(x))
  cd <- cd[apply(cd[,cCols],1,allFinite),]
  
  # set up data frame for descriptive stats
  cn[[1]] <- data.frame(climVar=climVarNames,mean=NA,sd=NA,q05=NA,median=NA,q95=NA,min=NA,max=NA)
  names(cn) <- 'climStats'
  
  i=1
  for (i in 1:nClim) {
    cc <- cCols[i]
    zCol <- paste(climVarNames[i],'.z',sep='')
    zs <- (cd[,cc]-mean(cd[,cc]))/sd(cd[,cc])
    cd <- cbind(cd,zs)
    names(cd)[length(names(cd))] <- zCol
    
    cn$climStats$mean[i] <- mean(cd[,cc],na.rm=T)
    cn$climStats$sd[i] <- sd(cd[,cc],na.rm=T)
    cn$climStats$q05[i] <- quantile(cd[,cc],probs=0.05,na.rm=T)
    cn$climStats$median[i] <- quantile(cd[,cc],probs=0.50,na.rm=T)
    cn$climStats$q95[i] <- quantile(cd[,cc],probs=0.95,na.rm=T)
    cn$climStats$min[i] <- min(cd[,cc],na.rm=T)
    cn$climStats$max[i] <- max(cd[,cc],na.rm=T)
  }
  
  # Calculate bivariate normal distribution and bivariate z-score of each site
  zCols <- 1:nClim+(3+nClim)
  zCov <- cov(cd[,zCols])
  
  # now I want probability density at each point in MV space - this is not the right command
  tmp <- dmnorm(cd[,zCols],mean=rep(0,nClim),varcov=zCov)
  dm.max <- dmnorm(rep(0,nClim),mean=rep(0,nClim),varcov=zCov)
  cd$mnScaledProb <- tmp/dm.max
  
  cn[[2]] <- zCov
  names(cn)[2] <- 'covMatrix'
  
  # PCA
  pc <- princomp(cd[,zCols],cor = F)
  cd <- cbind(cd,pc$scores)
  
  cn[[3]] <- pc
  names(cn)[3] <- 'pcResults'

  # add updated climate data to cn object
  cn[[4]] <- cd
  names(cn)[4] <- 'climData'
  return(cn)
}