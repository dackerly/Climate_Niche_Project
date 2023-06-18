# climate niche functions
library('mvtnorm')

# temporary for code development
cd <- read.csv('data/test_data/test_climDat.csv',row.names = 1)

# cd is a point data set of climate data for a species
# col1: species name
# col2-3: latitude,longitude (in any units)
# col4+: optional - climate data (any number of variables); if raster stack is provided then climate data are extracted
# st: optional: raster stack with climate variable means; if not available then niche values are calculated on directly provided climVars
climNiche <- function(cd,st=NULL,onePerCell=F,trunc=0,comp.cases=T) {
  
  # if user supplied stack of climate rasters is available, then data.frame is trimmed to first three columns (species, long, lat) and climate values are extracted from provided raster stack. This ensures that names and order of climate variables in stack and in data.frame are aligned. If user does not provide raster stack, then climVars must be in cd dataframe and descriptive stats are calculated on data provided.  
  if (!is.null(st)) {
    cd <- cd[,1:3]
    for (i in 1:nlayers(st)) cd <- cbind(cd,extract(st[[i]],cd[,2:3]))
    names(cd)[-c(1:3)] <- names(st)
  }
  
  # set up empty list to capture results
  cn <- list()
  
  # utility vars: names, number, and position of climate variables in cd dataframe
  climVarNames <- names(cd)[4:length(names(cd))]
  nClim <- length(climVarNames)
  cCols <- 1:nClim+3
  
  #remove records with incomplete climate data in any variable; only multivariate complete records are used for all calculations. 
  row.names(cd) <- 1:nrow(cd)
  allFinite <- function(x) all(is.finite(x))
  cd <- cd[apply(cd[,cCols],1,allFinite),]
  
  # set up data frame for descriptive stats, as first object in cn results list
  cn[[1]] <- data.frame(climVar=climVarNames,mean=NA,sd=NA,q05=NA,median=NA,q95=NA,min=NA,max=NA)
  names(cn) <- 'climStats'
  
  #step through climate variables and calculate descriptive stats
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
  
  # Calculate covariance matrix and bivariate probability score of each site; this can only be done on complete cases
  zCols <- 1:nClim+(3+nClim)
  zCov <- cov(cd[,zCols])
  
  # calculate probability density at each point in MV space. Then normalize by density at origin to get relative probability values, as a measure of multivariate normal 'suitability'
  tmp <- dmnorm(cd[,zCols],mean=rep(0,nClim),varcov=zCov)
  dm.max <- dmnorm(rep(0,nClim),mean=rep(0,nClim),varcov=zCov)
  cd$mnScaledProb <- tmp/dm.max
  
  cn[[2]] <- zCov
  names(cn)[2] <- 'covMatrix'
  
  # Calculate principal components and provide PC object with results; add PC scores to cd data.frame
  pc <- princomp(cd[,zCols],cor = F)
  cd <- cbind(cd,pc$scores)
  
  cn[[3]] <- pc
  names(cn)[3] <- 'pcResults'
  
  # add updated climate data to cn object
  cn[[4]] <- cd
  names(cn)[4] <- 'climData'
  
  # if user has supplied climate data, then calculate climate means based on proportion of available data space used
  # if onePerCell, then trim records to one per raster cell
  # if trunc>0, then truncate the requested proportion of records from either end of the climate values distribution. This is provided because one or a few extreme values that may fall in very rare climate types can skew the mean value a lot. Compare results with trunc=0 and trunc=0.01 or 0.05 to see the effect
  if (!is.null(st)) {
    cd$cell <- cellFromXY(st[[1]],cd[,c('longitude','latitude')])
    if (onePerCell) cdu <- cd[match(unique(cd$cell),cd$cell),] else cdu <- cd
    
    cn$climStats$wtd.mean <- NA
    i=2
    for (i in 1:nClim) {
      hClim <- hist(getValues(st[[i]]),breaks=20,plot=F)
      
      vals <- sort(cdu[,cCols[i]])
      if (trunc>0) {
        if (floor(length(vals)*trunc/2) < 1) tmp <- 1 else tmp <- floor(length(vals)*trunc/2)
        lt <- 1:tmp
        ut <- length(vals)+1-lt
        vals <- vals[-c(lt,ut)]
      }
      
      oClim <- hist(vals,breaks=hClim$breaks,plot=F)
      cn$climStats$wtd.mean[i] <- weighted.mean(hClim$mids,oClim$counts/hClim$counts)
    }
  }
  return(cn)
}