# climate niche functions
library('mvtnorm')
library('mnormt')
require('terra')

# temporary for code development
#cd <- read.csv('data/test_data/expanded_test_data.csv',row.names = 1)

# cd is a point data set of climate data for a species
# col1: species name
# col2: presence/absence (and -1 for GriddedEnv)
# col3-5: latitude,longitude (in any units)
# col5+: optional - climate data (any number of variables)
#   variable names must match variable in fit model object if predicted values desired
climNiche3 <- function(cd,vCols=NULL,trunc=0,comp.cases=T,model=NULL) {
  # set up empty list to capture results
  cn <- list()
  
  # utility vars: names, number, and position of climate variables in cd dataframe
  # if vCols isn't specified, then clim variables are all columns after long,lat. Specifying vCols is an option for the user to restrict the columns used for analysis, or to exclude quadratic variables that may be in the model fit object
  if (is.null(vCols)) vCols <- 5:ncol(cd)
  climVarNames <- names(cd)[vCols]
  nClim <- length(climVarNames)
  
  #remove records with incomplete climate data in any variable; only multivariate complete records are used for all calculations. Rownames created so that returned cd object shows which rows were included
  row.names(cd) <- 1:nrow(cd)
  allFinite <- function(x) all(is.finite(x))
  cd <- cd[apply(cd[,vCols],1,allFinite),]
  
  # set up data frame for descriptive stats, as first object in cn results list
  cn[[1]] <- data.frame(climVar=climVarNames,mean=NA,sd=NA,q05=NA,median=NA,q95=NA,min=NA,max=NA)
  names(cn) <- 'climStats'
  
  # p is rows for presence data
  # a is absence of pseudoabsence
  # e is optional, other environmental setting to project into
  p <- which(cd$pa==1)
  a <- which(cd$pa==0)
  pa <- which(cd$pa %in% c(0,1))
  e <- which(cd$pa==c(-1))
  
  #step through climate variables and calculate descriptive stats
  i=2
  for (i in 1:nClim) {
    cc <- vCols[i]
    zCol <- paste(climVarNames[i],'.z',sep='')
    
    # calculate z-score; calculated for all rows, but only on mean and sd for presence data
    zs <- (cd[,cc]-mean(cd[p,cc]))/sd(cd[p,cc])
    cd <- cbind(cd,zs)
    names(cd)[length(names(cd))] <- zCol
    
    cn$climStats$mean[i] <- mean(cd[p,cc],na.rm=T)
    cn$climStats$sd[i] <- sd(cd[p,cc],na.rm=T)
    cn$climStats$q05[i] <- quantile(cd[p,cc],probs=0.05,na.rm=T)
    cn$climStats$median[i] <- quantile(cd[p,cc],probs=0.50,na.rm=T)
    cn$climStats$q95[i] <- quantile(cd[p,cc],probs=0.95,na.rm=T)
    cn$climStats$min[i] <- min(cd[p,cc],na.rm=T)
    cn$climStats$max[i] <- max(cd[p,cc],na.rm=T)
  }
  
  # Calculate covariance matrix and bivariate probability score of each site; this can only be done on complete cases
  zCols <- (ncol(cd)-(nClim-1)):ncol(cd)
  zCov <- cov(cd[p,zCols])
  
  # calculate probability density at each point in MV space. Then normalize by density at origin to get relative probability values, as a measure of multivariate normal 'suitability'
  tmp <- dmnorm(cd[,zCols],mean=rep(0,nClim),varcov=zCov)
  dm.max <- dmnorm(rep(0,nClim),mean=rep(0,nClim),varcov=zCov)
  cd$mnScaledProb <- tmp/dm.max
  
  cn[[2]] <- zCov
  names(cn)[2] <- 'covMatrix'
  
  # Calculate principal components and provide PC object with results; add PC scores to cd data.frame
  pc <- princomp(cd[p,zCols],cor = F)
  cn[[3]] <- pc
  names(cn)[3] <- 'pcResults'
  
  # would be nice to add here projection of all rows onto principal components calculated only on presence data
  {}
  
  # if user has supplied absence or pseudoabsence data, then calculate climate means based on proportion of available data space used
  # if trunc>0, then truncate the requested proportion of records from either end of the climate values distribution. This is provided because one or a few extreme values that may fall in very rare climate types can skew the mean value a lot. Compare results with trunc=0 and trunc=0.01 or 0.05 to see the effect
  
  if (length(a)>0) {
    cn$climStats$wtd.mean <- NA
    i=2
    if (FALSE) {   # histogram based wtd mean 
      for (i in 1:nClim) {
        hClim <- hist(cd[pa,vCols[i]],breaks=20,plot=F)
        vals <- sort(cd[p,vCols[i]])
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
    if (TRUE) {
      for (i in 1:nClim) {
        mn <- min(cd[pa,vCols[i]],na.rm=T)
        mx <- max(cd[pa,vCols[i]],na.rm=T)
        hDen <- density(cd[a,vCols[i]],from=mn,to=mx)
        vals <- sort(cd[p,vCols[i]])
        if (trunc>0) {
          if (floor(length(vals)*trunc/2) < 1) tmp <- 1 else tmp <- floor(length(vals)*trunc/2)
          lt <- 1:tmp
          ut <- length(vals)+1-lt
          vals <- vals[-c(lt,ut)]
        }
        oDen <- density(vals,from=mn,to=mx)
        finDat <- which(is.finite(oDen$y/hDen$y))
        cn$climStats$wtd.mean[i] <- weighted.mean(hDen$x[finDat],oDen$y[finDat]/hDen$y[finDat])
      }
    }
  }
  
  # convex hull centroid
  sv <- vect(as.matrix(cd[p,vCols]),type='points')
  #plot(sv)
  ch <- convHull(sv)
  #lines(ch)
  cc <- centroids(ch)
  cn$climStats$chc <- ext(cc)[c(1,3)]
  
  # model fit
  if (!is.null(model)) {
    cn$climStats$pwt.mean <- NA
    cd$pVal <- predict(model,cd,type='response')
    
    # calculate weighted mean of the observed climate values, weighted by the predicted value from the model
    for (i in 1:nClim) cn$climStats$pwt.mean[i] <- stats::weighted.mean(cd[p,vCols[i]],cd$pVal[p])
    
    # find the presence point with the highest predicted value, and assign the climate values to mpt. If there is more than one that's identical, choose the one that is closest to the bivariate mean.
    mpt <- which.max(cd$pVal[p])
    cn$climStats$mpt <- as.numeric(cd[p[mpt],vCols])
    mpt.eq <- which(cd$pVal[p]==cd$pVal[p[mpt]])
    if (length(mpt.eq)>1) mpt <- mpt.eq[which.max(cd$mnScaledProb[p[mpt.eq]])]
    cn$climStats$mpt <- as.numeric(cd[pa[mpt],vCols])
    
    # if user has supplied absences, find the point within the union of presence and absences with highest predicted value. Same as above if there are multiple points with equivalent maximum predicted values. 
    if (length(a)>0) {
      mat <- which.max(cd$pVal[pa])
      mat.eq <- which(cd$pVal[pa]==cd$pVal[pa[mat]])
      if (length(mat.eq)>1) mat <- mat.eq[which.max(cd$mnScaledProb[pa[mat.eq]])]
      cn$climStats$mat <- as.numeric(cd[pa[mat],vCols])
    }
    
    # if user has supplied a gridded environment (or other env values), find the point with highest predicted value. Same as above if there are multiple points with equivalent maximum predicted values.  
    if (length(e)>0) {
      opt <- which.max(cd$pVal[e])
      opt.eq <- which(cd$pVal[e]==cd$pVal[e[opt]])
      if (length(opt.eq)>1) mat <- opt.eq[which.max(cd$mnScaledProb[pa[opt.eq]])]
      cn$climStats$opt <- as.numeric(cd[e[opt],vCols])
    }
  }
  # add updated climate data to cn object
  cn[[4]] <- cd
  names(cn)[4] <- 'climData'
  return(cn)
}