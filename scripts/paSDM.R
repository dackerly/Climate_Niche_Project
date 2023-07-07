# quick function to run a quadratic logistic model on presence points and equal number of pseudoabsences from underlying rasters

# cd is a data.frame with:
# col 1: species
# col 2: longitude
# col 3: latitude
# col 4: env var 1
# col 5: env var 2
# col 6 (optional): pa - presence/absence
library('lme4')
pasdm <- function(cd,st,npa=10000) {
  vnames <- names(cd)[4:5]
  
  # if presence/absence points already present, skip directly to running model
  if (!'pa' %in% names(cd)) {
    cd$pa <- 1 # create pa and set all points as 'presence'
    
    # sample equal number of random pseudoabsence points
    stv <- values(st)
    dim(stv)
    stv <- stv[which(is.finite(stv[,1])),]
    dim(stv)
    head(stv)
    
    rsamp <- sample(1:nrow(stv),npa)
    abs <- stv[rsamp,] # absence points
    head(abs)
    
    # put absence points in new data.frame with same structure as cd
    cda <- data.frame(x1=rep(NA,npa),x2=NA,x3=NA,x4=NA,x5=NA,pa=0)
    names(cda) <- names(cd)
    cda[,4:5] <- abs
    head(cda)
    
    # combine into one data.frame with presence and absence points
    cd2 <- rbind(cd,cda)
    
    # squared terms for quadratics
    cd2$aet2 <- cd2$aet^2
    cd2$cwd2 <- cd2$cwd^2
  }
  # fit logistic
  fit <- glm(pa~aet+cwd+aet2+cwd2,data=cd2,family='binomial')
  
  # calculate predicted values
  cd2$pVal <- predict(fit,cd2)
  
  # find point within presence data with highest pVal
  mpt <- which.max(cd2$pVal[which(cd2$pa==1)])
  
  # create artificial data, extending env vars from 0 to 2X
  aetVals <- seq(0,2*max(cd$aet),length.out=100)
  cwdVals <- seq(0,2*max(cd$cwd),length.out=100)
  
  xx <- expand.grid(aetVals,cwdVals)
  nd <- data.frame(aet=xx[,1],cwd=xx[,2])
  nd$aet2 <- nd$aet^2
  nd$cwd2 <- nd$cwd^2
  
  nd$pVal <- predict(fit,nd)
  opt <- which.max(nd$pVal)
  return(list(cd2,nd,mpt,opt,npa))
}

