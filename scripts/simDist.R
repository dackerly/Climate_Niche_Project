# mnR: data.frame with nlr rows, cols: min,max for means
# varR: data.frame with lnr rows, cols: min,max for variances
# covR: 2 vals, min,max for range of covariance
library(terra)
simDist <- function(rast=rastS,mnR,varR,covR,maxPts=20000) {
  cV <- values(rast)
  hasV <- which(complete.cases(cV))
  mnV <- apply(cV[hasV,],2,mean)
  sdV <- apply(cV[hasV,],2,sd)
  
  nlr <- nlyr(rastS)
  mns <- rep(NA,nlr)
  for (i in 1:nlr) mns[i] <- c(runif(1,rVz[i,1],rVz[i,2]))
  
  vcv <- matrix(NA,2,2)
  for (i in 1:nlr) vcv[i,i] <- c(runif(1,rVar[i,1],rVar[i,2]))
  vcv[1,2] <- runif(1,rCov[1],rCov[2])
  vcv[2,1] <- vcv[1,2]
  
  zTran <- function(x,mn,sd) (x-mn)/sd
  cVz <- cV
  for (i in 1:nlr) cVz[hasV,i] <- zTran(cV[hasV,i],mnV[i],sdV[i])
  
  mnprob <- dmnorm(cVz[hasV,],mean=mns,varcov=vcv)
  binV <- rbinom(length(mnprob),1,mnprob)
  presV <- which(binV==1)
  if (length(presV)>maxPts) {
    presV <- sample(presV,maxPts)
    binV <- rep(0,length(binV))
    binV[presV] <- 1
  }
  
  binRasV <- cV[,1]
  binRasV[hasV] <- binV
  distRas <- setValues(rast[[1]],binRasV)
  
  simMeans <- mnV+(mns*sdV)
  
  res <- list(simMeans,distRas)
  names(res) <- c('means','ras')
  return(res)
}
