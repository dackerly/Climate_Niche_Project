library('terra')

# set values for function testing
id <- cch
idNameVar <- 'current_name_binomial'
latCol <- 'latitude'
longCol <- 'longitude'
spName <- 'Quercus douglasii'
rst=rastS
npa <- 10000
psa.source <- c('ras')
gridEnv <- T
minFA <- 'A'
minGE <- 0
maxFA <- 'F'
maxGE <- 2
verbose <- T

prepareSpData <- function(id=cch,idNameVar='name',latCol='latitude',longCol='longitude',spName='Quercus douglasii',rst=rStack,npa=10000,psa.source=c('id'),gridEnv=T,minFA='A',minGE=0,maxFA='F',maxGE=2,verbose=F) {
  
  cd <- id[which(id[,idNameVar]==spName),c(idNameVar,longCol,latCol)]
  names(cd) <- c('name','longitude','latitude')
  if (!'pa' %in% names(cd)) cd$pa <- 1
  cd$name <- as.character(cd$name)
  cd$longitude <- as.numeric(cd$longitude)
  cd$latitude <- as.numeric(cd$latitude)
  cd <- cd[,c('name','pa','longitude','latitude')]
  if (verbose) {
    print(dim(cd))
    print(head(cd))
    print(tail(cd))
  }
  
  ## input data.frame should have minimum of 4 columns:
  # c1: species name
  # c2: presence/absence - all 1s for original obs data
  # c3-4: longitude/latitude
  
  # add pseudoabsence values, sampled in this case from CCH to reflectsampling bias
  if (npa>0) {
    if (psa.source=='id') {
      rsamp <- sample(1:nrow(cch),npa)
      psa <- data.frame(name=rep('PsAb',npa),pa=0,longitude=id[rsamp,longCol],latitude=id[rsamp,latCol])
      if (verbose) print(head(psa))
    } else {
      ll <- xyFromCell(rst,cells(rst))
      rsamp <- sample(1:nrow(ll),npa)
      psa <- data.frame(name=rep('PsAb',npa),pa=0,longitude=ll[rsamp,1],latitude=ll[rsamp,2])
    }
    cd <- rbind(cd,psa)
  }
  if (verbose) {
    print(head(cd))
    print(tail(cd))
  }
  # sample env data across all data points
  cvals <- data.frame(row=1:nrow(cd))
  for (i in 1:nlyr(rst)) {
    tmp <- terra::extract(rst[[i]],cd[,c('longitude','latitude')])[,2]
    cvals <- cbind(cvals,tmp)
    names(cvals)[ncol(cvals)] <- names(rst[[i]])
  }
  cd <- cbind(cd,cvals[,-1])
  
  if (verbose) {
    print(head(cd));print(tail(cd))
    
    # look at presence data
    op <- par(mfrow=c(1,2))
    plot(aet)
    points(cd[cd$pa==1,c('longitude','latitude')])
    plot(cwd~aet,data=cd[cd$pa==1,])
    par(op)
  }
  
  ## LUNCH HERE
  
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
  return(cd)
}