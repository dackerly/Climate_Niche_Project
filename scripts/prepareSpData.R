library('terra')

# set values for function testing; id and rst have to be commented out when sourcing from scratch since cch and rastS don't exist yet
#id <- cch
idNameVar <- 'current_name_binomial'
latCol <- 'latitude'
longCol <- 'longitude'
#rst=rastS # raster stack of environmental data
npa <- 10000 # number of pseudoabsences (0 to not create)
psa.source <- c('ras') # 'id' means pseudoabsences are provided in input data; 'ras' means to select them from raster
gridEnv <- T # create gridded env data
# next four variables can be set separately for each env layer
minFA <- c('set','set') # set=specify lower limit; adif=subtract amount minGE from min; rdif=divide min by minGE
minGE <- c(0,0)
maxFA <- c('rdif','rdif') # set=specify lower limit; adif=add amount maxGE to max; rdif=multiply max by maxGE
maxGE <- c(2,2)
verbose <- T

prepareSpData <- function(id=NULL,spName=NULL,rst=NULL,npa=10000,psa.source=c('id'),gridEnv=T,minFA=rep('set',2),minGE=rep(0,2),maxFA=rep('rdif',2),maxGE=rep(2,2),verbose=F) {
  
  ## input data.frame should have 3 columns
  # c1: species name
  # c3-4: longitude/latitude
  # can be large data set with other species if these are being used as the background data; or just the selected species if background data will be sampled from raster
  
  id.all <- id
  names(id.all) <- c('name','longitude','latitude')
  id <- id.all[id.all$name==spName,]
  
  # add presence/absence column, with presences
  id$pa <- 1
  id$name <- as.character(id$name)
  id$longitude <- as.numeric(id$longitude)
  id$latitude <- as.numeric(id$latitude)
  
  # warning if more than one species; but could be intentional by user, e.g. if analyzing a genus
  if (length(table(id$name))>1) print('WARNING: more than one species in input data')
  
  # rearrange columns
  id <- id[,c('name','pa','longitude','latitude')]
  
  # lots of verbose printouts for troubleshooting
  if (verbose) {
    print(dim(id))
    print(head(id))
    print(tail(id))
  }
  
  # add pseudoabsence values, sampled in this case from CCH to reflect sampling bias
  if (npa>0) {
    if (psa.source=='id') {
      rsamp <- sample(1:nrow(id.all),npa)
      psa <- data.frame(name=rep('PsAb',npa),pa=0,longitude=id.all[rsamp,'longitude'],latitude=id.all[rsamp,'latitude'])
      if (verbose) print(head(psa))
    } else {
      ll <- xyFromCell(rst,cells(rst))
      rsamp <- sample(1:nrow(ll),npa)
      psa <- data.frame(name=rep('PsAb',npa),pa=0,longitude=ll[rsamp,1],latitude=ll[rsamp,2])
    }
    id <- rbind(id,psa)
  }
  if (verbose) {
    print(head(id))
    print(tail(id))
  }
  # sample env data across all data points
  cvals <- data.frame(row=1:nrow(id))
  for (i in 1:nlyr(rst)) {
    tmp <- terra::extract(rst[[i]],id[,c('longitude','latitude')])[,2]
    cvals <- cbind(cvals,tmp)
    names(cvals)[ncol(cvals)] <- names(rst[[i]])
  }
  id <- cbind(id,cvals[,-1])
  vCols <- 5:ncol(id)
  
  if (verbose) {
    print(head(id));print(tail(id))
    
    # look at presence data
    op <- par(mfrow=c(1,2))
    plot(rst[[1]])
    points(id[id$pa==1,c('longitude','latitude')])
    plot(id[id$pa==1,vCols[1:2]])
    par(op)
  }
  
  
  # create gridded env data, if desired to extrapolate optimal conditions beyond environmental background
  minV <- rep(NA,nlyr(rst))
  maxV <- rep(NA,nlyr(rst))
  
  i=1
  for (i in 1:nlyr(rst)) {
    if (minFA[i]=='set') minV[i] <- minGE[i] else 
      if (minFA[i]=='adif') minV[i] <- min(id[,vCols[i]],na.rm=T) - minGE[i] else 
        if (minFA[i]=='rdif') minV[i] <- min(id[,vCols[i]],na.rm=T)/minGE[i]
        
        if (maxFA[i]=='set') maxV[i] <- maxGE[i] else 
          if (maxFA[i]=='adif') maxV[i] <- max(id[,vCols[i]],na.rm=T) + maxGE[i] else 
            if (maxFA[i]=='rdif') maxV[i] <- max(id[,vCols[i]],na.rm=T)*maxGE[i]
  }
  
  gVals <- list()
  i=1
  for (i in 1:nlyr(rst)) gVals[[i]] <- seq(minV[i],maxV[i],length.out=100)
  xx <- expand.grid(gVals)
  nd <- data.frame(name=rep('GriddedEnv',nrow(xx)),pa=c(-1),longitude=NA,latitude=NA)
  nd <- cbind(nd,xx)
  names(nd) <- names(id)
  id <- rbind(id,nd)
  
  if (verbose) head(nd)
  if (verbose) tail(nd)
  if (verbose) table(id$name)
  if (verbose) tail(id)
  
  # add quadratic values
  for (i in 1:length(vCols)) {
    sval <- id[,vCols[i]]^2
    id <- cbind(id,sval)
    names(id)[ncol(id)] <- paste(names(id)[vCols[i]],2,sep='')
  }
  if (verbose) head(id)
  return(id)
}