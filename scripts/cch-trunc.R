rm(list=ls())
library('terra')
library('dismo')
source('scripts/climNiche3.R')
source('scripts/prepareSpData.R')

# read in climate data using terra functions, and create stack
aet <- rast('data/gis_data/CAaet.tiff')
cwd <- rast('data/gis_data/CAcwd.tiff')
tmn <- rast('data/gis_data/CAtmn.tiff')
rastS <- c(cwd,aet)
names(rastS) <- c('cwd','aet')
nlr <- nlyr(rastS)

# also create raster class stack for maxent
rasterS <- stack(rastS)
plot(rastS)

# read in CA plant biodiversity database
cch <- readRDS('big_data/CCH_clean_data/California_Species_clean_All_epsg_3310.Rdata')
cch$current_name_binomial <- as.character(cch$current_name_binomial)
names(cch)

cval.all <- data.frame(values(rastS))
cval <- cval.all[which(complete.cases(cval.all)),]
dim(cval)
rsamp <- sample(1:nrow(cval),50000)
cvch <- terra::convHull(cval)

# insert alpha hull and contour density

# extract climate for cch
cch$aet <- extract(aet,cch[,c('longitude','latitude')])[,2]
cch$cwd <- extract(cwd,cch[,c('longitude','latitude')])[,2]
vCols <- c(which(names(cch)=='aet'),which(names(cch)=='cwd'))

sp <- 'Sequoia sempervirens'
spr <- which(cch$current_name_binomial==sp)
#spr <- spr[-153]
spch <- terra::convHull(cch[spr,c('cwd','aet')])
head(cch[spr,])
str(spch)
length(spr)


# function returns T/F indicating whether each row of matrix x (or matrix y, if specified) 
# has a kernel density higher than probability p, with respect to the rows of x
kdp <- function(x, # 2-column matrix used to fit kernel density model
                y = NULL, # optional matrix for which density is to be predicted; if not supplied, predictions are made for x
                p = 0.05, # probability threshold 
                ... # other arguments passed to MASS::kde2d
){
  kd <- MASS::kde2d(x[,1], x[,2], ...)
  kdx <- fields::interp.surface(kd, x)
  q <- quantile(kdx, p)
  if(is.null(y)){
    kdy <- kdx
  } else {
    kdy <- fields::interp.surface(kd, y)
  }
  kdy >= q
}

# demonstrate use of kernel density function
cali <- cval[rsamp,]
redwood <- cch[spr, c("cwd", "aet")]
1 - mean(kdp(cali, redwood, p = 0.1)) # proportion of redwood occurrences that are outside the 90% CA climate space contour
plot(cali, pch=19, cex=0.1, xlim=c(0,2000),
     col=c("lightblue", "dodgerblue")[as.integer(kdp(cali, p = 0.1))+1])
points(redwood, pch=19, cex=0.5, 
       col=c("darkred", "tomato")[as.integer(kdp(cali, redwood, p = 0.1))+1])



aetAll <- read.csv('results/cchAET.csv')
cwdAll <- read.csv('results/cchCWD.csv')
names(aetAll)
nr <- which(aetAll$name==sp)
aetAll[nr,]
cwdAll[nr,]

plot(cval[rsamp,],pch=19,cex=0.1,col='lightblue',xlim=c(0,2000))
plot(cvch,add=T)
points(cch$cwd[spr],cch$aet[spr],pch=19,cex=1,col='darkgreen')
plot(spch,add=T,lwd=2)
points(cwdAll$pmn[nr],aetAll$pmn[nr],cex=2,pch=19)
points(cwdAll$mpt[nr],aetAll$mpt[nr],cex=2,pch=19,col='orange')


op=par(mfrow=c(2,1))
plot(aet)
points(cch[spr,c('longitude','latitude')],cex=0.5,pch=19,col='darkgreen')
# color ramp for cwd?
plot(cwd)
points(cch[spr,c('longitude','latitude')],cex=0.5,pch=19,col='darkgreen')
par(op)

par(mfrow=c(1,1))
plot(cch$cwd[spr],cch$aet[spr],pch=19,cex=1,col='darkgreen')
points(cwdAll$pmn[nr],aetAll$pmn[nr],cex=2,pch=19)

sv <- vect(as.matrix(cch[spr,]),type='points')
ch <- convHull(sv)

summary(cval[,2])

hist(cwdAll$sd)
hist(aetAll$sd)
