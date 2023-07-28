rm(list=ls())
library('terra')
library('dismo')
library('fields')
library('alphahull')
library('maptools')
library('rgeos')
library('ggplot2')
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
dim(cval.all)
xy <- xyFromCell(aet,1:nrow(cval.all))
dim(xy)
cval.all <- cbind(xy,cval.all)
head(cval.all)

cval <- cval.all[which(complete.cases(cval.all)),]
dim(cval)
rsamp <- sample(1:nrow(cval),50000)
cvch <- terra::convHull(cval[,3:4])

# insert alpha hull and contour density

# extract climate for cch
cch$aet <- extract(aet,cch[,c('longitude','latitude')])[,2]
cch$cwd <- extract(cwd,cch[,c('longitude','latitude')])[,2]
vCols <- c(which(names(cch)=='aet'),which(names(cch)=='cwd'))sp <- 'Sequoia sempervirens'
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
1 - mean(kdp(cali[,c(3:4)], redwood, p = 0.1)) # proportion of redwood occurrences that are outside the 90% CA climate space contour

head(cali)
plot(cali[,3:4], pch=19, cex=0.1, xlim=c(0,2000),
     col=c("lightblue", "dodgerblue")[as.integer(kdp(cali[,3:4], p = 0.1))+1])
points(redwood, pch=19, cex=0.5, 
       col=c("darkred", "tomato")[as.integer(kdp(cali[,3:4], redwood, p = 0.1))+1])

# where are the 10% points?
plot(cali[,1:2], pch=19, cex=0.1,
     col=c("lightblue", "dodgerblue")[as.integer(kdp(cali[,3:4], p = 0.1))+1])

# function to calculate distance to alpha hull
alpha_dist <- function(x, # 2-column matrix used to define alpha hull
                       y = NULL, # optional matrix for which hull distance is to be calculated; if not supplied, distances are calcualted for x
                       alpha = 0.5 # radius of cookie cutter (in standard deviations of x)
){
  
  # unique, scaled values
  xs <- scale(x)
  dupe <- duplicated(xs)
  xu <- xs[!dupe,]
  if(is.null(y)) y <- x
  for(i in 1:2) y[,i] <- (y[,i] - mean(x[,i])) / sd(x[,i])
  
  # alpha hull
  ah <- alphahull::ahull(xu, alpha = alpha)
  
  # point-hull distances
  p <- as.data.frame(y)
  coordinates(p) <- colnames(p)
  source("https://raw.githubusercontent.com/matthewkling/range-edges/master/code/utilities.R")
  sp <- ah2sp(ah)
  as.vector(rgeos::gDistance(p, as(sp, "SpatialLines"), byid=T))
}

# demonstrate alpha_dist function with two different radii
cali$alpha_0.5 <- alpha_dist(cali[,3:4], alpha = .5)
cali$alpha_0.05 <- alpha_dist(cali[,3:4], alpha = .05)
dim(cali)
head(cali)

cal <- as.data.frame(cali)  
cal <- tidyr::gather(cal, alpha, edge_dist, alpha_0.05, alpha_0.5)
head(cal)
hist(cal$edge_dist)

ggplot(cal, aes(cwd, aet, color = 1.5-edge_dist)) +
  facet_wrap(~alpha) +
  geom_point() +
  scale_color_viridis_c()

ggplot(cal, aes(x, y, color = 1.5-edge_dist)) +
  facet_wrap(~alpha) +
  geom_point() +
  scale_color_viridis_c()

hist(cali$alpha_0.05)
hist(cali$alpha_0.5)
plot(alpha_0.05~alpha_0.5,data=cali)

rw_dist_0.5 <- alpha_dist(cali[,3:4], redwood, alpha = .5)
rw_dist_0.05 <- alpha_dist(cali[,3:4], redwood, alpha = .05)

al <- cali$alpha_0.5
rw <- rw_dist_0.5

head(cali)
plot(y~x,data=cali,col=c('lightblue','dodgerblue')[cut(al,breaks=c(min(al),quantile(al,0.9),max(al)))],asp=1)

## redwood example was helpful demo. Now do it for all species
cch2 <- cch[which(complete.cases(cch[,c('cwd','aet')])),]
cch2$alpha_dist_0.5 <- alpha_dist(cali[,3:4], cch2[,c('cwd','aet')], alpha = .5)
saveRDS(cch2$alpha_dist_0.5,'data/cch2_alpha_dist_0.5')

allSp <- sort(unique(cch2$current_name_binomial))
spMar <- data.frame(name=allSp,N=NA,Nm_0.1=NA)
length(allSp)

spMar$N <- tapply(cch2$current_name_binomial,cch2$current_name_binomial,length)

ltal <- function(x) length(which(x<quantile(al,0.1)))
spMar$Nm_0.1 <- tapply(cch2$alpha_dist_0.5,cch2$current_name_binomial,ltal)

spMar$Pm_0.1 <- spMar$Nm_0.1/spMar$N
head(spMar)
hist(spMar$Pm_0.1)
dim(spMar)

# number of species that have >= 5% of their observations in the 10% percentile edge of the alpha hull
length(which(spMar$Pm_0.1>=0.05))

# number of species that have >0 and <5% of their observations in the 10% percentile edge of the alpha hull
length(which(spMar$Pm_0.1<0.05 & spMar$Pm_0.1>0))

# number of species that have 0% of their observations in the 10% percentile edge of the alpha hull
length(which(spMar$Pm_0.1==0))

# show that hull can be fit to one data set and distance calculated for a second dataset

str(rw_dist)
str(cal)

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
