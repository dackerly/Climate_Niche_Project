## climate space
rm(list=ls())
library(terra)

extFromVect <- function(x) {
  x1 <- min(geom(x)[,3])
  x2 <- max(geom(x)[,3])
  y1 <- min(geom(x)[,4])
  y2 <- max(geom(x)[,4])
  return(c(x1,x2,y1,y2))
}

aet <- rast('big_data/CHELSA_climate/CHELSA_annual_AET_1979_2013.tiff')
plot(aet)

cwd <- rast('big_data/CHELSA_climate/CHELSA_annual_CWD_1979_2013.tiff')
plot(cwd)

ppt <- rast('big_data/CHELSA_climate/CHELSA_bio10_12.tiff')
plot(ppt)

tav <- rast('big_data/CHELSA_climate/CHELSA_bio10_01.tiff')/10
plot(tav)

# aggregate up 10x10
aet <- aggregate(aet,10)
cwd <- aggregate(cwd,10)
ppt <- aggregate(ppt,10)
tav <- aggregate(tav,10)

cr <- c(tav,ppt,cwd,aet)
names(cr) <- c('tav','ppt','cwd','aet')

cv <- values(cr)
cv <- data.frame(cv)
dim(cv)

gsel <- which(complete.cases(cv))
cv <- cv[gsel,]
dim(cv)
head(cv)

CA <- vect("data/gis_data/LatLong/CA.shp")
CAext <- ext(extFromVect(CA))

crs(cr)
crs(CA) <- crs(cr)

CAcr <- crop(cr,CAext)
CAcr <- mask(CAcr,CA)
CAcv <- data.frame(values(CAcr))
csel <- complete.cases(CAcv)
CAcv <- CAcv[csel,]
dim(CAcv)

rsamp <- sample(1:nrow(cv),20000)
#pairs(cv[rsamp,])
plot(aet~cwd,data=cv[rsamp,])
points(aet~cwd,data=CAcv,pch=1,col='blue')

dry.low <- which(cv$aet<500 & cv$cwd<500)
x <- rep(0,nrow(cv))
x[dry.low] <- 1
length(x)

xx <- rep(NA,nrow(values(cr)))
xx[gsel] <- x
table(xx,useNA = 'always')
length(xx)

xr <- setValues(aet,xx)
plot(xr)
