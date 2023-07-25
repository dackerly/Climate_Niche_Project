# PWD plot mean relationships
rm(list=ls())
library(terra)
library(raster)
cwd <- raster('data/pwd/cwd1981_2010_ave_PWDX.asc')
projection(cwd) <- crs('+proj=aea +datum=NAD83 +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0')

cwd <- terra::rast('data/pwd/cwd1981_2010_ave_PWDX.asc')
#crs(cwd) <- '+proj=aea +datum=NAD83 +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0'
crs(cwd) <- '+proj=aea +datum=NAD83 +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000'
plot(cwd)

pInfo <- read.csv('data/pwd/plotInfo.csv')
head(pInfo)
tail(pInfo)

xx <- vect(pInfo[,c('UTM.E','UTM.N')],geom=c('UTM.E','UTM.N'))
crs(xx) <- "+proj=utm +zone=10 +datum=WGS84"
plot(xx)

xxa <- project(xx,cwd)
plot(xxa)
plot(cwd)
points(xxa)

xcwd <- extract(cwd,xxa)
pInfo$cwd <- xcwd[,2]
write.csv(pInfo,'data/pwd/plotInfo.csv')
