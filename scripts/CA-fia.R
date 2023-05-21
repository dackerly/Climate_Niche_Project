rm(list=ls())
library(raster)
library(picante)

# check how matrix2sample works
data(phylocom)
head(phylocom)
tmp <- matrix2sample(phylocom$sample)
head(tmp)
#####

## read in plot locations for CA FIA plots
plot <- read.csv('data/fia_data/CA_PLOT.csv')

# this covers entire state
plot(plot$LON,plot$LAT,asp=1,cex=0.2)

# read in plots with tree species data
st <- read.csv('data/fia_data/CA_SITETREE.csv')
dim(st)
head(st)
st$ba <- pi*(st$DIA/2)^2

head(sort(unique(st$PLOT)))
tail(sort(unique(st$PLOT)))

stpba <- tapply(st$ba,list(st$PLOT,st$SPCD),sum)
dim(stpba)
head(stpba)
head(rownames((stpba)))
tail(rownames((stpba)))

pba <- matrix2sample(stpba)
head(pba)
tail(pba)
hist(table(pba$plot))
dim(pba)

pba$lon <- plot$LON[match(pba$plot,plot$PLOT)]
pba$lat <- plot$LAT[match(pba$plot,plot$PLOT)]
plot(pba$lon,pba$lat,asp=1)

## read in cwd data and transfer to st
cwd <- raster('big_data/CHELSA_climate/CHELSA_annual_CWD_1979_2013.tiff')
cwd

pba$cwd <- extract(cwd,data.frame(pba$lon,pba$lat))
hist(pba$cwd)

# calculate simple mean cwd by species
sp_cwd <- tapply(pba$cwd,pba$id,mean)
head(sp_cwd)
tail(sp_cwd)
length(sp_cwd)

# now transfer sp_cwd back to st
pba$sp_cwd <- sp_cwd[match(pba$id,as.numeric(names(sp_cwd)))]
head(pba)

plot(pba$cwd,pba$sp_cwd)
abline(0,1)

# throwing an error
plot_cwm_cwd <- tapply(pba$sp_cwd,pba$plot,weighted.mean,pba$abund)

  
############### sandbox
tmp <- read.csv('big_data/CAFIA/CA_TREE.csv')
dim(tmp)
names(tmp)
sort(unique(tmp$SPCD))
table(tmp$INVYR)
length(unique(tmp$PLT_CN))



# confirm all plots are in PLOT database
all(st$PLOT %in% plot$PLOT)

# transfer lon, lat data from plot to st
st$LON <- plot$LON[match(st$PLOT,plot$PLOT)]
st$LAT <- plot$LAT[match(st$PLOT,plot$PLOT)]
st$lonlat <- paste(st$LON,st$LAT)
plot(st$LON,st$LAT,asp=1)
length(unique(st$PLOT))
length(unique(st$lonlat))
hist(table(st$PLOT))

sort(unique(st$SPCD))
write.csv(sort(unique(st$SPCD)),'data/fia_data/CAspcd.csv')

