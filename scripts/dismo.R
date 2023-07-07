## NOW add species distribution modeling to calculate climatic optima
library(dismo)
#help(dismo)
#vignette('sdm', 'dismo')
library('raster')
#library('rJava')

aet <- raster('data/gis_data/CAaet.tiff')
cwd <- raster('data/gis_data/CAcwd.tiff')
evar <- stack(aet,cwd)
head(cd)

#mx <- maxent(evar,cd[,2:3])
bc <- bioclim(evar,cd[,2:3])