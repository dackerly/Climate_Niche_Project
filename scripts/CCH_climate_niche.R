###### FILE LOCATIONS NOT UPDATED YET - CLIMATE DATA NOT IN PROJECT

library(raster)

#read in specimen data
cch <- readRDS('/Users/david/Documents/Projects/CalDiversity/PhylodiversityProject/PD_project/Nov15_spatial/caldiv_data/California_Species_clean_All_epsg_3310.Rdata')
dim(cch)
names(cch)
cch.spp <- unique(cch$current_name_binomial)
cch$current_name_binomial <- as.character(cch$current_name_binomial)

#read in climate layers
jja <- raster('/Users/david/Documents/Data/BCM/CA2014/WY30/jja1951_1980W.grd')
djf <- raster('/Users/david/Documents/Data/BCM/CA2014/WY30/djf1951_1980W.grd')
cwd <- raster('/Users/david/Documents/Data/BCM/CA2014/WY30/cwd1951_1980_ave_HST.asc')
cwd[cwd<0] <- NA

# TRY ONE SPECIES
spp <- unique(cch$current_name_binomial)[1]
spp

cch1 <- subset(cch,cch$current_name_binomial==spp)
dim(cch1)
plot(djf)
points(cch1[,c('x_epsg_3310','y_epsg_3310')],pch='.')

djf1 <- extract(djf,cch1[,c('x_epsg_3310','y_epsg_3310')])
head(djf1)
cwd1 <- extract(cwd,cch1[,c('x_epsg_3310','y_epsg_3310')])
head(cwd1)
# GOOD

# Run individual data sets
# GLORIA
gl <- read.delim('/Users/david/Documents/Projects/CalDiversity/ClimaticNicheAnalysis/CNG_2016/species_lists/Gloria/gloria-species-list.txt',header=T,as.is=T)
head(gl)

gl.match <- gl$species %in% cch.spp
table(gl.match)
gl$species[!gl.match]

# find names in cch species list
cch.spp[grep('Carex d',cch.spp)]

cch_gl <- subset(cch,cch$current_name_binomial %in% gl$species)
dim(cch_gl)
head(cch_gl)
djfx <- extract(djf,cch_gl[,c('x_epsg_3310','y_epsg_3310')])
jjax <- extract(jja,cch_gl[,c('x_epsg_3310','y_epsg_3310')])
cwdx <- extract(cwd,cch_gl[,c('x_epsg_3310','y_epsg_3310')])

djf_spp <- tapply(djfx,cch_gl$current_name_binomial,mean,na.rm=T)
jja_spp <- tapply(jjax,cch_gl$current_name_binomial,mean,na.rm=T)
cwd_spp <- tapply(cwdx,cch_gl$current_name_binomial,mean,na.rm=T)

clim_spp <- data.frame(species=names(djf_spp),djf=djf_spp,jja=jja_spp,cwd=cwd_spp)
clim_spp[1:10,]

write.csv(clim_spp,'/Users/david/Documents/Projects/CalDiversity/ClimaticNicheAnalysis/CNG_2016/species_lists/Gloria/clim_means.csv')

# PEPPERWOOD
pwdh <- read.csv('/Users/david/Documents/Projects/CalDiversity/ClimaticNicheAnalysis/CNG_2016/species_lists/PWD_native/PWD_native.csv',as.is=T,row.names = 1)
head(pwdh)

pwdh.match <- pwdh$binomial %in% cch$current_name_binomial
table(pwdh.match)
pwdh$binomial[!pwdh.match]

#write.csv(cbind(pwdh[,1:6],pwdh.match),'/Users/david/Documents/Projects/CalDiversity/ClimaticNicheAnalysis/CNG_2016/species_lists/PWD_native/PWD_herbaceous.csv')

cch_pwdh <- subset(cch,cch$current_name_binomial %in% pwdh$binomial)
dim(cch_pwdh)
djfx <- extract(djf,cch_pwdh[,c('x_epsg_3310','y_epsg_3310')])
jjax <- extract(jja,cch_pwdh[,c('x_epsg_3310','y_epsg_3310')])
cwdx <- extract(cwd,cch_pwdh[,c('x_epsg_3310','y_epsg_3310')])

djf_spp <- tapply(djfx,cch_pwdh$current_name_binomial,mean,na.rm=T)
jja_spp <- tapply(jjax,cch_pwdh$current_name_binomial,mean,na.rm=T)
cwd_spp <- tapply(cwdx,cch_pwdh$current_name_binomial,mean,na.rm=T)

codes <- pwdh$code[match(names(djf_spp),pwdh$binomial)]
codes

clim_spp <- data.frame(code=codes,djf=djf_spp,jja=jja_spp,cwd=cwd_spp)
head(clim_spp)

write.csv(clim_spp,'/Users/david/Documents/Projects/CalDiversity/ClimaticNicheAnalysis/CNG_2016/species_lists/PWD_native/clim_means.csv')

# TRY ENTIRE DATA SET
t1 <- Sys.time()
djfx <- extract(djf,cch[1:1000,c('x_epsg_3310','y_epsg_3310')])
t2 <- Sys.time()
t2-t1
