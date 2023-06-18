getwd()
dir()

ppt_ramp <- read.csv('TBC3_climate_color_ramps_csv_files/TBC3_color_ramps_PPT.csv',comment.char='#')
head(ppt_ramp)
dim(ppt_ramp)
plot(rep(1,27),1:27,cex=2,pch=15,col=rev(paste('#',ppt_ramp$Hex,sep='')))

