# This analysis uses David Ackerly's species-level summary statistics
# of CWD and AET niches (e.g., various kinds of niche means, optima, etc.)
# to compute community weighted climatic niche means for tree assemblages
# in California FIA plots used in Rosenblad et al. (2023) PNAS.
# Then we see which niche mean version is best predicted by macroclimate.

options(scipen=9999)

library(foreach)
library(doParallel)
detectCores() # this tells you how many are available in your current environment
cores <- 10 # ADJUST TO THE NUMBER OF PROCESSOR CORES YOU WANT TO USE FOR PARALLELIZATION

# read previously prepped subplot data
subplots <- readRDS("./data/kr_fia_data_products/current_subplots.RDS")

# zoom out to plot level
plots <- unique(subplots[c("LAT_LON", "LAT", "LON", "STATECD")])
rm(subplots)
gc()

# keep california only
plots <- subset(plots, STATECD==6 & !is.na(STATECD))

# read previously prepped tree data
data <- readRDS("./data/kr_fia_data_products/data2.RDS")
data <- subset(data, LAT_LON%in%plots$LAT_LON)
gc()

# load in david ackerly's niche means 
cwdniches <- read.csv("./data/climNicheData/cchCWD.csv")
aetniches <- read.csv("./data/climNicheData/cchAET.csv")

# restructure so that the vector for each climate variable's summary
# statistic has the climate variable in the vector name
cwdniches$climVar <- NULL
names(cwdniches)[names(cwdniches)!="name"] <- paste("cwd",
                                                names(cwdniches)[names(cwdniches)!="name"],
                                                sep="_")
aetniches$climVar <- NULL
names(aetniches)[names(aetniches)!="name"] <- paste("aet",
                                                names(aetniches)[names(aetniches)!="name"],
                                                sep="_")

data$name <- paste(data$GENUS, data$SPECIES, sep=" ")
data <- merge(data, cwdniches, by="name", all.x=TRUE, all.y=FALSE, sort=FALSE)
data <- merge(data, aetniches, by="name", all.x=TRUE, all.y=FALSE, sort=FALSE)

# calculate basal area from diameter
data$basal_area <- (pi/(4*144)) * ((data$DIA)^2)

# filter out trees that fall
# within a macroplot but not within a standard subplot. this is an
# issue in the pnw but not the rm region.
data <- subset(data, DIST<=24)

# filter out a small subset of trees subject to procedural errors
# and changes:
data <- subset(data, RECONCILECD==1 | RECONCILECD==2 | is.na(RECONCILECD))

# identify vector names in "data" for which to compute community weighted values
library(stringr)
cwdnames <- names(data)[str_starts(string=names(data), pattern="cwd")]
aetnames <- names(data)[str_starts(string=names(data), pattern="aet")]
varnames <- c(cwdnames, aetnames)

# compute community weighted niche means (and optima, etc.)
registerDoParallel(cores) # AGAIN, THE NUMBER OF CORES MAY NEED TO BE ADJUSTED
d <- foreach(i=1:nrow(plots), .combine=rbind) %dopar% {
  trees <- data[data$LAT_LON==plots[i,"LAT_LON"],]
  trees <- subset(trees, STATUSCD==1 & INVYR==max(trees$INVYR))
  row <- c()
  for(j in 1:length(varnames)){
    row[j] <- sum(trees[c(varnames[j])]*trees$basal_area)/sum(trees$basal_area)
  }
  row
}
d <- as.data.frame(d)
names(d) <- varnames
plots <- cbind(plots, d)

rm(d, data)
gc()

# some plots had no living trees in the most recent census, so exclude these
plots <- na.exclude(plots)

# load and merge on chelsa climate data for each plot
library(raster)
coordinates(plots) <- c("LON", "LAT")
crs(plots) <- CRS("+proj=longlat +datum=NAD83")

cwd <- raster("./data/gis_data/CAcwd.tiff")
plots <- spTransform(plots, crs(cwd))
plots$cwd <- extract(cwd, plots)

aet <- raster("./data/gis_data/CAaet.tiff")
plots <- spTransform(plots, crs(aet))
plots$aet <- extract(aet, plots)

plots <- as.data.frame(plots)


# which community cwd statistic is best predicted by cwd?
cwdsumm <- data.frame(commstat=cwdnames, ols_es=NA, ols_r2=NA, mse1to1=NA)
for(i in 1:nrow(cwdsumm)){
  y <- plots[cwdsumm[i,"commstat"]]
  mod <- lm(y[,1] ~ plots$cwd)
  cwdsumm[i, "ols_es"] <- mod$coefficients[2]
  cwdsumm[i, "ols_r2"] <- summary(mod)$adj.r.squared
  cwdsumm[i, "mse1to1"] <- mean((y[,1]-plots$cwd)^2)
}

# X and N look like they're not meant to be niche means
cwdsumm <- subset(cwdsumm, commstat!= "cwd_N" & commstat!= "cwd_X")

# prep summary data for plotting
library(stringr)
cwdsumm$commstat <- str_remove_all(string=cwdsumm$commstat,
                                   pattern="cwd_")


# which community aet statistic is best predicted by aet?
aetsumm <- data.frame(commstat=aetnames, ols_es=NA, ols_r2=NA, mse1to1=NA)
for(i in 1:nrow(aetsumm)){
  y <- plots[aetsumm[i,"commstat"]]
  mod <- lm(y[,1] ~ plots$aet)
  aetsumm[i, "ols_es"] <- mod$coefficients[2]
  aetsumm[i, "ols_r2"] <- summary(mod)$adj.r.squared
  aetsumm[i, "mse1to1"] <- mean((y[,1]-plots$aet)^2)
}

# X and N look like they're not meant to be niche means
aetsumm <- subset(aetsumm, commstat!= "aet_N" & commstat!= "aet_X")

# prep summary data for plotting
aetsumm$commstat <- str_remove_all(string=aetsumm$commstat,
                                   pattern="aet_")

# order by effect size in plots
cwdsumm$commstat <- factor(cwdsumm$commstat, levels=unique(cwdsumm$commstat[order(-cwdsumm$ols_es)]), ordered=TRUE)
aetsumm$commstat <- factor(aetsumm$commstat, levels=unique(aetsumm$commstat[order(-aetsumm$ols_es)]), ordered=TRUE)

### plots
library(ggplot2)
cwdsumm_plot <- ggplot(cwdsumm, aes(x=commstat, y=ols_es, color=ols_r2, size=mse1to1))+
  ggtitle("CWD")+
  geom_hline(yintercept=1, lty="dotted")+
  geom_hline(yintercept=0)+
  xlab("Community Weighted Niche Metric")+
  ylab("Effect Size of Macroclimatic Predictor")+
  geom_point()+
  scale_color_viridis_c(name="R²",
                        limits=c(
                          min(c(cwdsumm$ols_r2, aetsumm$ols_r2)),
                          max(c(cwdsumm$ols_r2, aetsumm$ols_r2))
                        ))+
  scale_y_continuous(limits=c(
    min(c(cwdsumm$ols_es, aetsumm$ols_es)),
    max(c(cwdsumm$ols_es, aetsumm$ols_es))
  ))+
  scale_size_continuous(name="1:1\nMSE")+
  theme_bw()+
  theme(legend.position="none")

aetsumm_plot <- ggplot(aetsumm, aes(x=commstat, y=ols_es, color=ols_r2, size=mse1to1))+
  ggtitle("AET")+
  geom_hline(yintercept=1, lty="dotted")+
  geom_hline(yintercept=0)+
  xlab("Community Weighted Niche Metric")+
  geom_point()+
  scale_color_viridis_c(name="R²",
                        limits=c(
                          min(c(cwdsumm$ols_r2, aetsumm$ols_r2)),
                          max(c(cwdsumm$ols_r2, aetsumm$ols_r2))
                        ))+
  scale_y_continuous(limits=c(
    min(c(cwdsumm$ols_es, aetsumm$ols_es)),
    max(c(cwdsumm$ols_es, aetsumm$ols_es))
  ))+
  scale_size_continuous(name="1:1\nMSE")+
  theme_bw()+
  theme(axis.title.y=element_blank())

cwdsumm_plot
ggsave("./results/kr_fia_results/cch_niches_fia_plots_summary_cwd.png", height=5, width=5)
aetsumm_plot
ggsave("./results/kr_fia_results/cch_niches_fia_plots_summary_aet.png", height=5, width=5)


### scatterplots for each community mean statistic
cwdnames <- subset(cwdnames, cwdnames!="cwd_X" & cwdnames!="cwd_N")
cwdplots <- plots[c("cwd", cwdnames)]

aetnames <- subset(aetnames, aetnames!="aet_X" & aetnames!="aet_N")
aetplots <- plots[c("aet", aetnames)]

for(i in 2:ncol(cwdplots)){
  ggplot(plots, aes(x=cwd, y=cwdplots[,i]))+
    geom_point()+
    geom_smooth(method=lm)+
    xlab("CWD (mm)")+
    ylab(names(cwdplots)[i])+
    theme_bw()
  ggsave(paste("./results/kr_fia_results/cch_niches_", names(cwdplots)[i], "_vs_cwd_scatter.png", sep=""), width=5, height=5)
}

for(i in 2:ncol(aetplots)){
  ggplot(plots, aes(x=aet, y=aetplots[,i]))+
    geom_point()+
    geom_smooth(method=lm)+
    xlab("aet (mm)")+
    ylab(names(aetplots)[i])+
    theme_bw()
  ggsave(paste("./results/kr_fia_results/cch_niches_", names(aetplots)[i], "_vs_aet_scatter.png", sep=""), width=5, height=5)
}
