# PWD plot mean relationships
rm(list=ls())
pInfo <- read.csv('data/pwd/plotInfo.csv',row.names=1)[1:50,]
head(pInfo)
dim(pInfo)

tAll <- read.csv('data/pwd/tAll.csv',row.names=1)
head(tAll)
tAll$Basal.Area.13 <- as.numeric(tAll$Basal.Area.13)

cn <- read.csv('data/climNicheData/cchCWD.csv',row.names = 1)
sn <- read.csv('data/pwd/spnames.csv',row.names=1)

names(tAll)
plotList <- sort(unique(tAll$Plot.13))
s2t <- match(tAll$Species.13,sn$Species.13)

head(cn)

i=4
r2 <- c()
slp <- c()
mse <- c()
for (i in 4:ncol(cn)) {
  plotMeans <- rep(NA,50)
  cvar <- cn[s2t,i]
  #length(cvar)
  #length(tAll$Plot.13)
  #length(tAll$Basal.Area.13)
  p=1
  for (p in 1:50) {
    selR <- which(tAll$Plot.13==plotList[p] & is.finite(tAll$Basal.Area.13))
    tmp <- tAll[selR,]
    #length(cvar[selR])==length(tmp$Basal.Area.13)
    plotMeans[p] <- weighted.mean(cvar[selR],tmp$Basal.Area.13,na.rm=T)
  }
  plotMeans
  plot(pInfo$cwd,plotMeans,main=names(cn)[i])
  abline(0,1)
  fit <- lm(plotMeans~pInfo$cwd)
  slp <- c(slp,coef(fit)[2])
  r2 <- c(r2,summary(fit)$r.squared)
  mse <- c(mse,sum(plotMeans-pInfo$cwd)^2/50)
}
names(slp) <- names(cn)[-c(1:3)]
names(r2) <- names(cn)[-c(1:3)]
names(mse) <- names(cn)[-c(1:3)]
slp
r2
mse
plot(mse[c(1,4,8:13)])

## northness
r2 <- c()
slp <- c()
for (i in 4:ncol(cn)) {
  plotMeans <- rep(NA,50)
  cvar <- cn[s2t,i]
  #length(cvar)
  #length(tAll$Plot.13)
  #length(tAll$Basal.Area.13)
  p=1
  for (p in 1:50) {
    selR <- which(tAll$Plot.13==plotList[p] & is.finite(tAll$Basal.Area.13))
    tmp <- tAll[selR,]
    #length(cvar[selR])==length(tmp$Basal.Area.13)
    plotMeans[p] <- weighted.mean(cvar[selR],tmp$Basal.Area.13,na.rm=T)
  }
  plotMeans
  plot(pInfo$northness,plotMeans,main=names(cn)[i])
  abline(0,1)
  fit <- lm(plotMeans~pInfo$cwd)
  slp <- c(slp,coef(fit)[2])
  r2 <- c(r2,summary(fit)$r.squared)
}
names(slp) <- names(cn)[-c(1:3)]
names(r2) <- names(cn)[-c(1:3)]
slp
r2

plot(cwd~northness,data=pInfo)

