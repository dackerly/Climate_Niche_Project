rm(list=ls())
aet <- read.csv('data/quercus/quercusAET.csv')
cwd <- read.csv('data/quercus/quercusCWD.csv')
p50 <- read.csv('data/quercus/quercusP50.csv')

all(aet$name==cwd$name)
all(aet$name==p50$name)

ymn <- range(c(p50$Stem.P50,p50$Leaf.P50),na.rm=T)
which.min(p50$Stem.P50)
p50[which.min(p50$Stem.P50),]

names(cwd)
nct <- c(5,7,8,9,12:17)
# -15 removes Q. pacifica outlier
for (i in 1:length(nct)) {
  mvar <- nct[i]
  cval <- cor(cwd[-15,mvar],p50$Stem.P50[-15],use='pair')
  print(cval)
  plot(cwd[,mvar],p50$Stem.P50,pch=1,col='blue',ylim=ymn,main=paste(names(cwd)[mvar],round(cval,3)),xlab='cwd niche')
  points(cwd[-15,mvar],p50$Stem.P50[-15],pch=19,col='blue')
  abline(lm(p50$Stem.P50[-15]~cwd[-15,mvar]))
  #points(cwd[,mvar],p50$Leaf.P50,pch=19,col='red')
}

# -15 removes Q. pacifica outlier
for (i in 1:length(nct)) {
  mvar <- nct[i]
  cval <- cor(aet[-15,mvar],p50$Stem.P50[-15],use='pair')
  print(cval)
  plot(aet[,mvar],p50$Stem.P50,pch=19,col='blue',ylim=ymn,main=paste(names(cwd)[mvar],round(cval,3)),xlab='aet niche')
  #points(cwd[,mvar],p50$Leaf.P50,pch=19,col='red')
}

{
  op=par(mfrow=c(2,2),mar=c(5,5,2,2)) 
  nct <- c(5,15)
  # -15 removes Q. pacifica outlier
  i=1
  for (i in 1:length(nct)) {
    mvar <- nct[i]
    cval <- cor(cwd[-15,mvar],p50$Stem.P50[-15],use='pair')
    print(cval)
    plot(cwd[,mvar],p50$Stem.P50,pch=1,col='blue',ylim=ymn,main=paste(names(cwd)[mvar],round(cval,3)),xlab='cwd niche',cex=2)
    points(cwd[-15,mvar],p50$Stem.P50[-15],pch=19,col='blue',cex=2)
    abline(lm(p50$Stem.P50[-15]~cwd[-15,mvar]))
    #points(cwd[,mvar],p50$Leaf.P50,pch=19,col='red')
  }
  
  # -15 removes Q. pacifica outlier
  for (i in 1:length(nct)) {
    mvar <- nct[i]
    cval <- cor(aet[-15,mvar],p50$Stem.P50[-15],use='pair')
    print(cval)
    plot(aet[,mvar],p50$Stem.P50,pch=1,col='blue',ylim=ymn,main=paste(names(cwd)[mvar],round(cval,3)),xlab='aet niche',cex=2)
    points(aet[-15,mvar],p50$Stem.P50[-15],pch=19,col='blue',cex=2)
    abline(lm(p50$Stem.P50[-15]~aet[-15,mvar]))
    #points(cwd[,mvar],p50$Leaf.P50,pch=19,col='red')
  }
  par(op)
}

# check how my clim niches compare to Rob's
for (i in 1:length(nct)) {
  mvar <- nct[i]
  cval <- cor(cwd[,mvar],p50$CWD_rs,use='pair')
  print(cval)
  plot(cwd[,mvar],p50$CWD_rs,pch=19,col='blue',main=paste(names(cwd)[mvar],round(cval,3)),xlab='cwd niche')
  #points(cwd[,mvar],p50$Leaf.P50,pch=19,col='red')
}

cor(p50$CWD_rs,p50$Stem.P50,use='pair')
cor(p50$CWD_rs[-15],p50$Stem.P50[-15],use='pair')
