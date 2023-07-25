# overlapping niches, and mean niche values at one location (community)
rm(list=ls())

# assume equal niche breadth, and Guassian niche

evmin=0
evmax=10

nb <- 1
mn <- seq(0,10,by=0.5)
length(mn)
ev <- seq(evmin,evmax,by=0.05)
eva <- seq(0,evmin,by=0.05)
lev <- length(ev)
nreps <- 100

res <- data.frame(rep=1:nreps,pmn=NA,mwm=NA,mwa=NA,mpt=NA,mat=NA,opt=NA)

i=1
for (i in 1:nreps) {
  np <- dnorm(ev,mn[11],nb)
  npa <- dnorm(eva,mn[11],nb)
  pa <- rbinom(lev,1,dnorm(ev,mn[11],nb))
  p <- which(pa==1)
  psa <- runif(100,evmin,evmax)
  
  plot(ev,np,type='l',ylim=c(0,0.6),xlim=c(0,10),lwd=2)
  points(eva,npa,type='l',ylim=c(0,0.6),xlim=c(2,8),lty=2)
  points(ev[p],rep(0.6,length(p)),pch=19,col='red')
  points(psa,rep(0.05,length(psa)),pch=19)
  abline(v=mn[11])
  
  # now estimate niche with gaussian model
  allpa <- data.frame(e=c(ev[p],psa),pa=0)
  dim(allpa)
  allpa
  allpa$pa[1:length(p)] <- 1
  allpa$e2 <- allpa$e^2
  #plot(allpa)
  fit <- glm(pa~e+e2,data=allpa,family = 'binomial')
  allpa$pval <- fitted(fit,type='response')
  points(pval~e,data=allpa[order(allpa$e),],pch=10,type='l')
  
  (res$pmn[i] <- mean(ev[p]))
  (res$mwm[i] <- weighted.mean(ev[p],allpa$pval[allpa$pa==1]))
  (res$mwa[i] <- weighted.mean(allpa$e,allpa$pval))
  (res$mpt[i] <- ev[p[which.max(allpa$pval[allpa$pa==1])]])
  (res$mat[i] <- allpa$e[which.max(allpa$pval)])
  (res$opt[i] <- mn[11])
}
res
head(res)
round(apply(res[,-1],2,mean,na.rm=T),2)
round(apply(res[,-1],2,sd,na.rm=T),3)
### community figure

cc <- data.frame(mn=mn,dens=NA)
plot(ev,np,type='n',ylim=c(0,0.5))
for (i in 1:length(mn)) {
  cc$dens[i] <- dnorm(ev,mn[i])
  points(ev,dnorm(ev,mn[i]),nb,type='l')
}
abline(v=c(0.5,5))
weighted.mean(cc$mn,cc$dens)
