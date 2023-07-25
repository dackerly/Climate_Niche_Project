ct <- read.csv('data/citations/citations.csv')
plot(ct,type='b',log='y',pch=19,lwd=3)
abline(v=2006)
