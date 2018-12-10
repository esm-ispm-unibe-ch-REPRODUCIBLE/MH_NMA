library(netmeta)


### Fit MH-NMA, continuity correction=0 #####
MH.res=list()
biasMH=c()
coverageMH=c()
for (i in 1:N.sim){
  print(paste("Fitting NMA number ", i))
  MH=netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "MH", incr=0, data=data1[[i]])
  MH.res[[i]]=data.frame(mean=MH$TE.fixed[2:NT,1], lowerCI=MH$lower.fixed[2:NT,1],upperCI=MH$upper.fixed[2:NT,1])
  biasMH=c(biasMH, (MH.res[[i]]$mean-logOR))
  MH.res[[i]]$cover=(MH.res[[i]]$lowerCI<logOR)&(MH.res[[i]]$upperCI>logOR)
  coverageMH=c(coverageMH, MH.res[[i]]$cover)
}
mean(biasMH)
mean(abs(biasMH))
median(abs(biasMH))
mean(coverageMH)
#####################


### Fit IV-NMA, continuity correction=0.5 #######
IV.res=list()
biasIV=c()
coverageIV=c()
for (i in 1:N.sim){
  print(paste("Fitting NMA number ", i))
  IV=netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "Inverse", incr=0.5, data=data1[[i]], sm="OR")
  IV.res[[i]]=data.frame(mean=IV$TE.fixed[2:NT,1], lowerCI=IV$lower.fixed[2:NT,1],upperCI=IV$upper.fixed[2:NT,1])
  biasIV=c(biasIV, (IV.res[[i]]$mean-logOR))
  IV.res[[i]]$cover=(IV.res[[i]]$lowerCI<logOR)&(IV.res[[i]]$upperCI>logOR)
  coverageIV=c(coverageIV, IV.res[[i]]$cover)
}
mean(biasIV)
mean(abs(biasIV))
median(abs(biasIV))
mean(coverageIV)
###################


### Fit NCH, continuity correction=0 #####
NCH.res=list()
biasNCH=c()
coverageNCH=c()
for (i in 1:N.sim){
  print(paste("Fitting NMA number ", i))
  NCH=netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "NCH", incr=0, data=data1[[i]])
  NCH.res[[i]]=data.frame(mean=NCH$TE.fixed[2:NT,1], lowerCI=NCH$lower.fixed[2:NT,1],upperCI=NCH$upper.fixed[2:NT,1])
  biasNCH=c(biasNCH, (NCH.res[[i]]$mean-logOR))
  NCH.res[[i]]$cover=(NCH.res[[i]]$lowerCI<logOR)&(NCH.res[[i]]$upperCI>logOR)
  coverageNCH=c(coverageNCH, NCH.res[[i]]$cover)
}
mean(biasNCH)
mean(abs(biasNCH))
median(abs(biasNCH))
mean(coverageNCH)
###################

