library(parallel)
library(beepr)

# Initiate cluster
no_cores <- detectCores()
cl <- makeCluster(no_cores)


X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=data1[[i]],"logOR"=logOR[[i]]) }

################################################################
###### Mantel Haenszel NMA, no continuity correction ###########
MH=function(X)
{
  biasMH<-c()
  coverageMH<-c()
  MH1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "MH", incr=0, data=X$data)
  MH.res<-data.frame(mean=MH1$TE.fixed[2:NT,1], lowerCI=MH1$lower.fixed[2:NT,1],upperCI=MH1$upper.fixed[2:NT,1])
  biasMH<-c(biasMH, (MH.res$mean-X$logOR))
  MH.res$cover<-(MH.res$lowerCI<X$logOR)&(MH.res$upperCI>X$logOR)
  coverageMH<-c(coverageMH, MH.res$cover)
return(list("bias"=biasMH,"cov"=coverageMH))
  
}
clusterExport(cl,"X1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "MH")
clusterEvalQ(cl, {library(netmeta)})
l2=parLapply(cl,1:N.sim, function(x) MH(X1[[x]]))

biasMH=c()
for (i in 1:N.sim){  biasMH=c(biasMH, l2[[i]]$bias)}
mean(biasMH)

coverageMH=c()
for (i in 1:N.sim){  coverageMH=c(coverageMH, l2[[i]]$cov)}
mean(coverageMH)
beep(sound=3)
################################################################


################################################################
###### NCH NMA, no continuity correction ###########
NCH=function(X)
{
  biasNCH<-c()
  coverageNCH<-c()
  NCH1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "NCH", incr=0, data=X$data)
  NCH.res<-data.frame(mean=NCH1$TE.fixed[2:NT,1], lowerCI=NCH1$lower.fixed[2:NT,1],upperCI=NCH1$upper.fixed[2:NT,1])
  biasNCH<-c(biasNCH, (NCH.res$mean-X$logOR))
  NCH.res$cover<-(NCH.res$lowerCI<X$logOR)&(NCH.res$upperCI>X$logOR)
  coverageNCH<-c(coverageNCH, NCH.res$cover)
  return(list("bias"=biasNCH,"cov"=coverageNCH))
  
}
clusterExport(cl,"X1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "NCH")
clusterEvalQ(cl, {library(netmeta)})
l3=parLapply(cl,1:N.sim, function(x) NCH(X1[[x]]))


biasNCH=c()
for (i in 1:N.sim){  biasNCH=c(biasNCH, l3[[i]]$bias)}
mean(biasNCH)

coverageNCH=c()
for (i in 1:N.sim){  coverageNCH=c(coverageNCH, l3[[i]]$cov)}
mean(coverageNCH)
beep(sound=3)
################################################################


################################################################
###### IV NMA, 0.5 continuity correction ###########
IV=function(X)
{
  biasIV.FE<-c()
  coverageIV.FE<-c()
  biasIV.RE<-c()
  coverageIV.RE<-c()
  
  IV1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "Inverse", incr=0.5 ,cc.pooled=T, allstudies=T ,sm="OR",data=X$data)
  IV.FE.res<-data.frame(mean=IV1$TE.fixed[2:NT,1], lowerCI=IV1$lower.fixed[2:NT,1],upperCI=IV1$upper.fixed[2:NT,1])
  biasIV.FE<-c(biasIV.FE, (IV.FE.res$mean-X$logOR))
  IV.FE.res$cover<-(IV.FE.res$lowerCI<X$logOR)&(IV.FE.res$upperCI>X$logOR)
  coverageIV.FE<-c(coverageIV.FE, IV.FE.res$cover)
  
  IV.RE.res<-data.frame(mean=IV1$TE.random[2:NT,1], lowerCI=IV1$lower.random[2:NT,1],upperCI=IV1$upper.random[2:NT,1])
  biasIV.RE<-c(biasIV.RE, (IV.RE.res$mean-X$logOR))
  IV.RE.res$cover<-(IV.RE.res$lowerCI<X$logOR)&(IV.RE.res$upperCI>X$logOR)
  coverageIV.RE<-c(coverageIV.RE, IV.RE.res$cover)
  
  return(list("biasFE"=biasIV.FE,"biasRE"=biasIV.RE,"covFE"=coverageIV.FE, "covRE"=coverageIV.RE))
  
}
clusterExport(cl,"X1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "IV")
clusterEvalQ(cl, {library(netmeta)})
l1=parLapply(cl,1:N.sim, function(x) IV(X1[[x]]))


biasIV.FE=c()
for (i in 1:N.sim){  biasIV.FE=c(biasIV.FE, l1[[i]]$biasFE)}
mean(biasIV.FE)

coverageIV.FE=c()
for (i in 1:N.sim){  coverageIV.FE=c(coverageIV.FE, l1[[i]]$covFE)}
mean(coverageIV.FE)

biasIV.RE=c()
for (i in 1:N.sim){  biasIV.RE=c(biasIV.RE, l1[[i]]$biasRE)}
mean(biasIV.RE)

coverageIV.RE=c()
for (i in 1:N.sim){  coverageIV.RE=c(coverageIV.RE, l1[[i]]$covRE)}
mean(coverageIV.RE)

stopCluster(cl)
beep(sound=3)


################################################################


################
# RESULTS
################
# IV.FE
round(mean(biasIV.FE), digits=2)
round(mean(abs(biasIV.FE)), digits=2)
 round(mean(coverageIV.FE), digits=3)*100
# IV.RE
round(mean(biasIV.RE), digits=2)
round(mean(abs(biasIV.RE)), digits=2)
round(mean(coverageIV.RE), digits=3)*100
# MH
round(mean(biasMH), digits=2)
round(mean(abs(biasMH)), digits=2)
round(mean(coverageMH), digits=3)*100
# NCH
round(mean(biasNCH), digits=2)
round(mean(abs(biasNCH)), digits=2)
round(mean(coverageNCH), digits=3)*100


