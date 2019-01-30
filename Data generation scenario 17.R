set.seed(42) # the answer to the ultimate question of life, the universe and everything
N.sim=1000
library(MASS)


###### Scenario 17 ######
# FIXED effects three arm studies
name="Scenario-17.RData"
data1=list()
data2=list()
logOR=list()
OR=list()
NT=3 #### number of treatments in the network
NS=8 #### number of studies per comparison
tau=0 ####  heterogeneity SD
Npmin=30 #### minimum number of patients per arm
Npmax=60 #### maximum number of patients per arm

### define treatment indices
t1=c(1:3)
t=(rep(t1,NS))
stud=rep(1:NS, each=3)
for (i in 1:N.sim)
{ 
  data1[[i]]=data.frame(stud,t)
}
### define patients per treatment arm
for (i in 1:N.sim)
{   
  data1[[i]]$n=rep(round(runif(NS,Npmin,Npmax)),each=NT)
}

narms=length(data1[[1]]$stud)


### define logOR
for (i in 1:N.sim)
{   
  logOR[[i]]=seq(from =1/(NT-1), to = 1, by = 1/(NT-1))  ### true logOR across studies 
  OR[[i]]=c(1,exp(logOR[[i]]))
}


#### define odds per treatment, per study arm
for (i in 1:N.sim)
{  
    data1[[i]]$p.ref=rep(runif(NS,0.01,0.10),each=3)
   data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
   for(j in 1:narms){
   data1[[i]]$odds.2[j]=data1[[i]]$odds.ref[j]*OR[[i]][2]
   data1[[i]]$odds.3[j]=data1[[i]]$odds.ref[j]*OR[[i]][3]
    data1[[i]]$p.t2[j]=data1[[i]]$odds.2[j]/(1+data1[[i]]$odds.2[j])
    data1[[i]]$p.t3[j]=data1[[i]]$odds.3[j]/(1+data1[[i]]$odds.3[j])
  }}


#### generate the events
for (i in 1:N.sim){
   for(j in 1:narms){
    data1[[i]]$r1[j]=-100
    data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.ref[j]) *(data1[[i]]$t[j]==1)
  }
}


for (i in 1:N.sim)
{  
  for(j in 1:narms){
    data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.t2[j]) *(data1[[i]]$t[j]==2)
    data1[[i]]$r3[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.t3[j]) *(data1[[i]]$t[j]==3)
  }}

for (i in 1:N.sim){   data1[[i]]$r=data1[[i]]$r1+data1[[i]]$r2+data1[[i]]$r3}
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(4:12)]}

library(netmeta)
for(i in 1:N.sim){
 data2[[i]] =pairwise(treat=t, event=r, n=n, studlab=stud, data=data1[[i]])
}
data1=data2
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(2:6)]}
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(3,5,6,7,8)]}
##################
