set.seed(42) # the answer to the ultimate question of life, the universe and everything
N.sim=1000
library(MASS)



###### Scenario 1 ######
# FIXED effects
name="Scenario-1.RData"
data1=list()
logOR=list()
OR=list()
NT=5 #### number of treatments in the network
NS=2 #### number of studies per comparison
tau=0 ####  heterogeneity SD
Npmin=30 #### minimum number of patients per arm
Npmax=60 #### maximum number of patients per arm

### define treatment indices
t1=c()
t2=c()
for (i in 1:(NT-1)){
  for (k in (i+1):NT){
    for(j in 1:NS){
      t1=c(t1,i)
      t2=c(t2,k)      }}}
N.stud=length(t1)

### define patients per treatment arm
for (i in 1:N.sim)
{   
  logOR[[i]]=seq(from =1/(NT-1), to = 1, by = 1/(NT-1))  ### true logOR across studies 
  OR[[i]]=c(1,exp(logOR[[i]]))
  
  
  data1[[i]]=data.frame(t1,t2)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$n1=data1[[i]]$n2=round(runif(N.stud,Npmin,Npmax))}


#### define probabilities per treatment, per study arm
for (i in 1:N.sim)
{  
  data1[[i]]$p.ref=runif(N.stud,0.03,0.05) #### study-specific probability of an event in treatment 1
  data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
}

#### define probabilities per treatment, per study arm
for (i in 1:N.sim)
{  
  for(j in 1:N.stud){
    data1[[i]]$odds.t1[j]=data1[[i]]$odds.ref[j]*OR[[i]][data1[[i]]$t1[j]]
    data1[[i]]$odds.t2[j]=data1[[i]]$odds.ref[j]*OR[[i]][data1[[i]]$t2[j]]
    data1[[i]]$p.t1[j]=data1[[i]]$odds.t1[j]/(1+data1[[i]]$odds.t1[j])
    data1[[i]]$p.t2[j]=data1[[i]]$odds.t2[j]/(1+data1[[i]]$odds.t2[j])
  }}


#### generate the data
for (i in 1:N.sim)
{  
  for(j in 1:N.stud){
    data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n1[j],data1[[i]]$p.t1[j]) 
    data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n2[j],data1[[i]]$p.t2[j]) 
  }}
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(6:11)]}
##################

