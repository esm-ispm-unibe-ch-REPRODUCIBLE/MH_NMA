set.seed(42)
N.sim=1000


###### Scenario 1 ######
# Fixed effects
# 5 treatments, 1 study per comparison

data1=list()
NT=5 #### number of treatments in the network
NS=2 #### number of studies per comparison
tau=0 ####  heterogeneity SD
p.ref=runif(N.sim,0.01,0.05) ### probability of an event, reference arm
logOR=sort(runif(NT-1,0,1))  ### true logOR across studies 
OR=c(1,exp(logOR))
Npmin=20 #### minimum number of patients per arm
Npmax=40 #### maximum number of patients per arm

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
Np=c(0,rep(N.stud))
for (j in 1:N.stud){Np[j]=round(runif(1,Npmin,Npmax))}

#### define probabilities per treatment, per study arm
for (i in 1:N.sim)
{  
odds=c(p.ref[i]/(1-p.ref[i]),rep(0,NT-1))
for (i in 2:(NT)){odds[i]=odds[1]*OR[i]}
p.event=odds/(odds+1)
}
#### generate the data
for (i in 1:N.sim)
{  
data1[[i]]=data.frame(t1,t2)
data1[[i]]$n1=data1[[i]]$n2=Np
data1[[i]]$studlab=c(1:(N.stud))
data1[[i]]$t1=t1
data1[[i]]$t2=t2
for(j in 1:N.stud){
 data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n1[j],p.event[data1[[i]]$t1][j]) 
 data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n2[j],p.event[data1[[i]]$t2][j]) 
}}
##################
