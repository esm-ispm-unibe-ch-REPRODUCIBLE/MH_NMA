  ###################################################
  #       Master file for simulations               #
  ###################################################

# This is the master file for running the simulations of the paper 
# "A Mantel-Haenszel model for network meta-analysis of rare events"

  
#empty memory  
rm(list=ls())
  
  
# install packages if needed

#install.packages("netmeta")
#install.packages("beepr")


# load packages
library(MASS)
library(parallel)
library(beepr)  
library(netmeta)




# set working directory
setwd("C:/Users/efthimiou/Google Drive/PROJECT/SNF/paper 2-NMA/Simulations/MH_NMA")

# choose scenario (N.scen=1,2...20)
N.scen=12

# generate the data for scenario k
source(paste("Data generation scenario ",N.scen,".R",sep=""))

# fit the models
source("Fitting the models.R")

# print results 
results=matrix(c(rep(0,12)), nrow=4, byrow=T)
colnames(results)=c("mean.bias", "mean.absolute.bias" ,"mean.coverage")
rownames(results)=c("IV.FE", "IV.RE", "MH","NCH")
results[1,]=c(round(mean(biasIV.FE), digits=2),round(mean(abs(biasIV.FE)), digits=2),round(mean(coverageIV.FE), digits=3)*100)
results[2,]=c(round(mean(biasIV.RE), digits=2),round(mean(abs(biasIV.RE)), digits=2),round(mean(coverageIV.RE), digits=3)*100)
results[3,]=c(round(mean(biasMH), digits=2),round(mean(abs(biasMH)), digits=2),round(mean(coverageMH), digits=3)*100)
results[4,]=c(round(mean(biasNCH), digits=2),round(mean(abs(biasNCH)), digits=2),round(mean(coverageNCH), digits=3)*100)
results

# save a file with results
write.csv(results, file = paste("results scenario ",N.scen,".csv",sep=""))



