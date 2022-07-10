library(rstan)
library(boot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))
library(doParallel)

simulate.trial=source("simulatetrial.R")
logistic=stan_model("logistic.stan")
ibetabinomialpost=source("ibetabinomialpost.R")
ibetabinomial.post=ibetabinomialpost[[1]]

funsimt=simulate.trial[[1]]

ntrials = 10000
ns = list(seq(25,225,25),seq(45,225,45),seq(75,225,75)) #group size of 15 for a maximum sample size of 150
null.response.probs = 0.20


pi.star=0.5
pess=2

model="tlr"
max.deviation=3
max.ar=0.80
rand.type=c("Coin","Urn")

PP.bound.grid=seq(0.97,0.997,0.0001)

library(parallel)

#cl <- makeCluster(40)
registerDoParallel(cores = 40)
res2={}
res1={}
pstop={}
RMSE.resultsall={}
pstopscreening=function(ns,PP.bound.grid,res){
  Type.I.errorgridlist={}
  for (i in 1:3){
    Type.I.errorgrid={}
    for (z in 1:length(PP.bound.grid)){
      bounds= c(rep(PP.bound.grid[z],length(ns[[i]])-1),PP.bound.grid[z])  
      Type.I.Error = data.frame(rep(0,length(res)))
      Type.I.Error=sapply(res[[i]],function(x) max(x[,1] > bounds | x[,1] < (1-bounds)) == 1)
      Type.I.errorgrid=c(Type.I.errorgrid,mean(Type.I.Error))
    }
    Type.I.errorgridlist=rbind(Type.I.errorgridlist,Type.I.errorgrid)
  }
  return(Type.I.errorgridlist)
}

for (i in 1:3){
  res1[[i]]=foreach(icount(10000)) %dopar% funsimt(response.probs = c(null.response.probs,null.response.probs,null.response.probs), 
                                                  ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[1], 
                                                  max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                  beta0_prior_mu=0,beta1_prior_mu=0,
                                                  beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                  beta0_df=7,beta1_df=7,
                                                  ibetabinomial.post=ibetabinomial.post,logisticmodel = logistic,replicates = F)
}


for (i in 1:3){
  res2[[i]]=foreach(icount(10000)) %dopar% funsimt(response.probs = c(null.response.probs,null.response.probs,null.response.probs), 
                                                  ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[2], 
                                                  max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                  beta0_prior_mu=0,beta1_prior_mu=0,
                                                  beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                  beta0_df=7,beta1_df=7,
                                                  ibetabinomial.post=ibetabinomial.post,logisticmodel = logistic,replicates = F)
}

save(data,file = "nulllogisticvague.RData")
