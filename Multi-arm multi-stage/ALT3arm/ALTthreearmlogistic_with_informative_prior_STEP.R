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
PP.bound = c(0.9976,
             0.9964,
             0.9953,
             0.9976,
             0.9968,
             0.9952
)

powercalc=function(Alt.Scenario,PP.bound,ns){
  Power={}
  for (i in 1:3){
    bounds= c(rep(PP.bound[i],length(ns[[i]])-1),PP.bound[i])
    pw=mean(sapply(Alt.Scenario[[i]],function(x) max(x[,1] > bounds | x[,1] < 1-bounds) == 1))
    Power = c(Power,pw)
  }
  return(Power)
}


response.probs=c(0.2,0.2,0.4)
K=length(response.probs)
#cl <- makeCluster(40)
registerDoParallel(cores = 40)
Alt.ScenarioLFC1={}
Alt.ScenarioLFC2={}
power={}
Alt.Resultsall={}
RMSE.resultsall={}

Alt.ScenarioSTEP1={}
Alt.ScenarioSTEP2={}
response.probs=c(0.2,0.3,0.4)
K=length(response.probs)
for (i in 1:3){
  Alt.ScenarioSTEP1[[i]] = foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                                  ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[1], 
                                                                  max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                                  beta0_prior_mu=log(0.2/0.8),beta1_prior_mu=0,
                                                                  beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                                  beta0_df=7,beta1_df=7,
                                                                  ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                                  BARmethod="Thall",replicates=F)
}

for (i in 1:3){
  Alt.ScenarioSTEP2[[i]] = foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                                  ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[2], 
                                                                  max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                                  beta0_prior_mu=log(0.2/0.8),beta1_prior_mu=0,
                                                                  beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                                  beta0_df=7,beta1_df=7,
                                                                  ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                                  BARmethod="Thall",replicates=F)
}
save(Alt.ScenarioSTEP1,Alt.ScenarioSTEP2,file = "STEPlogisticinfThall.RData")
