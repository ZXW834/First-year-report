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
ns = seq(30,150,30) #group size of 15 for a maximum sample size of 150
null.response.probs = 0.12



max.deviation=3
max.ar=0.75
rand.type=c("Coin","Urn","Block")

PP.bound=c(0.9863,0.9928,0.9964,0.9925, 0.9856,0.9921,0.9966,0.9928, 0.9847,0.9921,0.9961,0.9928)


powercalc=function(Alt.Scenario,PP.bound,ns){
  Power={}
  for (i in 1:3){
    bounds= c(rep(PP.bound[i],length(ns[[i]])-1),PP.bound[i])
    pw=mean(sapply(Alt.Scenario[[i]],function(x) max(x[,1] > bounds | x[,1] < 1-bounds) == 1))
    Power = c(Power,pw)
  }
  return(Power)
}

ALT=function(Alt.Scenario,PP.bound,ns){
  resultall={}
  for (i in 1:3){
    alt.foo.15 = sapply(Alt.Scenario[[i]],function(x){
      bounds= c(rep(PP.bound[i],length(ns[[i]])-1),PP.bound[i])
      stop.at = which.max(x[,1:(K-1)] > c(bounds[-length(ns)],0)| x[,1:(K-1)] < (1-bounds)) 
      #determine when trial exceeds posterior probability stopping boundary for the first time
      n = sum(x[stop.at,seq(K,(3*K-2),2)]) #total number of people in trial before stopping, i.e., nET+nC
      narm=x[stop.at,seq(K,(3*K-2),2)] #total number of people on each arm including control before stopping
      
      nE=narm[2:K]
      nC=narm[1]
      nE.minus.nC = nE - nC #treatment arm sample size difference
      nE.less.nC = nE<nC
      result=c(n,nE.minus.nC,nE.less.nC)
      return(result)
    })
    resultall=rbind(resultall,alt.foo.15)
  }
  
  return(resultall) #returns desired parameters: n, nE-nC, and if nE-nC < 0
}

# function to generate Squared Error
SE = function(x,ns=ns,bounds=bounds){
  stop.at = which.max(x[,1] > c(bounds[-length(ns)],0)| x[,1] < (1-bounds)) #determine when trial exceeded PP stopping boundary
  if (stop.at == 1){rand.probs = 0.5} #a randomization probability of 0.5 was used to generate treatment assignments for first sequential group
  else {
    rand.probs = c(0.5,x[1:(stop.at-1),1]) #determine randomization probabilities used to generate treatment assignments in each sequential group
    for (b in 1:length(rand.probs)){ #max and min randomization probabilities are 0.75 and 0.25
      if (rand.probs[b] < 0.25){rand.probs[b] = 0.25}
      if (rand.probs[b] > 0.75){rand.probs[b] = 0.75}}}
  expected.nE = sum(rand.probs*(150/length(ns))) #determine expected number of subjects assigned to treatment group
  actual.nE = x[stop.at,4] #determine actual number of subjects assinged to treatment group
  nE.Error = expected.nE - actual.nE #calculate difference between expected and actual number of treatment assignments
  SE = (nE.Error)^2 #compute Squared Error
  return(SE)
}

RMSE=function(RMSE.of.interest,nslist=ns,PP.boundliist=PP.bound){
  RMSE.results = data.frame(rep(0,3)) #create object to store RMSE results
  for(i in 1:3) {
    PP.bound=PP.boundliist[i]
    ns=nslist[[i]]
    bounds= c(rep(PP.bound,length(ns)-1),PP.bound)
    RMSE.results[i,] = sqrt(mean(sapply(RMSE.of.interest[[i]], function(x) SE(x,ns=ns,bounds=bounds))))
  }
  return(RMSE.results)
}

EminusC = function(Alt.Scenario,nslist=ns,bounds){
  res={}
  PP.bound=bounds
  ns=nslist
  PP.bound= c(rep(PP.bound,length(ns)-1),PP.bound) 
  x=Alt.Scenario
  neminusnc = sapply(x,function(x){
    stop.at=which.max(x[,1] > c(PP.bound[-length(ns)],0)| x[,1] < (1-PP.bound))
    nE=x[stop.at,4] 
    nC=x[stop.at,2]  
    rest=nE-nC
    return(rest)
  }) #determine when trial exceeded PP stopping boundary
  res=neminusnc
  return(res)
}
sumtest=matrix(rep(0,10000*12),ncol = 12)

#cl <- makeCluster(40)
registerDoParallel(cores = 5)
Alt.Scenario={}
power={}
Alt.Resultsall={}
RMSE.resultsall={}
response.probs=c(0.12,0.37)
K=length(response.probs)

model="ibb"
pi.star=0.5
pess=2
Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                       ns=ns, max.ar=max.ar, rand.type=rand.type[1], 
                                                       max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                       beta0_prior_mu=0,beta1_prior_mu=0,
                                                       beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                       beta0_df=7,beta1_df=7,
                                                       ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                       BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[1])
sumtest[,1]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                    ns=ns, max.ar=max.ar, rand.type=rand.type[2], 
                                                    max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                    beta0_prior_mu=0,beta1_prior_mu=0,
                                                    beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                    beta0_df=7,beta1_df=7,
                                                    ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                    BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[5])
sumtest[,5]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                    ns=ns, max.ar=max.ar, rand.type=rand.type[3], 
                                                    max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                    beta0_prior_mu=0,beta1_prior_mu=0,
                                                    beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                    beta0_df=7,beta1_df=7,
                                                    ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                    BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[9])
sumtest[,9]=test
#------------------------------------------------------------------------------
model="ibb"
pi.star=0.12
pess=2
Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                    ns=ns, max.ar=max.ar, rand.type=rand.type[1], 
                                                    max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                    beta0_prior_mu=0,beta1_prior_mu=0,
                                                    beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                    beta0_df=7,beta1_df=7,
                                                    ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                    BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[3])
sumtest[,3]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                    ns=ns, max.ar=max.ar, rand.type=rand.type[2], 
                                                    max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                    beta0_prior_mu=0,beta1_prior_mu=0,
                                                    beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                    beta0_df=7,beta1_df=7,
                                                    ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                    BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[7])
sumtest[,7]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                    ns=ns, max.ar=max.ar, rand.type=rand.type[3], 
                                                    max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                    beta0_prior_mu=0,beta1_prior_mu=0,
                                                    beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                    beta0_df=7,beta1_df=7,
                                                    ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                    BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[11])
sumtest[,11]=test

#------------------------------------------------------------------------------
model="tlr"
Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                     ns=ns, max.ar=max.ar, rand.type=rand.type[1], 
                                                     max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                     beta0_prior_mu=0,beta1_prior_mu=0,
                                                     beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                     beta0_df=7,beta1_df=7,
                                                     ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                     BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[2])
sumtest[,2]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                     ns=ns, max.ar=max.ar, rand.type=rand.type[2], 
                                                     max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                     beta0_prior_mu=0,beta1_prior_mu=0,
                                                     beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                     beta0_df=7,beta1_df=7,
                                                     ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                     BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[6])
sumtest[,6]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                     ns=ns, max.ar=max.ar, rand.type=rand.type[3], 
                                                     max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                     beta0_prior_mu=0,beta1_prior_mu=0,
                                                     beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                     beta0_df=7,beta1_df=7,
                                                     ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                     BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[10])
sumtest[,10]=test

#------------------------------------------------------------------------------
model="tlr"
Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                     ns=ns, max.ar=max.ar, rand.type=rand.type[1], 
                                                     max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                     beta0_prior_mu=log(0.12/0.88),beta1_prior_mu=0,
                                                     beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                     beta0_df=7,beta1_df=7,
                                                     ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                     BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[4])
sumtest[,4]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                     ns=ns, max.ar=max.ar, rand.type=rand.type[2], 
                                                     max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                     beta0_prior_mu=log(0.12/0.88),beta1_prior_mu=0,
                                                     beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                     beta0_df=7,beta1_df=7,
                                                     ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                     BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[8])
sumtest[,8]=test

Alt.Scenario= foreach(icount(10000)) %dopar% funsimt(response.probs = response.probs, 
                                                     ns=ns, max.ar=max.ar, rand.type=rand.type[3], 
                                                     max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                     beta0_prior_mu=log(0.12/0.88),beta1_prior_mu=0,
                                                     beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                     beta0_df=7,beta1_df=7,
                                                     ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                     BARmethod="Thall",replicates=T)

test=EminusC(Alt.Scenario,nslist=ns,bounds=PP.bound[12])
sumtest[,12]=test

colnames(sumtest)=c("Beta binomial 0.5","Logistic 0.5","Beta binomial 0.12","Logistic 0.12",
                    "Beta binomial 0.5","Logistic 0.5","Beta binomial 0.12","Logistic 0.12",
                    "Beta binomial 0.5","Logistic 0.5","Beta binomial 0.12","Logistic 0.12")

library(reshape)
library(ggplot2)
# PNG device
png("30Pcoin.png",width = 960)

# Code

test1=melt(sumtest[,1:4],value.name="cohort")
colnames(test1)=c("trial","cohortsize","gap")
ggplot(data=test1,aes(gap,color=cohortsize)) + stat_ecdf(geom = "step")

# Close device
dev.off()

# PNG device
png("30Pmass.png",width = 960)

# Code
test2=melt(sumtest[,5:8],value.name="cohort")
colnames(test2)=c("trial","cohortsize","gap")
ggplot(data=test2,aes(gap,color=cohortsize)) + stat_ecdf(geom = "step")

# Close device
dev.off()

# PNG device
png("30Pblock.png",width = 960)

# Code
test3=melt(sumtest[,9:12],value.name="cohort")
colnames(test3)=c("trial","cohortsize","gap")
ggplot(data=test3,aes(gap,color=cohortsize)) + stat_ecdf(geom = "step")

# Close device
dev.off()


