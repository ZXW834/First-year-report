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
ns = list(seq(15,150,15),seq(30,150,30),seq(50,150,50)) #group size of 15 for a maximum sample size of 150
null.response.probs = 0.12


pi.star=0.5
pess=2

model="tlr"
max.deviation=3
max.ar=0.75
rand.type=c("Coin","Urn","Block")
PP.bound = c( 0.9941,0.9925,0.9904,0.9940,0.9928,0.9908,0.9964,0.9952,0.9937  )

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

K=length(response.probs)
#cl <- makeCluster(40)
registerDoParallel(cores = 5)
Alt.Scenario={}
power={}
Alt.Resultsall={}
RMSE.resultsall={}
response.probs=c(0.12,0.37)

for (i in 1:3){
Alt.Scenario[[i]] = foreach(icount(1000)) %dopar% funsimt(response.probs = response.probs, 
                                                       ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[1], 
                                                       max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                       beta0_prior_mu=log(0.12/0.88),beta1_prior_mu=0,
                                                       beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                       beta0_df=7,beta1_df=7,
                                                       ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                       BARmethod="Thall",replicates=T)
}

powerresult=powercalc(Alt.Scenario=Alt.Scenario,PP.bound = PP.bound[1:3],ns=ns)
power=c(power,powerresult)

alt.foo.15=ALT(Alt.Scenario=Alt.Scenario,PP.bound=PP.bound[1:3],ns=ns)

Alt.Results <- t(data.frame(rowMeans(alt.foo.15))); #returns average n, nE-nC, and proportion of trials in which nE-nC<0
Alt.Resultsall=c(Alt.Resultsall,Alt.Results)

#compute RMSE of interest
RMSE.of.interest <- Alt.Scenario #scenarios generated from hypothesized null (12%) and alternative (37%, 12%) response rates

RMSE.results=t(RMSE(RMSE.of.interest=RMSE.of.interest,nslist=ns,PP.boundliist=PP.bound[1:3]))
RMSE.resultsall=cbind(RMSE.resultsall,RMSE.results)

for (i in 1:3){
  Alt.Scenario[[i]] = foreach(icount(1000)) %dopar% funsimt(response.probs = response.probs, 
                                                            ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[2], 
                                                            max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                            beta0_prior_mu=log(0.12/0.88),beta1_prior_mu=0,
                                                            beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                            beta0_df=7,beta1_df=7,
                                                            ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                            BARmethod="Thall",replicates=T)
}

powerresult=powercalc(Alt.Scenario=Alt.Scenario,PP.bound = PP.bound[4:6],ns=ns)
power=c(power,powerresult)

alt.foo.15=ALT(Alt.Scenario=Alt.Scenario,PP.bound=PP.bound[4:6],ns=ns)

Alt.Results <- t(data.frame(rowMeans(alt.foo.15))); #returns average n, nE-nC, and proportion of trials in which nE-nC<0
Alt.Resultsall=c(Alt.Resultsall,Alt.Results)

#compute RMSE of interest
RMSE.of.interest <- Alt.Scenario #scenarios generated from hypothesized null (12%) and alternative (37%, 12%) response rates

RMSE.results=t(RMSE(RMSE.of.interest=RMSE.of.interest,nslist=ns,PP.boundliist=PP.bound[4:6]))
RMSE.resultsall=cbind(RMSE.resultsall,RMSE.results)

for (i in 1:3){
  Alt.Scenario[[i]] = foreach(icount(1000)) %dopar% funsimt(response.probs = response.probs, 
                                                            ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[3], 
                                                            max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                            beta0_prior_mu=log(0.12/0.88),beta1_prior_mu=0,
                                                            beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                            beta0_df=7,beta1_df=7,
                                                            ibetabinomial.post=ibetabinomial.post,logisticmodel=logistic,
                                                            BARmethod="Thall",replicates=T)
}

powerresult=powercalc(Alt.Scenario=Alt.Scenario,PP.bound = PP.bound[7:9],ns=ns)
power=c(power,powerresult)

alt.foo.15=ALT(Alt.Scenario=Alt.Scenario,PP.bound=PP.bound[7:9],ns=ns)

Alt.Results <- t(data.frame(rowMeans(alt.foo.15))); #returns average n, nE-nC, and proportion of trials in which nE-nC<0
Alt.Resultsall=c(Alt.Resultsall,Alt.Results)

#compute RMSE of interest
RMSE.of.interest <- Alt.Scenario #scenarios generated from hypothesized null (12%) and alternative (37%, 12%) response rates

RMSE.results=t(RMSE(RMSE.of.interest=RMSE.of.interest,nslist=ns,PP.boundliist=PP.bound[7:9]))
RMSE.resultsall=cbind(RMSE.resultsall,RMSE.results)
#-------------------------Summary results-----------------------------------

names(power) = c("15PowerCoin","30PowerCoin","50PowerCoin",
                 "15PowerUrn","30PowerUrn","50PowerUrn",
                 "15PowerBlock","30PowerBlock","50PowerBlock")
power

names(Alt.Resultsall) <- c("Average Sample Size15Coin", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size30Coin", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size50Coin", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size15Urn", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size30Urn", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size50Urn", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size15Block", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size30Block", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"),
                           "Average Sample Size50Block", paste0("Average nE",1:(K-1),"-nC"),paste0("mean(nE",1:(K-1),"-nC<0)"))

Alt.Resultsall



rownames(RMSE.resultsall) = c("alternative")
colnames(RMSE.resultsall) =c("15PowerCoin","30PowerCoin","50PowerCoin",
                             "15PowerUrn","30PowerUrn","50PowerUrn",
                             "15PowerBlock","30PowerBlock","50PowerBlock")

RMSE.resultsall
