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


pi.star=0.12
pess=2

model="ibb"
max.deviation=3
max.ar=0.75
rand.type=c("Coin","Urn","Block")

PP.bound.grid=seq(0.97,0.999,0.0001)

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

library(parallel)

#cl <- makeCluster(40)
registerDoParallel(cores = 40)
res={}
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
  res[[i]]=foreach(icount(10000)) %dopar% funsimt(response.probs = c(null.response.probs,null.response.probs), 
                                                  ns=ns[[i]], max.ar=max.ar, rand.type=rand.type[3], 
                                                  max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
                                                  beta0_prior_mu=0,beta1_prior_mu=0,
                                                  beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                                                  beta0_df=7,beta1_df=7,
                                                  ibetabinomial.post=ibetabinomial.post,logisticmodel = logistic,replicates = T)
}

Type.I.errorgrid=pstopscreening(ns=ns,PP.bound.grid = PP.bound.grid,res = res)
for (i in 1:dim(Type.I.errorgrid)[1]){
  pstop=c(pstop,min(PP.bound.grid[which(Type.I.errorgrid[i,]<=0.05)]))
}

RMSE.of.interest <- res #scenarios generated from hypothesized null (12%) and alternative (37%, 12%) response rates

RMSE.results=t(RMSE(RMSE.of.interest=RMSE.of.interest,nslist=ns,PP.boundliist=pstop))
RMSE.resultsall=cbind(RMSE.resultsall,RMSE.results)

names(pstop)<-c("15Block","30Block","50Block")
pstop

rownames(RMSE.resultsall) = c("null")
colnames(RMSE.resultsall) =c("15PowerBlock","30PowerBlock","50PowerBlock")

RMSE.resultsall
