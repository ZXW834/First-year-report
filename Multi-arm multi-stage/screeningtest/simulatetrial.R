#Two arms situation simulation function using stan

#response.probs: first element from control; second element from treatment

simulate.trial <- function(ii,response.probs, ns, max.ar=0.75, 
                            rand.type=c("Coin","Urn","Block"), max.deviation=3, 
                            model=c("tlr","ibb"),
                            pi.star=0.5, pess=2, 
                            beta0_prior_mu=0,beta1_prior_mu=0,
                            beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
                            beta0_df=7,beta1_df=7,
                            logisticmodel=logisticmodel,
                            ibetabinomial.post=ibetabinomial.post,
                            BARmethod="Thall",replicates=TRUE){
  
  #Initialize Data 
  K=length(response.probs)
  n = rep(0,K)
  y1 = rep(0,K) 
 
  
  rand.prob = 1/K
  randomprob=rep(rand.prob,K)
  
  rpblocktwoarm=rand.prob
  
  z <- array(0.0, 0)
  y <- array(0.0, 0)
  
  
  #Storage object for posterior probabilities
  stats = matrix(NA,nrow=length(ns),ncol=K-1+K*2)
  #K-1 posterior probability better than control and K columns number of success + K columns number of patient

  #Use name function to assign the names of each column of stats matrix
  namefunc=function(K){
    name1={}
    name2={}
    for (n in 1:(K-1)){
      name1=c(name1,paste0("PP",n,"C"))
      name2=c(name2,paste0("nE",n),paste0("yE",n))
    }
    return(c(name1,"nC","yC",name2))
  }
  
  
  statsname=namefunc(K)
  colnames(stats) = statsname
  rownames(stats) = 1:length(ns)
  
  post.prob.best.mat=matrix(0,length(ns),K)
  
  #Generate Data under Block AR scheme and true response probabilities
  for(group in 1:nrow(stats)){
    #Number of new enrollees during current group
    n.new = c(0,ns)[group+1]-c(0,ns)[group]
    
    #Generate assignments and outcomes for current interval
    if(rand.type == "Coin"){
      
      #Assuming K arms including Control (K-1 treatment vs 1 Control)
      
      #Randomisation to each k arm including control arm, where rand.prob=c(AR1,AR2,AR3,...,ARK-1) for K arm allocation
      #rand.prob was suggested by Trippa et.al (2012; 2014)
      #The first element is the control group, the elements 2 to K are treatment groups
      randomsample=sample(K,n.new,replace = T,prob = randomprob)
      randcount=table(factor(randomsample,levels = seq(1,K)))
      nstage=randcount
      
      #Simulate the response of randomised patients
      ystage=rbinom(K,nstage,response.probs)
      
    }
    if(rand.type == "Urn"){
      allocation = NULL; rand.prob.temp = randomprob; count=matrix(rep(0,K),ncol = K)
      for(i in 1:n.new){
        urnprob=randomprob*max.deviation-as.vector(count)+(i-1)*randomprob
        
        maxcon=ifelse(urnprob>0,urnprob,0)
        rand.prob.temp=maxcon/sum(maxcon)
        # allocation[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
        allocation[i]=sample(K,1,T,prob = rand.prob.temp)
        count[allocation[i]]= count[allocation[i]]+1
      }
      randomsample=allocation
      nstage=as.vector(count)
      ystage=rbinom(K,nstage,response.probs)
    }
    if(rand.type == "Block"){
      u = rbinom(1,1,n.new*rpblocktwoarm - floor(n.new*rpblocktwoarm))
      nE.new = u*ceiling(n.new*rpblocktwoarm)+(1-u)*floor(n.new*rpblocktwoarm)
      nC.new = n.new-nE.new
      yE.new = rbinom(1,nE.new,response.probs[2])
      yC.new = rbinom(1,nC.new,response.probs[1])
      ystage=c(yC.new,yE.new)
      nstage=c(nC.new,nE.new)
      nsequence=c(rep(1,nstage[1]),rep(2,nstage[2]))
      randomsample=sample(nsequence)
    }
    #Update dataset
    
    #Use sample here will make movement in seed sequence, therefore different Type I will be generated compared with example of origin
    #When using seed 1994, the result is slightly different but the code is correct.
    n=n+nstage
    y1=y1+ystage
    stats2=as.vector(matrix(c(n, y1), nrow = 2, byrow = TRUE))
    
    znew=randomsample
    ynew=rep(0,sum(nstage))
    
    m=seq(1,K)
    y.outcome.for.each.k=sapply(m, function(m) {ynew[which(znew==m)[sample(ystage[m])]]<-1; return(ynew)})
    ynew=rowSums(y.outcome.for.each.k)
    
    #These two command will lead to some mistakes when length(which(znew==2)) is 1    
    #For example, which(znew==2) is 14, sample(which(znew==2),yE.new) will not give the value of 14 which means that the 14th patient
    #uses treatment. It will give a random number from 1 to 14. Use debug(simulate.trial2) and set.seed(1)
    # re10=lapply(1:10,function(trial) simulate.trial2(response.probs = c(0.95,0.95), ns=ns, max.ar=max.ar, rand.type=rand.type, 
    #                     max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess,
    #                     beta0_prior_mu=0,beta1_prior_mu=0,
    #                     beta0_prior_sigma=2.5,beta1_prior_sigma=2.5,
    #                     beta0_df=7,beta1_df=7))
    #to check the result of 10th iteration (4th cohort).
    #Debugged at 20:55 on 27/03/2022 by Ziyan Wang
    
    
    z=c(z,znew)
    y=c(y,ynew)
    N = length(z)
    
    
    #Update posterior probability that CCL is superior to ICU and randomization probability for next group
    
    # When inputing n and y the first element is the data from control group, the rest data are data from treatment arms
    if(model=="ibb"){
      resultibb=ibetabinomial.post(n=n,y=y1,pi.star=pi.star, pess=pess)
      # resultibb2=ibetabinomial.post(n=c(nC,nE),y=c(yC,yE),pi.star=pi.star, pess=pess)
      # More than one posterior probability for more than one treatment arm situation
      if (K==2){
      post.prob.btcontrol = resultibb 
      randomprob = min(max.ar,max(1-max.ar,post.prob.btcontrol))
      rpblocktwoarm=randomprob
      
      post.prob.best=c(1-post.prob.btcontrol,post.prob.btcontrol)
      post.prob.best=c(1-randomprob,randomprob)
      randomprob=post.prob.best
      
      }
      else{
        stop("Error for ibb")
      }
      # samptwoarmtox= resultibb[[2]]
    }
    
    if(model=="tlr") {
      
      data=list(K = K, N = N,
                y = array(y, dim = N), z = array(z, dim = N),x=seq(0,K-1),
                beta0_prior_mu=beta0_prior_mu,beta1_prior_mu=beta1_prior_mu,beta0_prior_sigma=beta0_prior_sigma,
                beta1_prior_sigma=beta1_prior_sigma,beta0_nu=beta0_df,beta1_nu=beta1_df)
      
      fit <- rstan::sampling(logisticmodel,
                             data = data, chains = 1, refresh=0, warmup=2500, iter=5000)
      
      sampeff=rstan::extract(fit,'tox_prob')[[1]]
      control=matrix(sampeff[,1])
      treatment=matrix(sampeff[,-1],ncol = dim(sampeff)[2]-1)
      treatindex=seq(1,dim(treatment)[2])
      #posterior of each treatment better than control
      post.prob.btcontrol=colMeans(sapply(treatindex, function(treatindex) {postprob=treatment[,treatindex]>control; return(postprob)}))
      
      #posterior probability of each arm including control to be the best.
      #This probability is used in Thall's approach (2015)

      for (q in 1: K){
        post.prob.best.mat[group,q]=(sum(max.col(sampeff)==q))/2500
      }
      post.prob.best=post.prob.best.mat[group,]
      

      
      # samp=rstan::extract(fit,'beta1')[[1]]
      # samptox=rstan::extract(fit,'tox_prob')[[1]]
      # post.prob=mean(samp>0)
      
      if (BARmethod=="Trippa"){
#---------------------Trippa's approach---------------------     
    ##Tuning the paprameter using method mentioned in Trippa's paper (2014)
    gamma_stage=10*((group/dim(stats)[1]))^0.75
    eta_stage=0.25*(group/dim(stats)[1])
    ##Reweigh the allocation probability
    ###K>=1, treatment group
    allocate_trt=post.prob.btcontrol^gamma_stage/sum(post.prob.btcontrol^gamma_stage)
    ###k=0, control group
    allocate_control=1/K*(exp(max(n[-1])-n[1]))^eta_stage
    
    sum_pi=allocate_control+sum(allocate_trt)
    alloc.prob.btcontrol=c(allocate_control/sum_pi,allocate_trt/sum_pi)
    
    randomprob=alloc.prob.btcontrol
    }
#----------------------------------------------------------
    
    if (BARmethod=="Thall"){
#---------------------Thall's approach---------------------    
    ##Tuning parameter c for Thall's approach
    
    if (replicates==T){
     c=1
    }
    else{
     c=group/(2*dim(stats)[1]) 
    }
    
    ##Reweigh the allocation probability
    alloc.prob.best=post.prob.best^c/sum(post.prob.best^c)
    rpblocktwoarm=min(max.ar,max(1-max.ar,post.prob.btcontrol))

    randomprob=alloc.prob.best
    
    #------------------------Allocation bounds restriction (two arm)---------------
    lower=ifelse(alloc.prob.best<(1-max.ar),1-max.ar,alloc.prob.best)
    upper=ifelse(lower>max.ar,max.ar,lower)
    randomprob = upper

    #Power reduced
    }
#-------------------------------------------------------------   
      
    }
    
    # More than one posterior probability for more than one treatment arm situation.
    # Need to adjust for each probability
    
    
    # rand.prob = min(max.ar,max(1-max.ar,post.prob))
    
    #Store Relevant Statistics
    #If K arm situation this part should be revised, See Trippa et.al (2012)
    stats1=post.prob.btcontrol
    stats[group,] = c(stats1,stats2)
    # stats2[group,] = c(resultibb2,nE,yE,nC,yC)
  }
  
  return(stats)
}
