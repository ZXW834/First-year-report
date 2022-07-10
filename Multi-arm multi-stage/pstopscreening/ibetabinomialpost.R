ibetabinomial.post=function(n,y,pi.star=0.5, pess=2){
  
  #First element of n and y are from control
  K=length(n)
  # Posterior Probability of each arm better than the same control arm
  # See https://www.evanmiller.org/bayesian-ab-testing.html#cite1
  # Treatment_k~beta(a_k,b_k), Control~beta(a_1,b_1)
  
  # p.prior*ess.prior: prior success
  # (1-p.prior)*ess.prior: prior failure
  
  
  post.prob={}
  for (k in 2:K){
    
    post.prob[k-1] <- unlist(integrate(function(x)
      pbeta(x,y[k]+pi.star*pess,(n[k]-y[k])+(1-pi.star)*pess,lower.tail=FALSE)*dbeta(x,y[1]+pi.star*pess,(n[1]-y[1])+(1-pi.star)*pess),
      lower=0,upper=1))$value
    names(post.prob[k-1])= paste("Treatment",k-1,"vs", "Control", sep = " ")
  }
  
  #Slower version
  # # p.prior*ess.prior: prior success
  # # (1-p.prior)*ess.prior: prior failure
  # narm=length(n)
  # rn=matrix(rbeta(random.number*narm,
  #                 y+p.prior*ess.prior,
  #                 n-y+(1-p.prior)*ess.prior),
  #           random.number, byrow = TRUE)
  # 
  # controlrn=rn[,1]
  # treatmentrn=as.matrix(rn[,-1])
  # 
  # # rnT=rbeta(random.number,y[1]+p.prior*ess.prior,n[1]-y[1]+(1-p.prior)*ess.prior)
  # # rnC=rbeta(random.number,y[2]+p.prior*ess.prior,n[2]-y[2]+(1-p.prior)*ess.prior)
  # # postprob<-mean(rnT>rnC)
  # 
  # post.prob=colMeans(treatmentrn>controlrn)
  
  return(post.prob)
}