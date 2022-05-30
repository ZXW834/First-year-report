library(rstan)
library(boot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))
library(doParallel)

EffToxdesign=stan_model("Eff-tox.stan")

pval=1.1087016
input=c(0.4,0.75,0.6,0.3)
##True toxicity and efficacy
True_Tox<-c(0.2,0.25,0.35,0.6)
True_Eff<-c(0.25,0.35,0.55,0.6)

# True_Tox<-c(0.1,0.15,0.2,0.25)
# True_Eff<-c(0.3,0.45,0.6,0.8)
# 
# True_Tox<-c(0.3,0.55,0.7,0.8)
# True_Eff<-c(0.6,0.65,0.7,0.75)
# 
# True_Tox<-c(0.65,0.7,0.75,0.8)
# True_Eff<-c(0.5,0.6,0.65,0.7)


real_dose<-c(20,40,60,80)
##Number of real doses
K=length(real_dose)

##Number of patients
MAX_N=60
cohort=3

#Admissibility criteria
lowerpi_E<-0.4
upperpi_T<-0.45
cutoffE=0.1
cutoffT=0.1


#Prior
#The target ESS in the book is 0.9
#Thus the priors are

muT_mean=-7.9593
muT_sd=3.5487

muE_mean=0.7367
muE_sd=2.5423

betaT_mean=1.5482
betaT_sd=3.5018

betaE1_mean=3.4181
betaE1_sd=2.4406

betaE2_mean=0
betaE2_sd=0.2

phi_mean=0
phi_sd=1

finalassign={}

efftox=function(iii,muT_mean=-7.9593,muT_sd=3.5487,muE_mean=0.7367,muE_sd=2.5423,
                betaT_mean=1.5482,betaT_sd=3.5018,betaE1_mean=3.4181,betaE1_sd=2.4406,betaE2_mean=0,betaE2_sd=0.2,
                phi_mean=0,phi_sd=1,lowerpi_E=0.4,upperpi_T=0.45,cutoffE=0.1,cutoffT=0.1,MAX_N=60,cohort=3,
                real_dose=c(20,40,60,80),K=4,True_Tox=c(0.2,0.25,0.35,0.6), True_Eff=c(0.25,0.35,0.55,0.6),
                method=EffToxdesign,pval=1.1087016,input=c(0.4,0.75,0.6,0.3)){
  
  eff <- array(0.0, 0)
  tox<- array(0.0, 0)
  z <- array(0.0, 0)
  
  ##Generating dataset with True prob
  
  datasetE=matrix(rep(0),ncol = K+1,nrow = MAX_N) # Initializing the obs datasetE
  datasetT=matrix(rep(0),ncol = K+1,nrow = MAX_N) # Initializing the obs datasetT
  
  #Generating random data
  for (q in 1:(MAX_N/cohort)){
    for (k in 1:K){
      datasetE[(cohort*q-cohort+1):(cohort*q),k]=rbinom(cohort,1,True_Eff[k])
      datasetT[(cohort*q-cohort+1):(cohort*q),k]=rbinom(cohort,1,True_Tox[k])
    }
  }
  #label the cohort number
  datasetE[,K+1]=rep((1:(MAX_N/cohort)),each=cohort)
  datasetT[,K+1]=rep((1:(MAX_N/cohort)),each=cohort)
  
  start_dose=1
  greedyassign=start_dose
  
  TSselection=start_dose
  criteriaresult=0
  
  for (n in 1:(MAX_N/cohort)){
    
    if (n==1){
      z[(cohort*n-cohort+1):(cohort*n)]=rep(1,cohort)
      #-------------------Apply arm selected-----------------
      eff[(cohort*n-cohort+1):(cohort*n)]=datasetE[(cohort*n-cohort+1):(cohort*n),z[(cohort*n)]]
      tox[(cohort*n-cohort+1):(cohort*n)]=datasetT[(cohort*n-cohort+1):(cohort*n),z[(cohort*n)]]
    }
    else if (sum(criteriaresult)>0){
      greedyassign=which.max(ifelse(criteriaresult,desirabilitymean,NA))
      z[(cohort*n-cohort+1):(cohort*n)]=rep(greedyassign,cohort)
      #-------------------Apply arm selected-----------------
      eff[(cohort*n-cohort+1):(cohort*n)]=datasetE[(cohort*n-cohort+1):(cohort*n),z[(cohort*n)]]
      tox[(cohort*n-cohort+1):(cohort*n)]=datasetT[(cohort*n-cohort+1):(cohort*n),z[(cohort*n)]]
    } 
    else{
      greedyassign=NA
      z[(cohort*n-cohort+1):(cohort*n)]=rep(greedyassign,cohort)
      break()
    }
    
    
    ##-----------------------------Compute posterior distribution based on current data
    
    data=list(K = K, N = (cohort*n),
              eff = array(eff, dim = (cohort*n)), tox = array(tox, dim = (cohort*n)), z = array(z, dim = (cohort*n)),
              real_dose=real_dose,
              #---------------------P and userm defined contour from shiny app efftox contour-----------
              p=pval,
              eff1=input[1],tox2=input[2],eff3=input[3],tox3=input[4],
              #---------------------prior----------------------
              muT_mean=muT_mean, muT_sd=muT_sd,muE_mean=muE_mean,muE_sd=muE_sd,
              betaT_mean=betaT_mean,betaT_sd=betaT_sd,
              betaE1_mean=betaE1_mean,betaE1_sd=betaE1_sd,
              betaE2_mean=betaE2_mean,betaE2_sd=betaE2_sd,
              phi_mean=phi_mean,phi_sd=phi_sd)
    
    # data=list(K = K, N = (cohort*n-cohort),
    #          eff = array(eff, dim = (cohort*n-cohort)), tox = array(tox, dim = (cohort*n-cohort)), z = array(z, dim = (cohort*n-cohort)),
    #          real_dose=real_dose,
    # #---------------------P and userm defined contour from shiny app efftox contour-----------
    #          p=pval,
    #          eff1=input[1],tox2=input[2],eff3=input[3],tox3=input[4],
    # #---------------------prior----------------------
    #          muT_mean=muT_mean, muT_sd=muT_sd,muE_mean=muE_mean,muE_sd=muE_sd,
    #          betaT_mean=betaT_mean,betaT_sd=betaT_sd,
    #          betaE1_mean=betaE1_mean,betaE1_sd=betaE1_sd,
    #          betaE2_mean=betaE2_mean,betaE2_sd=betaE2_sd,
    #          phi_mean=phi_mean,phi_sd=phi_sd)
    
    fit <- rstan::sampling(method, data = data, chains = 1, refresh=0,
                           warmup=2500, iter=5000)
    
    toxicitysamples=rstan::extract(fit,'toxicity')[[1]]
    efficacysamples=rstan::extract(fit,'efficacy')[[1]]
    desirabilitysamples=rstan::extract(fit,'desirability')[[1]]
    desirabilitymean=colMeans(desirabilitysamples)
    #Constructing admissibility set
    
    mindose=min(data$z)
    maxdose=max(data$z)
    
    criteriaresult=((colMeans(efficacysamples-lowerpi_E>0)>cutoffE) & 
                      (colMeans(toxicitysamples-upperpi_T<0)>cutoffT) &
                      sapply(1:K, function(nextdose) nextdose>=(mindose-1) & nextdose<=(maxdose+1)))
    
    toxicitymean=colMeans(toxicitysamples)
    efficacymean=colMeans(efficacysamples)
    
    
    #-------------------Normal Eff-tox design selection-----------------------------
    # if (sum(criteriaresult)>0){
    #   greedyassign=which.max(ifelse(criteriaresult,desirabilitymean,NA))
    #   z[(cohort*n-cohort+1):(cohort*n)]=rep(greedyassign,cohort)
    # } else{
    #   greedyassign=NA
    #   z[(cohort*n-cohort+1):(cohort*n)]=rep(greedyassign,cohort)
    #   break()
    # }
    
    #-------------------Thompson sampling selection-----------------------
    # TSprob=table(factor(max.col(admissible_desirabilitysamples),
    #                     levels=1:ncol(admissible_desirabilitysamples)))/dim(admissible_desirabilitysamples)[1]
    # 
    # TSselection=sample(K,1,replace=T,TSprob)
    # z[n]=TSselection
    
    # #-------------------Apply arm selected-----------------
    # eff[(cohort*n-cohort+1):(cohort*n)]=datasetE[(cohort*n-cohort+1):(cohort*n),z[(cohort*n)]]
    # tox[(cohort*n-cohort+1):(cohort*n)]=datasetT[(cohort*n-cohort+1):(cohort*n),z[(cohort*n)]]
  }
  return(list(greedyassign,c(sum(z==1),sum(z==2),sum(z==3),sum(z==4))))
}


#cl <- makeCluster(40)

registerDoParallel(cores = 30)
res=foreach(icount(1000)) %dopar% efftox()
#res=clusterApply(cl = cl, x=x,fun = efftox,muT_mean=-7.9593,muT_sd=3.5487,muE_mean=0.7367,muE_sd=2.5423,
#                betaT_mean=1.5482,betaT_sd=3.5018,betaE1_mean=3.4181,betaE1_sd=2.4406,betaE2_mean=0,betaE2_sd=0.2,
 #               phi_mean=0,phi_sd=1,lowerpi_E=0.4,upperpi_T=0.45,cutoffE=0.1,cutoffT=0.1,MAX_N=60,cohort=3,
#                real_dose=c(20,40,60,80),K=4,True_Tox=c(0.2,0.25,0.35,0.6), True_Eff=c(0.25,0.35,0.55,0.6),
  #              EffToxdesign=EffToxdesign,pval=1.1087016,input=c(0.4,0.75,0.6,0.3))

#registerDoParallel(cores = 40)
#res={}

#res=foreach(icount(100)) %dopar% efftox()

resultassign={}
patient={}
for (i in 1:length(res)){
  resultassign=c(resultassign,res[[i]][[1]])
  patient=rbind(patient,res[[i]][[2]])
}
resultassign[is.na(resultassign)]="None"
patient[is.na(patient)]=0

out=matrix(rep(NA,10),ncol = 5)
colnames(out)=c("1","2","3","4","None")
out[2,1:4]=colMeans(patient)
for(j in 1:4){
  out[1,j]=sum(resultassign==j)/100*100
}
out[1,5]=sum(resultassign=="None")/100*100
out



