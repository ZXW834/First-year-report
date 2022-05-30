functions{
  //Define the bivariate binary likelihood
  real logit_jointpdf_EffTox(real muE, real muT,
                             real betaT,real betaE1,real betaE2,real phi,
                             int N, int[] eff, int[] tox, int[] z, 
                             real[] coded_dose,real[] coded_dosesqu){
    real log_likelihood;
    log_likelihood=0;
    for (n in 1:N){
      //for each patient in N, calculate the likelihood function and sum the log likelihood
      real efficacy;
      real toxicity;
      real p;
      toxicity=inv_logit(muT+betaT*coded_dose[z[n]]);
      efficacy=inv_logit(muE+betaE1*coded_dose[z[n]]+betaE2*coded_dosesqu[z[n]]);
      
      p=efficacy^eff[n]*(1-efficacy)^(1-eff[n])*
        toxicity^tox[n]*(1-toxicity)^(1-tox[n])+
        (-1)^(eff[n]+tox[n])*((exp(phi)-1)/(exp(phi)+1))*
        efficacy*(1-efficacy)*toxicity*(1-toxicity);
      log_likelihood+=log(p);
    }
    return log_likelihood;
  }
}

data{
  //Fixed mean an sd of hyperparameters.
  //Prior can be calculated based on prior effective sample size (not very clear at this moment)
  real muE_mean;
  real muE_sd;
  
  real muT_mean;
  real muT_sd;
  
  real betaT_mean;
  real betaT_sd;
  
  real betaE1_mean;
  real betaE1_sd;
  
  real betaE2_mean;
  real betaE2_sd;
  
  real phi_mean;
  real phi_sd;
  
  // The value of p is determined by three equally desirable outcome pairs.
  // This should be calculated before doing HMC.
  real p;
  
  // Three equally desirable outcome pairs
  real eff1; //pi_1*=(eff1,0)
  real tox2; //pi_2*=(1,tox2)
  //pi_3*=(eff3,tox3)
  real eff3; 
  real tox3;
  
  //number of real doses
  int<lower=1> K;
  //vector of real doses
  real<lower=0> real_dose[K]; // Eg. 10mg/l,20mg/l,30mg/l
  //number of patients
  int<lower=0> N;
  //efficacy outcome of n patients
  int<lower=0, upper=1> eff[N];
  //toxicity outcome of n patients
  int<lower=0, upper=1> tox[N];
  //treatment used on nth patient
  int<lower=0, upper=K> z[N];
}

transformed data{
  //Calculate the standardized doses
  real mean_log_dose;  
  real coded_dose[K];
  real coded_dosesqu[K];
  mean_log_dose=0.0;
  for (k in 1:K)
    mean_log_dose=mean_log_dose+log(real_dose[k]);
  mean_log_dose=mean_log_dose/K;

  for (k in 1:K){
  coded_dose[k]=log(real_dose[k])-mean_log_dose;
  coded_dosesqu[k]=coded_dose[k]^2;
  }
}

parameters {
  //hyperparamters
  real muE;
  real muT;
  real betaT;
  real betaE1;
  real betaE2;
  real phi;
}

transformed parameters{
  //toxicity and efficacy probability estimates of each dose
  real<lower=0,upper=1> efficacy[K];
  real<lower=0,upper=1> toxicity[K]; 
  
  //calculate the desirability for decision making
  real desirability[K];
  for (k in 1:K){
  
  //distance_to_most_desirable_point (1,0) which is p(efficacy)=1, p(toxicity)=0
  real distance_to_most_desirable_point;
  
  toxicity[k]=inv_logit(muT+betaT*coded_dose[k]);
  efficacy[k]=inv_logit(muE+betaE1*coded_dose[k]+betaE2*coded_dosesqu[k]);
  
  distance_to_most_desirable_point=
  (((1-efficacy[k])/(1-eff1))^p+((toxicity[k]-0)/(tox2-0))^p)^(1/p);
  
  //calculate desirability function
  desirability[k]=1-distance_to_most_desirable_point;
  }
}

model {
  // prior
  muE~normal(muE_mean,muE_sd);
  muT~normal(muT_mean,muT_sd);
  betaT~normal(betaT_mean,betaT_sd);
  betaE1~normal(betaE1_mean,betaE1_sd);
  betaE2~normal(betaE2_mean,betaE2_sd);
  phi~normal(phi_mean,phi_sd);
  
  //likelihood
  target += logit_jointpdf_EffTox(muE,muT,betaT,betaE1,betaE2,phi,
                                  N, eff, tox, z, coded_dose, coded_dosesqu);
}
