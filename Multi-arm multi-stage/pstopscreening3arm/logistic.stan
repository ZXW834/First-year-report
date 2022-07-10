data {
  // prior beta0 and beta1 follows t distribution with df = nu
  real beta0_prior_mu;
  real beta1_prior_mu;
  real <lower=0> beta0_prior_sigma;
  real <lower=0> beta1_prior_sigma;
  int <lower=0> beta0_nu;
  int <lower=0> beta1_nu;
  
  int <lower=1>K; // number of arms
  int <lower=0>N; // number of trials
  int <lower=0, upper=K> z[N]; // arm applied on trial n
  // for 2 arm: 1 means treatment, 0 means control
  int <lower=0, upper=1> y[N]; // outcome on trial n
  vector <lower=0, upper=K>[K] x;
  // matrix [N,K-1] x;
  // vector <lower=0, upper=K>[K] arm;
}

//transformed data {
  //vector[K] centered_x;
  //real mean_x;
  //mean_x=mean(x);
  //centered_x=x-mean_x;
  //int successes[K] = rep_array(0, K);
  //int trials[K] = rep_array(0, K);
  //for (n in 1:N) {
  //  trials[z[n]] += 1;
  //  successes[z[n]] += y[n];
  //}
//}
parameters{
  real beta0;
  vector[K] beta1;
  // vector[K-1] beta1;
}

//transformed parameters{
//  vector <lower=0,upper=1>[K] tox_prob;
  //inv_logit(beta_0+beta_1*codeddose[k]) represents the toxicity probability.
//  for (k in 1:K){
//  tox_prob[k]=inv_logit(beta0+beta1*x[k]);
//  }
//}

model{
  // Prior for beta_0 and beta_1
  // In this paper, the prior mean of beta_0 is 0, sd = 100.
  // Prior mean of beta_1 is 1.
  beta0 ~ student_t(beta0_nu,beta0_prior_mu,beta0_prior_sigma);
  beta1 ~ student_t(beta1_nu,beta1_prior_mu,beta1_prior_sigma);
  // The outcome of nth trial follows Bernoulli distribution
  // with toxicity probability of the arm used in the nth trial.
  // There will be a loop with size N for each y[n] ~ bernoulli(tox_prob[z[n]]),
  // where n is in 1:N.
  // target+=bernoulli_logit_lpmf( y | beta0 + x*beta1);
 {
  vector[K] tox_prob;
  for (k in 1:K)
  tox_prob[k] = inv_logit(beta0+beta1[k].* x[k]);
  // tox_prob=inv_logit(beta0 + x*beta1);
  y ~ bernoulli(tox_prob[z]);
  // y~bernoulli_logit(beta0 + x*beta1);

 }
}
generated quantities {
  vector[K] tox_prob;
  for (k in 1:K)
    tox_prob[k] = inv_logit(beta0+beta1[k].* x[k]);
}

