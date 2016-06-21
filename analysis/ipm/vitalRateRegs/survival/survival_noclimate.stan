data{
  // Training Data
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0,upper=1> Y[N]; // observation vector
  vector[N] X; // size vector
  matrix[N,2] W; // crowding matrix
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[2] w;
  vector[G] gint;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  real mu[N];
  vector[N] crowdEff;
  crowdEff <- W*w;
  for(n in 1:N){
    mu[n] <- inv_logit(a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n]);
  }
}
model{
  // Priors
  a_mu ~ normal(0,10);
  w ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  gint ~ normal(0, sig_G);
  a ~ normal(a_mu, sig_a);
  b1 ~ normal(b1_mu, sig_b1);
  
  // Likelihood
  Y ~ binomial(1,mu);
}