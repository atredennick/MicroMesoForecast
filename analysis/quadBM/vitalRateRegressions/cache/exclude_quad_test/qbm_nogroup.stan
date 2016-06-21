data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  real<lower=0,upper=1> Y[N]; // observation vector
  matrix[N,Covs] C; // climate matrix
  vector[N] X; // size vector
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Covs] b2;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> tau;
}
transformed parameters{
  real mu[N];
  vector[N] climEff;
  climEff <- C*b2;
  for(n in 1:N)
    mu[n] <- a[yid[n]] + b1[yid[n]]*X[n] + climEff[n];
}
model{
  // Priors
  a_mu ~ uniform(-300,300);
  b1_mu ~ uniform(-100,100);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  b2 ~ uniform(-10,10);
  a ~ normal(a_mu, sig_a);
  b1 ~ normal(b1_mu, sig_b1);
  tau ~ cauchy(0,5);

  //Likelihood
  Y ~ lognormal(mu, tau);
}
