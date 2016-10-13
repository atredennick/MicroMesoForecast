data{
  int<lower=0> N; // observations
  int<lower=0> npreds; // holdout observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  real<lower=0,upper=1> Y[N]; // observation vector
  real<lower=0> sd_clim; // prior sd on climate effects
  matrix[N,Covs] C; // climate matrix
  vector[N] X; // size vector
  matrix[npreds,Covs] C_out;
  vector[npreds] X_out;
  vector[npreds] y_holdout;
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Covs] b2;
  vector[G] gint;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
  real<lower=0> sigmaSq;
}
transformed parameters{
  real mu[N];
  vector[N] climEff;
  real<lower=0> tau;
  tau <- sqrt(sigmaSq);
  climEff <- C*b2;
  for(n in 1:N)
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + climEff[n];
}
model{
  // Priors
  a_mu ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  gint ~ normal(0, sig_G);
  b2 ~ normal(0,sd_clim);
  a ~ normal(a_mu, sig_a);
  b1 ~ normal(b1_mu, sig_b1);
  sigmaSq ~ inv_gamma(1, 1);
  
  //Likelihood
  Y ~ lognormal(mu, tau);
}
generated quantities {
  real intercept;
  real dens_dep;
  vector[npreds] climpred;
  real muhat[npreds]; // prediction vector
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  climpred <- C_out*b2;
  intercept <- normal_rng(a_mu, sig_a);
  dens_dep <- normal_rng(b1_mu, sig_b1);
  for(n in 1:npreds){
    muhat[n] <- intercept + dens_dep*X_out[n] + climpred[n];
    log_lik[n] <- lognormal_log(y_holdout[n], muhat[n], tau);
  }
}