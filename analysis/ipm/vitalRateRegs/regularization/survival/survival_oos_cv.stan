data{
  // Training Data
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0,upper=1> Y[N]; // observation vector
  matrix[N,Covs] C; // climate matrix
  vector[N] X; // size vector
  matrix[N,2] W; // crowding matrix
  real beta_tau; // prior sdev for climate effects
  
  // Out-of-Sample Data
  int<lower=0> npreds; // number of oos observations
  int<lower=0, upper=1> yhold[npreds]; // oos response vector
  matrix[npreds,Covs] Chold; // oos climate matrix
  vector[npreds] Xhold; // oos size vector
  matrix[npreds,2] Whold; // oos crowding matrix
  int<lower=0> gid_out[npreds]; // group id
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Covs] b2;
  vector[2] w;
  vector[G] gint;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  real mu[N];
  vector[N] climEff;
  vector[N] crowdEff;
  climEff <- C*b2;
  crowdEff <- W*w;
  for(n in 1:N){
    mu[n] <- inv_logit(a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n]);
  }
}
model{
  // Priors
  a_mu ~ normal(0,100);
  w ~ normal(0,100);
  b1_mu ~ normal(0,100);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  gint ~ normal(0, sig_G);
  b2 ~ normal(0, beta_tau);
  a ~ normal(a_mu, sig_a);
  b1 ~ normal(b1_mu, sig_b1);
  
  // Likelihood
  Y ~ binomial(1,mu);
}
generated quantities{
  vector[npreds] climpred;
  real muhat[npreds]; // prediction vector
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  climpred <- Chold*b2;
  real int_t;
  int_t <- normal_rng(a_mu, sig_a); // draw random year effect
  for(n in 1:npreds){
    muhat[n] <- inv_logit(int_t + gint[gid_out[n]] + b1_mu*Xhold[n] + climpred[n]);
    log_lik[n] <- bernoulli_log(yhold[n], muhat[n]);
  }
}