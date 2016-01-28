data{
  int<lower=0> N; // observations
  int<lower=0> npreds;
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  real<lower=0> tau_beta; // prior sdev for climate effects
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> gid_out[npreds]; // group id holdout
  vector[N] Y; // observation vector
  vector[npreds] y_holdout;
  matrix[N,Covs] C; // climate matrix
  matrix[npreds,Covs] Chold;
  vector[N] X; // size vector
  vector[npreds] Xhold;
  matrix[N,2] W; // crowding matrix
  matrix[npreds,2] Whold;
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Covs] b2;
  vector[2] w;
  real gint[G];
  real tau;
  real tauSize;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  vector[N] mu;
  real<lower=0> sigma[N];
  vector[N] climEff;
  vector[N] crowdEff;
  climEff <- C*b2;
  crowdEff <- W*w;
  for(n in 1:N){
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n];
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  
  }
}
model{
  // Priors
  a_mu ~ normal(0,100);
  w ~ normal(0,10);
  b1_mu ~ normal(0,100);
  tau ~ normal(0,100);
  tauSize ~ normal(0,100);
  sig_a ~ cauchy(0,2);
  sig_b1 ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  b2 ~ normal(0, tau_beta);
  for(g in 1:G)
    gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ normal(mu, sigma);
}
generated quantities {
  vector[npreds] climpred;
  vector[npreds] crowdhat;
  vector[npreds] sigmahat;
  vector[npreds] muhat;
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  climpred <- Chold*b2;
  crowdhat <- Whold*w;
  real int_t;
  int_t <- normal_rng(a_mu, sig_a); // draw random year effect
  for(n in 1:npreds){
    muhat[n] <- int_t + gint[gid_out[n]] + b1_mu*Xhold[n] + crowdhat[n] + climpred[n];
    sigmahat[n] <- sqrt((fmax(tau*exp(tauSize*muhat[n]), 0.0000001))); 
    log_lik[n] <- normal_log(y_holdout[n], muhat[n], sigmahat[n]);
  }
}