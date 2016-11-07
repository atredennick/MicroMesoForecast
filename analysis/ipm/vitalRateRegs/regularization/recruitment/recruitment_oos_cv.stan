data{
  int<lower=0> N; // observations
  int<lower=0> npreds;
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> gid_out[npreds]; // group id
  int<lower=0> Y[N]; // observation vector
  int<lower=0> y_holdout[npreds]; // observation vector
  matrix[N,Covs] C; // climate matrix
  matrix[npreds,Covs] C_out; // climate matrix, holdout
  vector[N] parents1; // crowding vector
  vector[N] parents2; // crowding vector
  vector[npreds] parents1_out; // crowding vector
  vector[npreds] parents2_out; // crowding vector
  real<lower=0> tau; // prior variance
}
parameters{
  real a_mu;
  vector[Yrs] a;
  vector[Covs] b2;
  real dd;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
}
transformed parameters{
real mu[N];
vector[N] climEff;
vector[N] trueP1;
vector[N] trueP2;
vector[N] lambda;
//vector[N] q;
climEff <- C*b2;
  for(n in 1:N){
    trueP1[n] <- parents1[n]*u + parents2[n]*(1-u);
    trueP2[n] <- sqrt(trueP1[n]);
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + dd*trueP2[n] + climEff[n]);
    lambda[n] <- trueP1[n]*mu[n];
    //q[n] <- lambda[n]*theta;
  }
}
model{
// Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,10);
  a_mu ~ normal(0,10);
  dd ~ uniform(-10,10);
  sig_a ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  for(g in 1:G)
   gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
  }
  //lambda1 ~ cauchy(0, 1);
  //tau ~ cauchy(0, 1);
  for(j in 1:Covs)
    //b2[j] ~ normal(0, lambda1[j] * tau);
    b2[j] ~ normal(0, tau);

// Likelihood
Y ~ neg_binomial_2(lambda, theta);
}
generated quantities {
  real int_t;
  real muhat[npreds];
  vector[npreds] climpred;
  vector[npreds] trueP1_pred;
  vector[npreds] trueP2_pred;
  vector[npreds] lambda_hat;
  vector[npreds] qpred;
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  climpred <- C_out*b2;
  int_t <- normal_rng(a_mu, sig_a); // draw random year effect
  for(n in 1:npreds){
    trueP1_pred[n] <- parents1_out[n]*u + parents2_out[n]*(1-u);
    trueP2_pred[n] <- sqrt(trueP1_pred[n]);
    muhat[n] <- exp(int_t + gint[gid_out[n]] + dd*trueP2_pred[n] + climpred[n]);
    lambda_hat[n] <- trueP1_pred[n]*muhat[n];
    //qpred[n] <- lambda_hat[n]*theta;
    log_lik[n] <- neg_binomial_2_log(y_holdout[n], lambda_hat[n], theta);
  } 
}