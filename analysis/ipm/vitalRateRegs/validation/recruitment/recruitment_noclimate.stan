data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Y[N]; // observation vector
  vector[N] parents1; // crowding vector
  vector[N] parents2; // crowding vector
}
parameters{
  real a_mu;
  vector[Yrs] a;
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
vector[N] q;
  for(n in 1:N){
    trueP1[n] <- parents1[n]*u + parents2[n]*(1-u);
    trueP2[n] <- sqrt(trueP1[n]);
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + dd*trueP2[n]);
    lambda[n] <- trueP1[n]*mu[n];
    q[n] <- lambda[n]*theta;
  }
}
model{
// Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,100);
  a_mu ~ normal(0,1000);
  dd ~ uniform(-100,100);
  sig_a ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  for(g in 1:G)
   gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
  }

// Likelihood
  Y ~ neg_binomial_2(q, theta);
}
