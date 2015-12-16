##
##  R Script to Test LPPD Regularization on Simulated Data
##

library(matrixStats)
library(rstan)
source("waic_fxns.R")


####
####  Simulate Data
####
b0 <- 10
b1 <- 1.5
b2 <- 0.4

nobs <- 200
x1 <- rnorm(nobs,0,1)
x2 <- rnorm(nobs,0,1)
y <- rnorm(nobs, b0+b1*x1+b2*x2, 2)
par(mfrow=c(1,2))
plot(x1,y)
plot(x2,y)
summary(lm(y~x1+x2))



####
####  Write Stan Model
####
stan_string <- "
data{
  // Training Data
  int<lower=0> ntrain; // number of observations
  int<lower=0> k; // number of covariates
  matrix[ntrain,k] Xtrain; // design matrix, minus intercept
  vector[ntrain] ytrain; // observations
  real beta_sdev; // prior sdev for covariate effects

  // Holdout Data
  int<lower=0> nhold; // number of holdout observations
  matrix[nhold,k] Xhold; // design matrix, minus intercept
  vector[nhold] yhold; // holdout observations
}
parameters{
  real b0;
  vector[k] b;
  real<lower=0.00001> sdev_mod;
}
model{
  // Priors
  b0 ~ normal(0, 1000);
  b ~ normal(0, beta_sdev);
  sdev_mod ~ cauchy(0,5);

  // Likelihood
  ytrain ~ normal(b0+Xtrain*b, sdev_mod);
}
generated quantities{
  vector[nhold] log_lik; // vector for computing log pointwise predictive density
  vector[nhold] muhat;
  vector[nhold] cov_eff;
  cov_eff <- Xhold*b;
  for(i in 1:nhold){
    muhat[i] <- b0+cov_eff[i];
    log_lik[i] <- normal_log(yhold[i], muhat[i], sdev_mod);
  }
}
"



####
####  Initialize Stan Model
####
trainids <- sample(c(1:length(y)), .75*length(y), replace = FALSE)
X <- data.frame(x1=x1, x2=x2)
ytrain <- y[trainids]
yhold <- y[-trainids]
Xtrain <- X[trainids,]
Xhold <- X[-trainids,]
datalist <- list(ytrain=ytrain, yhold=yhold, Xtrain=Xtrain, Xhold=Xhold,
                 ntrain=length(ytrain), nhold=length(yhold), k=ncol(X),
                 beta_sdev=1)
pars <- "b"
stan_init <- stan(model_code = stan_string, data = datalist, chains = 0, pars=pars)



sdev_vec <- seq(0.1, 2, length.out = 100)
lppd <- numeric(24)
for(i in 1:length(sdev_vec)){
  datalist <- list(ytrain=ytrain, yhold=yhold, Xtrain=Xtrain, Xhold=Xhold,
                   ntrain=length(ytrain), nhold=length(yhold), k=ncol(X),
                   beta_sdev=sdev_vec[i])
  pars <- c("b", "log_lik")
  fit <- stan(fit=stan_init, data = datalist, chains = 1, pars=pars,
              iter = 2000, warmup = 1000)
  waics <- waic(fit)
  lppd[i] <- waics[["total"]][["elpd_loo"]]
}

par(mfrow=c(1,1))
plot(sdev_vec, lppd, type="l")




