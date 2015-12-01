## Script to estimate parameters for the quad-based model using STAN

library(rstan)
library(parallel)
library(ggmcmc)
library(snowfall)
library(rlecuyer)
source("../../../waic_fxns.R")

##  STAN model
model_string <- "
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
  real<lower=0> tau;
}
transformed parameters{
  real mu[N];
  vector[N] climEff;
  climEff <- C*b2;
  for(n in 1:N)
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + climEff[n];
}
model{
  // Priors
  a_mu ~ normal(0,1000);
  b1_mu ~ normal(0,1000);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  gint ~ normal(0, sig_G);
  b2 ~ normal(0,sd_clim);
  a ~ normal(a_mu, sig_a);
  b1 ~ normal(b1_mu, sig_b1);
  tau ~ cauchy(0,5);

  //Likelihood
  Y ~ lognormal(mu, tau);
}
generated quantities {
  vector[npreds] climpred;
  real muhat[npreds]; // prediction vector
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  climpred <- C_out*b2;
  for(n in 1:npreds){
    muhat[n] <- a_mu + b1_mu*X_out[n] + climpred[n];
    log_lik[n] <- lognormal_log(y_holdout[n], muhat[n], tau);
  }
}
"

##  Read in data
#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)

backD <- data.frame(climYear=NA,
                    quad = NA,
                    year= NA,
                    totCover= NA,
                    Species= NA,
                    propCover= NA,
                    lag.cover= NA,
                    pptLag= NA,
                    ppt1= NA,
                    TmeanSpr1= NA,
                    ppt2= NA,
                    TmeanSpr2= NA,
                    TmeanSum1= NA,
                    TmeanSum2= NA,
                    yearID= NA,
                    group = NA,
                    percCover = NA,
                    percLagCover = NA)

#loop through species and remake data frame
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover variable
  tmp=sppD[,c("quad","year","totCover")]
  tmp$year=tmp$year+1
  names(tmp)[3]="lag.cover"
  sppD=merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear=sppD$year+1900-1  
  sppD=merge(sppD,climD,by.x="climYear",by.y="year")
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  growD$yearID <- growD$year #for random year offset on intercept
  growD$group <- substring(growD$quad, 1, 1)
  growD$percCover <- growD$totCover/10000
  growD$percLagCover <- growD$lag.cover/10000
  backD <- rbind(backD, growD)
}#end species loop
growD_all <- backD[2:nrow(backD),]


growD <- subset(growD_all, Species=="PASM" & year != 33)
clim_covs <- growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups <- as.numeric(as.factor(growD$group))
G <- length(unique(growD$group))
Yrs <- length(unique(growD$year))
yid <- as.numeric(as.factor(growD$year))

hold_data <- subset(growD_all, Species=="PASM" & year == 33)
clim_covs_out <- hold_data[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]


datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=growD$percCover, X=log(growD$percLagCover),
                 C=clim_covs, G=G, gid=groups, sd_clim=0.1,
                 y_holdout=hold_data$percCover, X_out=log(hold_data$percLagCover),
                 C_out=clim_covs_out, npreds=nrow(hold_data))
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "tau", "gint", "log_lik")

mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=1, iter = 100, warmup = 50)


waic_metrics <- waic(mcmc_samples)
lpd <- waic_metrics[["total"]]["lpd"]


####
####  Set up cross validation and regularization loop
####
grow_pasm <- subset(growD_all, Species=="PASM")
n.beta <- 24
sd_vec <- seq(.25,1.5,length.out = n.beta)
yrs.vec <- unique(subset(growD_all, Species=="PASM")$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)

clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")

cps=detectCores()
sfInit(parallel=TRUE, cpus=cps)
sfExportAll()
sfClusterSetupRNG()

cv.fcn <- function(i){
  library(rstan)
  library(matrixStats)
  k <- cv.s2.grid[i,2]
  fold.idx <- fold.idx.mat[,k] 
  yr.out <- yrs.vec[fold.idx]
  sd.now <- sd_vec[cv.s2.grid[i,1]]
  df_train <- subset(grow_pasm, year!=yr.out)
  df_hold <- subset(grow_pasm, year==yr.out) 
  clim_covs <- df_train[,clim_vars]
  groups <- as.numeric(as.factor(df_train$group))
  G <- length(unique(df_train$group))
  Yrs <- length(unique(df_train$year))
  yid <- as.numeric(as.factor(df_train$year))
  clim_covs_out <- df_hold[,clim_vars]
  
  datalist <- list(N=nrow(df_train), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=df_train$percCover, 
                   X=log(df_train$percLagCover),
                   C=clim_covs, G=G, gid=groups, sd_clim=sd.now,
                   y_holdout=df_hold$percCover, X_out=log(df_hold$percLagCover),
                   C_out=clim_covs_out, npreds=nrow(df_hold))
  pars <- c("log_lik")
  fit <- stan(fit = mcmc_samples, data=datalist,
              pars=pars, chains=1, iter = 200, warmup = 100)
  waic_metrics <- waic(fit)
  lpd <- waic_metrics[["total"]]["elpd_loo"]
  return(lpd)
}

sfExport("sd_vec","cv.s2.grid","cv.fcn",
         "fold.idx.mat","yrs.vec","grow_pasm")
tmp.time=Sys.time()
score.list=sfClusterApplySR(1:n.grid,cv.fcn,perUpdate=1)
time.2=Sys.time()-tmp.time
sfStop()

score.cv.mat=matrix(unlist(score.list),n.beta,K)


####
####  Shrinkage Trajectories and CV Score
####
score.cv.vec=apply(score.cv.mat,1,sum)
sd.beta.opt <- sd_vec[which(score.cv.vec==max(score.cv.vec))]
plot(sd_vec^2,score.cv.vec,type="l",lwd=3,ylab="c-v score",
     xlab=bquote(sigma[beta]^2))
abline(v=sd.beta.opt^2,col=8,lwd=2)

