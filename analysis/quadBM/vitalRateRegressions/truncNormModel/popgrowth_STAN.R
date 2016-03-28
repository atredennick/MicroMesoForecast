## Script to estimate parameters for the quad-based model using STAN

rm(list=ls())

library(rstan)
library(parallel)
library(ggmcmc)

priors_df <- read.csv("../../../all_maxlppds.csv")
priors <- subset(priors_df, vital=="cover")

##  STAN model
model_string <- "
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  real<lower=0,upper=1> Y[N]; // observation vector
  matrix[N,Covs] C; // climate matrix
  vector[N] X; // size vector
  real sdclim; // prior for climate effects
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
  b2 ~ normal(0,sdclim);
  a ~ normal(a_mu, sig_a);
  b1 ~ normal(b1_mu, sig_b1);
  tau ~ cauchy(0,5);

  //Likelihood
  Y ~ lognormal(mu, tau);
}
"

##  Read in data
#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")

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


growD <- subset(growD_all, Species=="BOGR")
##  Create and scale interaction covariates
growD$ppt1TmeanSpr1 <- growD$ppt1*growD$TmeanSpr1
growD$ppt2TmeanSpr2 <- growD$ppt2*growD$TmeanSpr2
growD$sizepptLag <- growD$pptLag*log(growD$percLagCover)
growD$sizeppt1 <- growD$ppt1*log(growD$percLagCover)
growD$sizeppt2 <- growD$ppt2*log(growD$percLagCover)
growD$sizeTmeanSpr1 <- growD$TmeanSpr1*log(growD$percLagCover)
growD$sizeTmeanSpr2 <- growD$TmeanSpr2*log(growD$percLagCover)
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                   "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
clim_covs <- growD[,clim_vars_all]
# Get scalers for climate covariates from training data
clim_means <- colMeans(clim_covs)
clim_sds <- apply(clim_covs, 2, FUN = sd)
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
groups <- as.numeric(as.factor(growD$group))
G <- length(unique(growD$group))
Yrs <- length(unique(growD$year))
yid <- as.numeric(as.factor(growD$year))

datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                 Covs=ncol(clim_covs), Y=growD$percCover, X=log(growD$percLagCover),
                 C=clim_covs, G=G, gid=groups, sdclim=0.1)
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "tau", "gint")

mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=0)

##  Loop through species and fit the model
for (do_species in sppList){
  print(paste("fitting model for", do_species))
  
  growD <- subset(growD_all, Species==do_species)
  prior_stddev <- as.numeric(subset(priors, species==do_species)["prior_stdev"])
  
  ##  Create and scale interaction covariates
  growD$ppt1TmeanSpr1 <- growD$ppt1*growD$TmeanSpr1
  growD$ppt2TmeanSpr2 <- growD$ppt2*growD$TmeanSpr2
  growD$sizepptLag <- growD$pptLag*log(growD$percLagCover)
  growD$sizeppt1 <- growD$ppt1*log(growD$percLagCover)
  growD$sizeppt2 <- growD$ppt2*log(growD$percLagCover)
  growD$sizeTmeanSpr1 <- growD$TmeanSpr1*log(growD$percLagCover)
  growD$sizeTmeanSpr2 <- growD$TmeanSpr2*log(growD$percLagCover)
  clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                     "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
  clim_covs <- growD[,clim_vars_all]
  # Get scalers for climate covariates from training data
  clim_means <- colMeans(clim_covs)
  clim_sds <- apply(clim_covs, 2, FUN = sd)
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  groups <- as.numeric(as.factor(growD$group))
  G <- length(unique(growD$group))
  Yrs <- length(unique(growD$year))
  yid <- as.numeric(as.factor(growD$year))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=c(0,0), sig_b1=0.5, sig_a=0.5, tau=0.5,
                     sig_G=0.5, b2=rep(0,length(clim_covs)))
  inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                     gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, tau=1,
                     sig_G=1, b2=rep(1,length(clim_covs)))
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                     gint=rep(0.5,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                     sig_G=0.1, b2=rep(-1,length(clim_covs)))
  
  datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=growD$percCover, X=log(growD$percLagCover),
                   C=clim_covs, G=G, gid=groups, sdclim=prior_stddev)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "tau", "gint")
  
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  long <- ggs(fit)
  saveRDS(long, paste("popgrowth_stanmcmc_", do_species, ".RDS", sep=""))
}

