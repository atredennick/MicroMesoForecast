## Script to estimate parameters for the quad-based model using STAN

rm(list=ls())

library(rstan)
library(parallel)
library(ggmcmc)

priors_df <- read.csv("../../../all_maxlppds.csv")
priors <- subset(priors_df, vital=="cover")

##  Read in data
growD_all <- readRDS("../../../processed_data/cover_with_weather.RDS")
sppList <- unique(growD_all$species)



####
####  COMPILE STAN MODEL
####
growD <- subset(growD_all, species=="BOGR")
##  Create and scale interaction covariates
clim_vars_all <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "ppt1TmeanSpr1", "ppt2TmeanSpr2")
clim_covs <- growD[,clim_vars_all]
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)

groups <- as.numeric(as.factor(growD$group))
G <- length(unique(growD$group))
Yrs <- length(unique(growD$year))
yid <- as.numeric(as.factor(growD$year))

datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                 Covs=ncol(clim_covs), Y=growD$propCover.t1, X=log(growD$propCover.t0),
                 C=clim_covs, G=G, gid=groups, sd_clim=0.1)
pars=c("a_mu", "a", "b1_mu",  "b1",
       "tau", "gint", "sig_a", "sig_b1", "sig_G")

mcmc_samples <- stan(file="qbm.stan", data=datalist, pars=pars, chains=0)

##  Loop through species and fit the model
for (do_species in sppList){
  print(paste("fitting model for", do_species))
  
  growD <- subset(growD_all, species==do_species)
  prior_stddev <- as.numeric(subset(priors, species==do_species)["prior_stdev"])
  
  ##  Create and scale interaction covariates
  clim_covs <- growD[,clim_vars_all]
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  groups <- as.numeric(as.factor(growD$group))
  G <- length(unique(growD$group))
  Yrs <- length(unique(growD$year))
  yid <- as.numeric(as.factor(growD$year))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=c(0,0), sig_b1=0.5, sig_a=0.5, sigmaSq=0.5,
                     sig_G=0.5, b2=rep(0,ncol(clim_covs)))
  inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                     gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, sigmaSq=1,
                     sig_G=1, b2=rep(1,ncol(clim_covs)))
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                     gint=rep(0.5,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, sigmaSq=0.1,
                     sig_G=0.1, b2=rep(-1,ncol(clim_covs)))
  
  
  datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                   Covs=ncol(clim_covs), Y=growD$propCover.t1, X=log(growD$propCover.t0),
                   C=clim_covs, G=G, gid=groups, sd_clim=prior_stddev)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "tau", "gint", "sig_a", "sig_b1", "sig_G")
  
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  long <- ggs(fit)
  r_hats <- summary(fit)$summary[,10] 
  write.csv(r_hats, paste("rhat_", do_species, ".csv", sep=""))
  saveRDS(long, paste("popgrowth_stanmcmc_", do_species, ".RDS", sep=""))
}

