#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all.names = TRUE))

#load libraries
library(rstan)
library(parallel)
library(reshape2)
library(ggmcmc)

##  Read in lppd scores and selected prior climate stddevs
priors_df <- read.csv("../../../all_maxlppds.csv")
priors <- subset(priors_df, vital=="recruitment")

sppList=sort(c("BOGR","HECO","PASM","POSE"))

data_path <- "../../../processed_data/"
allD <- readRDS(paste0(data_path,"rec_with_weather.RDS"))


## Compile model outside of loop
recD <- subset(allD, species==paste("R.",sppList[1],sep=""))
##  Create and scale interaction covariates
clim_vars_all <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2", "ppt1TmeanSpr1", "ppt2TmeanSpr2")
clim_covs <- recD[,clim_vars_all]
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))
yid <- as.numeric(as.factor(recD$year))

datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                 Covs=ncol(clim_covs), Y=recD$recruits, C=clim_covs, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups, tau=0.1)
pars=c("a_mu", "a", "u", "theta", "dd", "gint", "sig_a", "sig_G", "b2")
mcmc_samples <- stan(file="recruitment.stan", data=datalist, pars=pars, chains=0)


big_list <- list()
for(do_species in 1:length(sppList)){
  recD <- subset(allD, species==paste("R.",sppList[do_species],sep=""))
  clim_covs <- recD[,clim_vars_all]
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  groups <- as.numeric(recD$group)
  G <- length(unique(recD$group))
  Yrs <- length(unique(recD$year))
  yid <- as.numeric(as.factor(recD$year))

  # Get climate prior std dev
  prior_stddev <- as.numeric(subset(priors, species==sppList[do_species])["prior_stdev"])
  
  datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                   Covs=ncol(clim_covs), Y=recD$recruits, C=clim_covs, 
                   parents1=recD$parents1, parents2=recD$parents2,
                   G=G, gid=groups, tau=prior_stddev)
  pars=c("a_mu", "a", "u", "theta", "dd", "gint", "sig_a", "sig_G", "b2")
  
  inits=list()
  inits[[1]]=list(a=rep(4,Yrs), a_mu=4, sig_a=0.1,
                  gint=rep(0.01,G), sig_G=0.01, u=0.4,
                  dd=-1,theta=1, b2=rep(0,ncol(clim_covs))) 
  inits[[2]]=list(a=rep(6,Yrs), a_mu=6, sig_a=0.25,
                  gint=rep(0,G), sig_G=0.1,  u=0.7,
                  dd=-0.05,theta=1.5, b2=rep(0.5,ncol(clim_covs))) 
  inits[[3]]=list(a=rep(10,Yrs), a_mu=10, sig_a=1,
                  gint=rep(-0.1,G), sig_G=0.5,  u=0.5,
                  dd=-2,theta=0.1, b2=rep(-0.5,ncol(clim_covs))) 

  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  longfit <- ggs(fit) # convert StanFit --> data frame
  saveRDS(fit, paste0("recruitment_stanmcmc_",sppList[do_species],".RDS"))
  r_hats <- summary(fit)$summary[,10] 
  write.csv(r_hats, paste("rhat_", do_species, ".csv", sep=""))
}

