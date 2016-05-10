#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all.names=TRUE))

#load libraries
library(rstan)
library(parallel)
library(reshape2)
library(ggmcmc)

##  Read in data
allD <- readRDS("../../../processed_data/rec_with_weather.RDS")
sppList <- c("BOGR", "HECO", "PASM", "POSE")

## Compile model outside of loop
recD <- subset(allD, species==paste("R.",sppList[1],sep=""))

##  Create and scale interaction covariates
groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))
yid <- as.numeric(as.factor(recD$year))

datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                 Y=recD$recruits, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups)
pars=c("a_mu", "a", "u", "theta",
       "dd", "gint")
mcmc_samples <- stan(file="recruitment_noclimate.stan", data=datalist, pars=pars, chains=0)


for(do_species in 1:length(sppList)){
  recD <- subset(allD, species==paste("R.",sppList[do_species],sep=""))
  groups <- as.numeric(recD$group)
  G <- length(unique(recD$group))
  Yrs <- length(unique(recD$year))
  yid <- as.numeric(as.factor(recD$year))

  datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                   Y=recD$recruits, 
                   parents1=recD$parents1, parents2=recD$parents2,
                   G=G, gid=groups)
  pars=c("a_mu", "a", "u", "theta",
         "dd", "gint")
  
  inits=list()
  inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                  gint=rep(0,G), sig_G=1, u=0.4,
                  dd=-1,theta=1)
  inits[[2]]=list(a=rep(0.5,Yrs), a_mu=0.2, sig_a=10,
                  gint=rep(0,G), sig_G=0.1,  u=0.7,
                  dd=-0.05,theta=1.5)
  inits[[3]]=list(a=rep(1,Yrs), a_mu=0.5, sig_a=5,
                  gint=rep(-0.1,G), sig_G=0.5,  u=0.5,
                  dd=-1,theta=1)

  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  longfit <- ggs(fit) # convert StanFit --> data frame
  saveRDS(fit, paste0("recruitment_stanmcmc_noclimate_",sppList[do_species],".RDS"))
}

