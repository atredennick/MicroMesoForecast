#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all.names=TRUE))

#load libraries
library(rstan)
library(parallel)
library(ggmcmc)

##  Read in data
growD_all <- readRDS("../../../processed_data/grow_with_weather.RDS")
sppList <- sort(unique(growD_all$species))

## Compile model outside of loop
growD <- subset(growD_all, species==sppList[1]) # grab a species

##  Create objects for model effects
groups <- as.numeric(growD$Group)
G <- length(unique(growD$Group))
nyrs <- length(unique(growD$year))
W <- cbind(growD$W, growD$W*log(growD$area.t0))
yid <- as.numeric(as.factor(growD$year))

datalist <- list(N=nrow(growD), Yrs=nyrs, yid=yid,
                 Y=log(growD$area.t1), X=log(growD$area.t0),
                 W=W, G=G, gid=groups)
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "w", "gint", "tau", "tauSize", "sig_a", "sig_b1", "sig_G")
mcmc_samples <- stan(file="growth_noclimate.stan", data=datalist, pars=pars, chains=0)


## Loop through and fit each species' model
for(do_species in sppList){
  cat(paste("\n\nfitting model for", do_species, "\n\n"))
  
  growD <- subset(growD_all, species==do_species)
  
  ##  Create objects for other effects
  groups <- as.numeric(growD$Group)
  G <- length(unique(growD$Group))
  nyrs <- length(unique(growD$year))
  W <- cbind(growD$W, growD$W*log(growD$area.t0))
  yid <- as.numeric(as.factor(growD$year))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,nyrs), b1_mu=0.01, b1=rep(0.01,nyrs),
                     gint=rep(0,G), w=c(0,0), sig_b1=0.5, sig_a=0.5, tau=0.5, tauSize=0.5,
                     sig_G=0.5)
  inits[[2]] <- list(a_mu=1, a=rep(1,nyrs), b1_mu=1, b1=rep(1,nyrs),
                     gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, tau=1, tauSize=1,
                     sig_G=1)
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,nyrs), b1_mu=0.5, b1=rep(0.5,nyrs),
                     gint=rep(0.5,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                     sig_G=0.1)
  
  ##  Create data list for Stan
  datalist <- list(N=nrow(growD), Yrs=nyrs, yid=yid,
                   Y=log(growD$area.t1), X=log(growD$area.t0),
                   W=W, G=G, gid=groups)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "w", "gint", "tau", "tauSize", "sig_a", "sig_b1", "sig_G")
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  
  long <- ggs(fit) # convert StanFit --> dataframe
  outfile <- paste("growth_stanmcmc_noclimate_", do_species, ".RDS", sep="")
  saveRDS(long, outfile)
  r_hats <- summary(fit)$summary[,10] 
  write.csv(r_hats, paste("rhat_noclimate_", do_species, ".csv", sep=""))
} # end species loop




