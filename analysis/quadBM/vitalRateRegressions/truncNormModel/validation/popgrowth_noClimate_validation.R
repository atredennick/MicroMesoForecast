## Script to estimate parameters for the quad-based model using STAN

library(rstan)
library(parallel)
library(ggmcmc)

## Set leave_out_year for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
leave_out_year_id <- as.numeric(myargument)

sppList <- c("BOGR", "HECO", "PASM", "POSE")
growD_all <- readRDS("cover_with_weather.RDS")

year_ids <- unique(growD_all$year)
yearD_all <- subset(growD_all, year!=year_ids[leave_out_year_id])

growD <- subset(yearD_all, species=="BOGR")
groups <- as.numeric(as.factor(growD$group))
G <- length(unique(growD$group))
Yrs <- length(unique(growD$year))
yid <- as.numeric(as.factor(growD$year))

datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                 Y=growD$propCover.t1, X=log(growD$propCover.t0),
                 G=G, gid=groups)
pars=c("a_mu", "a", "b1_mu",  "b1",
       "tau", "gint")

mcmc_samples <- stan(file="qbm_noclimate.stan", data=datalist, pars=pars, chains=0)

##  Loop through species and fit the model
for (do_species in sppList){
  growD <- subset(yearD_all, species==do_species)
  groups <- as.numeric(as.factor(growD$group))
  G <- length(unique(growD$group))
  Yrs <- length(unique(growD$year))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=c(0,0), sig_b1=0.5, sig_a=0.5, tau=0.5,
                     sig_G=0.5)
  inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                     gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, tau=1,
                     sig_G=1)
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                     gint=rep(0.5,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                     sig_G=0.1)
  
  datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                   Y=growD$propCover.t1, X=log(growD$propCover.t0),
                   G=G, gid=groups)
  pars=c("a_mu", "a", "b1_mu",  "b1",
         "tau", "gint")
  
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  r_hats <- summary(fit)$summary[,10] 
  write.csv(r_hats, paste("rhat_leaveout", year_ids[leave_out_year_id], "_noclimate_", do_species, ".csv", sep=""))
  
  long <- ggs(fit)
  saveRDS(long, paste("popgrowth_stanmcmc_noclimate_", do_species, "_leaveout", year_ids[leave_out_year_id],".RDS", sep=""))
}
