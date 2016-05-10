#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all.names=TRUE))

## Set leave_out_year for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
leave_out_year_id <- as.numeric(myargument)

#load libraries
library(rstan)
library(parallel)
library(ggmcmc)

sppList=sort(c("BOGR","HECO","PASM","POSE"))

data_path <- "./" # on HPC
# data_path <- "../../../../processed_data/" # on PC
growD_all <- readRDS(paste0(data_path, "grow_with_weather.RDS"))

# Get vector of all possible years
year_ids <- sort(unique(growD_all$year))


## Compile model outside of loop
yearD_all <- subset(growD_all, year!=year_ids[leave_out_year_id]) # remove lo year
growD <- subset(yearD_all, species==sppList[1]) # get just one species

groups <- as.numeric(growD$Group)
G <- length(unique(growD$Group))
nyrs <- length(unique(growD$year))
W <- cbind(growD$W, growD$W*log(growD$area.t0))
yid <- as.numeric(as.factor(growD$year))

datalist <- list(N=nrow(growD), Yrs=nyrs, yid=yid,
                 Y=log(growD$area.t1), X=log(growD$area.t0),
                 W=W, G=G, gid=groups)
pars=c("a_mu", "a", "b1_mu",  "b1",
       "w", "gint", "tau", "tauSize")
mcmc_samples <- stan(file="growth_noclimate.stan", data=datalist, pars=pars, chains=0)


## Loop through and fit each species' model
for(do_species in sppList){
#   print(paste("fitting model for", do_species, sep=""))
  growD <- subset(yearD_all, species==do_species)
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
  
  yid <- as.numeric(as.factor(growD$year))
  datalist <- list(N=nrow(growD), Yrs=nyrs, yid=yid,
                   Y=log(growD$area.t1), X=log(growD$area.t0),
                   W=W, G=G, gid=groups, tau_beta=prior_stddev)
  pars=c("a_mu", "a", "b1_mu",  "b1",
         "w", "gint", "tau", "tauSize")
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  r_hats <- summary(fit)$summary[,10] 
  write.csv(r_hats, paste("rhat_leaveout_noclimate", year_ids[leave_out_year_id], "_", do_species, ".csv", sep=""))
  long <- ggs(fit)
  saveRDS(long, paste("growth_stanmcmc_noclimate_", do_species, "_leaveout", year_ids[leave_out_year_id],".RDS", sep=""))
} # end species loop




