##  The script takes 'do_species' as a command line prompt. So, e.g.,
##    run as: "R CMD BATCH -1 survivalAllSpp_STAN.R" for species 1.

#clear everything, just to be safe 
rm(list=ls(all.names=TRUE))

#load libraries
library(rstan)
library(parallel)

## Set do_year for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_species <- as.numeric(myargument)
do_species <- 1
sppList <- sort(c("BOGR","HECO","PASM","POSE"))

# data_path <- "../../../processed_data/" #on local machine
data_path <- "./" #on HPC server
survD_all <- readRDS(paste0(data_path,"surv_with_weather.RDS"))
survD <- subset(survD_all, species==sppList[do_species])

##  Create objects for other effects
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
nyrs <- length(unique(survD$year))
W <- cbind(survD$W, survD$W*log(survD$area))
yid <- as.numeric(as.factor(survD$year))

datalist <- list(N=nrow(survD), Yrs=nyrs, yid=yid,
                 Y=survD$survives, X=log(survD$area),
                 W=W, G=G, gid=groups)
pars=c("a_mu", "a", "b1_mu",  "b1",
       "w", "gint")

mcmc_samples <- stan(file="survival_noclimate.stan", data=datalist, pars=pars, chains=0)

## Set reasonable initial values for three chains
inits <- list()
inits[[1]] <- list(a_mu=0, a=rep(0,nyrs), b1_mu=0.01, b1=rep(0.01,nyrs),
                   gint=rep(0,G), w=c(-0.05,-0.05), sig_b1=0.5, sig_a=0.5,
                   sig_G=0.5, b2=rep(0,ncol(clim_covs)))
inits[[2]] <- list(a_mu=-0.1, a=rep(-0.1,nyrs), b1_mu=0.1, b1=rep(0.1,nyrs),
                   gint=rep(-0.1,G), w=c(-0.1,-0.1), sig_b1=0.2, sig_a=0.2,
                   sig_G=0.2, b2=rep(0.1,ncol(clim_covs)))
inits[[3]] <- list(a_mu=0.05, a=rep(0.05,nyrs), b1_mu=0.05, b1=rep(0.05,nyrs),
                   gint=rep(0.1,G), w=c(-0.05,-0.05), sig_b1=0.3, sig_a=0.3,
                   sig_G=0.3, b2=rep(0.05,ncol(clim_covs)))

##  Run MCMC in parallel
rng_seed <- 123
sflist <-
  mclapply(1:3, mc.cores=3,
          function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                           seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                           iter=2000, warmup=1000, init=list(inits[[i]])))
fit <- sflist2stanfit(sflist)
longfit <- ggs(fit) # convert StanFit --> data frame
outfile <- paste("survival_stanmcmc_noclimate_", sppList[do_species], ".RDS", sep="")
saveRDS(fit, outfile)
