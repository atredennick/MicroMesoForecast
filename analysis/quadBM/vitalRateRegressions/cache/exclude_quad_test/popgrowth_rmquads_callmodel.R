## Script to estimate parameters for the quad-based model using STAN

rm(list=ls())


library(rstan)
library(parallel)
library(ggmcmc)

####
####  Source preparation scripts
####
# source("popgrowth_STANmodel.R")
# source("popgrowth_STANmodel_nogroup.R")
source("popgrowth_read_data.R")

priors_df <- read.csv("../../../all_maxlppds.csv")
priors <- subset(priors_df, vital=="cover")


####
####  Loop over species
####
spps <- c("BOGR", "HECO", "PASM", "POSE")
for(do_spp in spps){
  ####
  ####  Subset dataframe for focal species (do_spp)
  ####
  growD_spp <- subset(growD_all, Species==do_spp)
  prior_stddev <- as.numeric(subset(priors, species==do_spp)["prior_stdev"])
  
  ####
  ####  Make list of quads to remove
  ####
  ##  Remove sets of 2, 5, 10, and 15 quadrats
  all_quads <- unique(growD_spp$quad)
  num_sets <- c(2,5,10,15)
  reps_per_removal <- 10
  quads_to_remove <- list()
  for(i in num_sets){
    tmp <- matrix(ncol=i, nrow=reps_per_removal)
    for(j in 1:reps_per_removal){
      tmp[j,] <- sample(all_quads, i, replace = FALSE)
    }
    quads_to_remove[[as.character(i)]] <- tmp
  }
  
  
  ####
  ####  Remove quads, fit models
  ####
  for(g in 1:length(quads_to_remove)){
    for(rp in 1:reps_per_removal){
      torm <- which(growD_spp$quad %in% quads_to_remove[[g]][rp,]) 
      growD <- growD_spp[-torm,]
      
      clim_covs <- growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
      clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
      clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
      clim_covs <- scale(clim_covs, center=TRUE, scale=TRUE)
      groups <- as.numeric(as.factor(growD$group))
      G <- length(unique(growD$group))
      Yrs <- length(unique(growD$year))
      yid <- as.numeric(as.factor(growD$year))
      
      datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                       Covs=ncol(clim_covs), Y=growD$percCover, X=log(growD$percLagCover),
                       C=clim_covs, G=G, gid=groups, sd_clim=prior_stddev)
      
      ##  Fit Model with Group Effect IF unique groups > 1
      if(G>1){
        pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
               "tau", "gint")
        inits <- list()
        inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                           gint=rep(0,G), w=c(0,0), sig_b1=0.5, sig_a=0.5, tau=0.5,
                           sig_G=0.5, b2=rep(0,ncol(clim_covs)))
        inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                           gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, tau=1,
                           sig_G=1, b2=rep(1,ncol(clim_covs)))
        inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                           gint=rep(0.5,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                           sig_G=0.1, b2=rep(-1,ncol(clim_covs)))
        mcmc_samples <- stan(file="qbm.stan", data=datalist, pars=pars, chains=0)
      } # end group effect IF
      if(G==1){
        pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
               "tau")
        inits <- list()
        inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                           w=c(0,0), sig_b1=0.5, sig_a=0.5, tau=0.5,
                           b2=rep(0,ncol(clim_covs)))
        inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                           w=c(0.5,0.5), sig_b1=1, sig_a=1, tau=1,
                           b2=rep(1,ncol(clim_covs)))
        inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                           w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                           b2=rep(-1,ncol(clim_covs)))
        mcmc_samples <- stan(file="qbm_nogroup.stan", data=datalist, pars=pars, chains=0)
      } # end no group effect IF
      
      ## Set reasonable initial values for three chains
      
      rng_seed <- 123
      sflist <-
        mclapply(1:3, mc.cores=3,
                 function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                                  seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                                  iter=2000, warmup=1000, init=list(inits[[i]])))
      fit <- sflist2stanfit(sflist)
      long <- ggs(fit)
      
      numout <- names(quads_to_remove)[g]
      outfile <- paste(do_spp,"_numout",numout,"_rep",rp, ".RDS", sep="")
      outdir <- paste("./results/",outfile,sep="")
      saveRDS(long, outdir)
    } # end rep loop
  } # end group drop loop
} # end species loop

