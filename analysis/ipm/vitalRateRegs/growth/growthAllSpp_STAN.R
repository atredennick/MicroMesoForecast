#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rstan)
library(parallel)
library(ggmcmc)

##  Read in lppd scores and selected prior climate stddevs
priors_df <- read.csv("../../../all_maxlppds.csv")
priors <- subset(priors_df, vital=="growth")

sppList=sort(c("BOGR","HECO","PASM","POSE"))

####
#### Read in data by species and make one long data frame -------------
####
outD <- data.frame(X=NA,
                   quad=NA,
                   year=NA,
                   trackID=NA,
                   area.t1=NA,
                   area.t0=NA,
                   age=NA,
                   allEdge=NA,
                   distEdgeMin=NA,
                   species=NA)

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/edited/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900

growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRgrowthCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOgrowthCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMgrowthCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEgrowthCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)

# Merge crowding and growth data
growD_all <- merge(growD, crowd, by=c("species", "X"))

## Compile model outside of loop
growD <- subset(growD_all, species==sppList[1]) # grab a species

##  Create and scale interaction covariates
growD$ppt1TmeanSpr1 <- growD$ppt1*growD$TmeanSpr1
growD$ppt2TmeanSpr2 <- growD$ppt2*growD$TmeanSpr2
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
clim_covs <- growD[,clim_vars_all]
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE) # center and scale

##  Create objects for other effects
groups <- as.numeric(growD$Group)
G <- length(unique(growD$Group))
nyrs <- length(unique(growD$year))
W <- cbind(growD$W, growD$W*log(growD$area.t0))
yid <- as.numeric(as.factor(growD$year))

datalist <- list(N=nrow(growD), Yrs=nyrs, yid=yid,
                 Covs=ncol(clim_covs), Y=log(growD$area.t1), X=log(growD$area.t0),
                 C=clim_covs, W=W, G=G, gid=groups, tau_beta=0.1)
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "w", "gint", "tau", "tauSize")
mcmc_samples <- stan(file="growth.stan", data=datalist, pars=pars, chains=0)


## Loop through and fit each species' model
for(do_species in sppList){
  print(paste("fitting model for", do_species))
  
  prior_stddev <- as.numeric(subset(priors, species==do_species)["prior_stdev"])
  
  growD <- subset(growD_all, species==do_species)
  ##  Create and scale interaction covariates
  growD$ppt1TmeanSpr1 <- growD$ppt1*growD$TmeanSpr1
  growD$ppt2TmeanSpr2 <- growD$ppt2*growD$TmeanSpr2
  clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
  clim_covs <- growD[,clim_vars_all]
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  
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
                     sig_G=0.5, b2=rep(0,ncol(clim_covs)))
  inits[[2]] <- list(a_mu=1, a=rep(1,nyrs), b1_mu=1, b1=rep(1,nyrs),
                     gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, tau=1, tauSize=1,
                     sig_G=1, b2=rep(1,ncol(clim_covs)))
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,nyrs), b1_mu=0.5, b1=rep(0.5,nyrs),
                     gint=rep(0.5,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                     sig_G=0.1, b2=rep(-1,ncol(clim_covs)))
  
  ##  Create data list for Stan
  datalist <- list(N=nrow(growD), Yrs=nyrs, yid=yid,
                   Covs=ncol(clim_covs), Y=log(growD$area.t1), X=log(growD$area.t0),
                   C=clim_covs, W=W, G=G, gid=groups, tau_beta=prior_stddev)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "w", "gint", "tau", "tauSize")
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  
  long <- ggs(fit) # convert StanFit --> dataframe
  outfile <- paste("growth_stanmcmc_", do_species, ".RDS", sep="")
  saveRDS(long, outfile)
} # end species loop




