## Script to estimate parameters for the quad-based model using STAN


####
####  Load Libraries and Subroutines
####
library(rstan)
library(plyr)
library(reshape2)
library(ggmcmc)
library(matrixStats)
library(parallel)



##  Read in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))
climD <- read.csv("../../weather/Climate.csv")
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



####
####  Compile Stan Model
####
trainD <- subset(growD_all, Species=="BOGR")
##  Create and scale interaction covariates
trainD$ppt1TmeanSpr1 <- trainD$ppt1*trainD$TmeanSpr1
trainD$ppt2TmeanSpr2 <- trainD$ppt2*trainD$TmeanSpr2
trainD$sizepptLag <- trainD$pptLag*log(trainD$percLagCover)
trainD$sizeppt1 <- trainD$ppt1*log(trainD$percLagCover)
trainD$sizeppt2 <- trainD$ppt2*log(trainD$percLagCover)
trainD$sizeTmeanSpr1 <- trainD$TmeanSpr1*log(trainD$percLagCover)
trainD$sizeTmeanSpr2 <- trainD$TmeanSpr2*log(trainD$percLagCover)
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                   "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
clim_covs <- trainD[,clim_vars_all]
# Get scalers for climate covariates from training data
clim_means <- colMeans(clim_covs)
clim_sds <- apply(clim_covs, 2, FUN = sd)
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
groups <- as.numeric(as.factor(trainD$group))
G <- length(unique(trainD$group))
Yrs <- length(unique(trainD$year))
yid <- as.numeric(as.factor(trainD$year))

##  Initialize STAN model
datalist <- list(N=nrow(trainD), Yrs=Yrs, yid=yid,
                 Covs=ncol(clim_covs), Y=trainD$percCover, 
                 X=log(trainD$percLagCover),
                 C=clim_covs, G=G, gid=groups,
                 sd_clim=0.1)
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "tau", "gint")
stan_mod <- stan(file="qbm_reg_cv.stan", data=datalist,
                     pars=pars, chains=0)



##  Loop through species and fit the model
big_list <- list()
for (do_species in sppList){
  trainD <- subset(growD_all, Species==do_species)
  ##  Create and scale interaction covariates
  trainD$ppt1TmeanSpr1 <- trainD$ppt1*trainD$TmeanSpr1
  trainD$ppt2TmeanSpr2 <- trainD$ppt2*trainD$TmeanSpr2
  trainD$sizepptLag <- trainD$pptLag*log(trainD$percLagCover)
  trainD$sizeppt1 <- trainD$ppt1*log(trainD$percLagCover)
  trainD$sizeppt2 <- trainD$ppt2*log(trainD$percLagCover)
  trainD$sizeTmeanSpr1 <- trainD$TmeanSpr1*log(trainD$percLagCover)
  trainD$sizeTmeanSpr2 <- trainD$TmeanSpr2*log(trainD$percLagCover)
  clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                     "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
  clim_covs <- trainD[,clim_vars_all]
  # Get scalers for climate covariates from training data
  clim_means <- colMeans(clim_covs)
  clim_sds <- apply(clim_covs, 2, FUN = sd)
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  groups <- as.numeric(as.factor(trainD$group))
  G <- length(unique(trainD$group))
  Yrs <- length(unique(trainD$year))
  yid <- as.numeric(as.factor(trainD$year))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), sig_b1=0.5, sig_a=0.5, tau=0.5,
                     sig_G=0.5, b2=rep(0,ncol(clim_covs)))
  inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                     gint=rep(1,G), sig_b1=1, sig_a=1, tau=1,
                     sig_G=1, b2=rep(1,ncol(clim_covs)))
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                     gint=rep(0.5,G), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                     sig_G=0.1, b2=rep(-1,ncol(clim_covs)))
  
  datalist <- list(N=nrow(trainD), Yrs=Yrs, yid=yid,
                   Covs=ncol(clim_covs), Y=trainD$percCover, 
                   X=log(trainD$percLagCover),
                   C=clim_covs, G=G, gid=groups,
                   sd_clim=0.1)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "tau", "gint")
  
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=stan_mod, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  long <- ggs(fit)
  saveRDS(long, paste("popgrowth_stanmcmc_", do_species, ".RDS", sep=""))
}