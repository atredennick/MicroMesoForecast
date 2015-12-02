##
##  R Script for survival statistical model (for IPM)
##
##  Includes Model Selection via Bayesian Regularization
##
##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com
##  Date created: 12-2-2015
##


# Clear the workspace
rm(list=ls())


####
####  Load libraries and subroutines
####
library(rstan); library(plyr)
library(reshape2); library(ggmcmc)
library(parallel); library(rlecuyer)
library(snowfall)
source("../../../../waic_fxns.R")



####
####  Load Climate Covariates and Survival Data
####
sppList <- sort(c("BOGR","HECO","PASM","POSE"))
do_species <- 3 # for PASM, for testing

##  Read in data by species and make one long data frame
outD <- data.frame(X=NA,
                   quad=NA,
                   year=NA,
                   trackID=NA,
                   area=NA,
                   survives=NA,
                   age=NA,
                   distEdgeMin=NA,
                   allEdge=NA,
                   species=NA)

data_path <- "../../../../speciesData/" #on local machine
# data_path <- "speciesData/" #on HPC server

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste(data_path, doSpp, "/edited/survD.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste(data_path, doSpp, "/survD.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

survD <- outD[2:nrow(outD),]

##  Read in climate data
climD <- read.csv("../../../../weather/Climate.csv") #on local machine
# climD <- read.csv("Climate.csv") #on HPC server
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
survD <- merge(survD,climD)
survD$Group=as.factor(substr(survD$quad,1,1))

## Read in previously estimated crowding indices
c1 <- read.csv("BOGRsurvCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOsurvCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMsurvCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEsurvCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)
colnames(crowd) <- c("W", "X", "species")

# Merge crowding and growth data
survD_all <- merge(survD, crowd, by=c("species", "X"))



####
####  Collate Data for Stan Model Initialization
####
survD <- subset(survD_all, species==sppList[do_species])
clim_covs <- survD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
Yrs <- length(unique(survD$year))
W <- cbind(survD$W, survD$W*log(survD$area))
yid <- as.numeric(as.factor(survD$year))

## Initialize shrinkage model with full data set
datalist <- list(N=nrow(survD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=survD$survives, X=log(survD$area),
                 C=clim_covs, W=W, G=G, gid=groups, beta_tau=1)
pars <- c("b2")
mcmc_reg <- stan(file="survival_reg_cv.stan", data=datalist, pars=pars, chains=0)



####
####  Use Parallel Regularization Fit with Full Data Sets 
####
####  Issue this command in shell before starting R: export OMP_NUM_THREADS=1 
####
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
n.beta <- 24
sd_vec <- seq(.1,1.5,length.out = n.beta)

cps=detectCores()
sfInit(parallel=TRUE, cpus=cps)
sfExportAll()
sfClusterSetupRNG()

reg.fcn <- function(i){
  library(rstan)
  library(ggmcmc)
  library(plyr)
  clim_covs <- survD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  groups <- as.numeric(survD$Group)
  G <- length(unique(survD$Group))
  Yrs <- length(unique(survD$year))
  W <- cbind(survD$W, survD$W*log(survD$area))
  yid <- as.numeric(as.factor(survD$year))
  sd.now <- sd_vec[i]
  datalist <- list(N=nrow(survD), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=survD$survives, X=log(survD$area),
                   C=clim_covs, W=W, G=G, gid=groups, beta_tau=sd.now)
  pars <- c("b2")
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=c(-0.05,-0.05), sig_b1=0.5, sig_a=0.5,
                     sig_G=0.5, b2=rep(0,length(clim_covs)))
  fit <- stan(fit = mcmc_reg, data=datalist, init=list(inits[[1]]),
              pars=pars, chains=1, iter = 200, warmup = 100)
  long <- ggs(fit)
  longagg <- ddply(ggs(fit), .(Parameter), summarise,
                   avg_value = mean(value))
  return(longagg$avg_value)
}

sfExport("sd_vec", "survD", "clim_vars")
tmp.time <- Sys.time()
beta.list <- sfClusterApplySR(1:n.beta,reg.fcn,perUpdate=1)
time.1 <- Sys.time()-tmp.time
sfStop()
time.1

beta.mat <- matrix(unlist(beta.list),ncol=ncol(clim_covs),byrow=TRUE)

