##
##  R script for survival logistic model as run on Utah State University's
##  High Performance Computing system.
##
##  Script takes command line argument 'do_grid' to run all levels of
##  regularization and each cross-validation set in parallel. Script launches
##  n*k models: 'n' prior standard deviations and 'k' validation sets.
##
##  Saved output:
##    1. Summed log pointwise predictive density based on 
##       out-of-sample predictions
##
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 12-6-2015
##

# Change species four letter code here...
do_species <- "PASM"

####
####  Set SD Prior and CV Set from Command Line Arguments
####
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_grid <- as.numeric(myargument)


####
####  Load Libraries and Subroutines
####
library(rstan)
library(plyr)
library(reshape2)
library(ggmcmc)
library(matrixStats)
source("waic_fxns.R")



####
####  Load Climate Covariates, Crowding Data, and Survival Data
####
climD <- read.csv("Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900

crowd <- read.csv(paste0(do_species,"survCrowding.csv"))[,2:3]
crowd$species <- do_species
colnames(crowd) <- c("W", "X", "species")

survD <- read.csv(paste0(do_species, "/survD.csv"))
survD$species <- do_species
survD <- merge(survD,climD)
survD$Group <- as.factor(substr(survD$quad,1,1))

# Merge crowding and growth data
survD <- merge(survD, crowd, by=c("species", "X"))



####
####  Compile Stan Model
####
holdyear <- unique(survD$year)[1] # grab first year as holdout
trainD <- subset(survD, year!=holdyear)
clim_covs <- trainD[,clim_vars]
groups <- as.numeric(trainD$Group)
G <- length(unique(trainD$Group))
Yrs <- length(unique(trainD$year))
W <- cbind(trainD$W, trainD$W*log(trainD$area))
yid <- as.numeric(as.factor(trainD$year))

holdD <- subset(survD, year==holdyear)
clim_covs_oos <- holdD[,clim_vars]
W_oos <- cbind(holdD$W, holdD$W*log(holdD$area))

datalist <- list(N=nrow(trainD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=trainD$survives, X=log(trainD$area),
                 C=clim_covs, W=W, G=G, gid=groups, beta_tau=1,
                 npreds=nrow(holdD), yhold=holdD$survives, Xhold=log(holdD$area),
                 Chold=clim_covs_oos, Whold=W_oos)
pars <- c("log_lik")
mcmc_oos <- stan(file="survival_oos_cv.stan", data=datalist, 
                 pars=pars, chains=0)



####
####  Set up CV x Regularization grid
####
n.beta <- 24
sd_vec <- seq(0.1,1.5,length.out = n.beta)
yrs.vec <- unique(survD$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)



####
####  Source Survival Model Function and Fit Model
####
source("survival_fxns.R")
out_lpd <- cv.fcn(do_grid)
saveRDS(out_lpd, paste0("oos_cv_dogrid_",do_grid,".RDS"))

