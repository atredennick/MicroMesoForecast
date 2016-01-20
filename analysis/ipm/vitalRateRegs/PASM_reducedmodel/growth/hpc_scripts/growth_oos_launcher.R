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
####  Load Climate Covariates, Crowding Data, and Growth Data
####
climD <- read.csv("Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
# climD_raw <- climD[,c("year",clim_vars)]
# climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
# colnames(climD_raw)[2:ncol(climD_raw)] <- paste0(colnames(climD_raw)[2:ncol(climD_raw)],"_raw")
climD$year <- climD$year-1900
# climD_raw$year <- climD_raw$year-1900

crowd <- read.csv(paste0(do_species,"growthCrowding.csv"))[,2:3]
crowd$species <- do_species
colnames(crowd) <- c("W", "X", "species")

growD <- read.csv(paste0(do_species, "/growDnoNA.csv"))
growD$species <- do_species
growD <- merge(growD,climD)
# growD <- merge(growD,climD_raw)
growD$Group <- as.factor(substr(growD$quad,1,1))

# Merge crowding and growth data
growD <- merge(growD, crowd, by=c("species", "X"))



####
####  Compile Stan Model
####
holdyear <- unique(growD$year)[1] # grab first year as holdout
trainD <- subset(growD, year!=holdyear)
##  Create and scale interaction covariates
trainD$ppt1TmeanSpr1 <- trainD$ppt1*trainD$TmeanSpr1
trainD$ppt2TmeanSpr2 <- trainD$ppt2*trainD$TmeanSpr2
trainD$sizepptLag <- trainD$pptLag*log(trainD$area.t0)
trainD$sizeppt1 <- trainD$ppt1*log(trainD$area.t0)
trainD$sizeppt2 <- trainD$ppt2*log(trainD$area.t0)
trainD$sizeTmeanSpr1 <- trainD$TmeanSpr1*log(trainD$area.t0)
trainD$sizeTmeanSpr2 <- trainD$TmeanSpr2*log(trainD$area.t0)
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                   "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
clim_covs <- trainD[,clim_vars_all]
# Get scalers for climate covariates from training data
clim_means <- colMeans(clim_covs)
clim_sds <- apply(clim_covs, 2, FUN = sd)
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
groups <- as.numeric(trainD$Group)
G <- length(unique(trainD$Group))
nyrs <- length(unique(trainD$year))
W <- cbind(trainD$W, trainD$W*log(trainD$area.t0))
yid <- as.numeric(as.factor(trainD$year))

# Holdout data
holdD <- subset(growD, year==holdyear)
holdD$ppt1TmeanSpr1 <- holdD$ppt1*holdD$TmeanSpr1
holdD$ppt2TmeanSpr2 <- holdD$ppt2*holdD$TmeanSpr2
holdD$sizepptLag <- holdD$pptLag*log(holdD$area.t0)
holdD$sizeppt1 <- holdD$ppt1*log(holdD$area.t0)
holdD$sizeppt2 <- holdD$ppt2*log(holdD$area.t0)
holdD$sizeTmeanSpr1 <- holdD$TmeanSpr1*log(holdD$area.t0)
holdD$sizeTmeanSpr2 <- holdD$TmeanSpr2*log(holdD$area.t0)
clim_covs_oos <- holdD[,clim_vars_all]
for(j in 1:ncol(clim_covs_oos)){
  clim_covs_oos[,j] <- (clim_covs_oos[,j] - clim_means[j])/clim_sds[j]
}
W_oos <- cbind(holdD$W, holdD$W*log(holdD$area.t0))

datalist <- list(N=nrow(trainD), Yrs=nyrs, yid=yid,
                 Covs=ncol(clim_covs), Y=log(trainD$area.t1), X=log(trainD$area.t0),
                 C=clim_covs, W=W, G=G, gid=groups, tau_beta=1,
                 npreds=nrow(holdD), y_holdout=log(holdD$area.t1), Xhold=log(holdD$area.t0),
                 Chold=clim_covs_oos, Whold=W_oos)
pars <- c("log_lik")
mcmc_oos <- stan(file="growth_oos_cv.stan", data=datalist, 
                 pars=pars, chains=0)



####
####  Set up CV x Regularization grid
####
n.beta <- 24
sd_vec <- seq(0.1,1.5,length.out = n.beta)
yrs.vec <- unique(growD$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)



####
####  Source Growth Model Function and Fit Model
####
source("growth_fxns.R")
out_lpd <- cv.fcn(do_grid)
saveRDS(out_lpd, paste0(do_species,"_oos_cv_dogrid_",do_grid,".RDS"))

