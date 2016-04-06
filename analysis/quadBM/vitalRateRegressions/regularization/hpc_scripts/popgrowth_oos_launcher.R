## Script to estimate parameters for the quad-based model using STAN


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



##  Read in data
allD <- read.csv("quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))
climD <- read.csv("Climate.csv")
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
  growD$percCover <- growD$totCover/10000 # proportional cover
  growD$percLagCover <- growD$lag.cover/10000 # proportional cover
  backD <- rbind(backD, growD)
}#end species loop
growD_all <- backD[2:nrow(backD),]



####
####  Compile Stan Model
####
holdyear <- unique(growD_all$year)[1] # grab first year as holdout
trainD <- subset(growD_all, Species==do_species & year != holdyear)
##  Create and scale interaction covariates
trainD$ppt1TmeanSpr1 <- trainD$ppt1*trainD$TmeanSpr1
trainD$ppt2TmeanSpr2 <- trainD$ppt2*trainD$TmeanSpr2
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
clim_covs <- trainD[,clim_vars_all]
# Get scalers for climate covariates from training data
clim_means <- colMeans(clim_covs)
clim_sds <- apply(clim_covs, 2, FUN = sd)
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
groups <- as.numeric(as.factor(trainD$group))
G <- length(unique(trainD$group))
Yrs <- length(unique(trainD$year))
yid <- as.numeric(as.factor(trainD$year))

holdD <- subset(growD_all, Species==do_species & year == holdyear)
holdD$ppt1TmeanSpr1 <- holdD$ppt1*holdD$TmeanSpr1
holdD$ppt2TmeanSpr2 <- holdD$ppt2*holdD$TmeanSpr2
clim_covs_oos <- holdD[,clim_vars_all]
for(j in 1:ncol(clim_covs_oos)){
  clim_covs_oos[,j] <- (clim_covs_oos[,j] - clim_means[j])/clim_sds[j] # scale by the training data
}

### Initialize Regularization MCMC
datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                 Covs=ncol(clim_covs), Y=growD$percCover, X=log(growD$percLagCover),
                 C=clim_covs, G=G, gid=groups, sd_clim=0.1)
pars=c("a_mu")
mcmc_reg <- stan(file="qbm_reg_cv.stan", data=datalist, pars=pars, chains=0)



####
####  Set up CV x Regularization grid
####
grow_data <- subset(growD_all, Species==do_species)
n.beta <- 24
sd_vec <- seq(0.1,1.5,length.out = n.beta)
yrs.vec <- unique(subset(growD_all, Species==do_species)$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)



####
####  Source Population Growth Model Function and Fit Model
####
source("popgrowth_fxns.R")
out_lpd <- cv.fcn(do_grid)
saveRDS(out_lpd, paste0(do_species,"_oos_cv_dogrid_",do_grid,".RDS"))




