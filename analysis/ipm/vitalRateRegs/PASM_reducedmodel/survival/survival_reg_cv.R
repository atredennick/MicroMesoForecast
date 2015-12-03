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
## Initialize shrinkage model with full data set
survD <- subset(survD_all, species==sppList[do_species])
clim_covs <- survD[,clim_vars]
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
Yrs <- length(unique(survD$year))
W <- cbind(survD$W, survD$W*log(survD$area))
yid <- as.numeric(as.factor(survD$year))

datalist <- list(N=nrow(survD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=survD$survives, X=log(survD$area),
                 C=clim_covs, W=W, G=G, gid=groups, beta_tau=1)
pars <- c("b2")
mcmc_reg <- stan(file="survival_reg_cv.stan", data=datalist, pars=pars, chains=0)

##  Initialize out-of-sample cross validation model
holdyear <- unique(survD$year)[1]
survD <- subset(survD_all, species==sppList[do_species] & year!=holdyear)
clim_covs <- survD[,clim_vars]
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
Yrs <- length(unique(survD$year))
W <- cbind(survD$W, survD$W*log(survD$area))
yid <- as.numeric(as.factor(survD$year))

holdD <- subset(survD_all, species==sppList[do_species] & year==holdyear)
clim_covs_oos <- holdD[,clim_vars]
W_oos <- cbind(holdD$W, holdD$W*log(holdD$area))

datalist <- list(N=nrow(survD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=survD$survives, X=log(survD$area),
                 C=clim_covs, W=W, G=G, gid=groups, beta_tau=1,
                 npreds=nrow(holdD), yhold=holdD$survives, Xhold=log(holdD$area),
                 Chold=clim_covs_oos, Whold=W_oos)
pars <- c("b2")
mcmc_oos <- stan(file="survival_oos_cv.stan", data=datalist, pars=pars, chains=0)



####
####  Use Parallel Regularization Fit with Full Data Sets 
####
####  Issue this command in shell before starting R: export OMP_NUM_THREADS=1 
####
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
  clim_covs <- survD[,clim_vars]
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
              pars=pars, chains=1, iter = 2000, warmup = 1000)
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



####
####  Fit Each Logit Survival Model using Parallelized C-V 
####
####  Issue this command in shell before starting R: export OMP_NUM_THREADS=1 
####
survD <- subset(survD_all, species==sppList[do_species])
n.beta <- 24
sd_vec <- seq(.1,1.5,length.out = n.beta)
yrs.vec <- unique(survD$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)

cps=detectCores()
sfInit(parallel=TRUE, cpus=cps)
sfExportAll()
sfClusterSetupRNG()

cv.fcn <- function(i){
  library(rstan)
  library(matrixStats)
  k <- cv.s2.grid[i,2]
  fold.idx <- fold.idx.mat[,k] 
  yr.out <- yrs.vec[fold.idx]
  sd.now <- sd_vec[cv.s2.grid[i,1]]
  df_train <- subset(survD, year!=yr.out)
  df_hold <- subset(survD, year==yr.out) 
  clim_covs <- df_train[,clim_vars]
  groups <- as.numeric(df_train$Group)
  G <- length(unique(df_train$Group))
  Yrs <- length(unique(df_train$year))
  W <- cbind(df_train$W, df_train$W*log(df_train$area))
  yid <- as.numeric(as.factor(df_train$year))
  sd.now <- sd_vec[cv.s2.grid[i,1]]
  
  clim_covs_oos <- df_hold[,clim_vars]
  W_oos <- cbind(df_hold$W, df_hold$W*log(df_hold$area))
  
  datalist <- list(N=nrow(df_train), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=df_train$survives, X=log(df_train$area),
                   C=clim_covs, W=W, G=G, gid=groups, beta_tau=sd.now,
                   npreds=nrow(df_hold), yhold=df_hold$survives, Xhold=log(df_hold$area),
                   Chold=clim_covs_oos, Whold=W_oos)
  pars <- c("b2")
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=c(-0.05,-0.05), sig_b1=0.5, sig_a=0.5,
                     sig_G=0.5, b2=rep(0,length(clim_covs)))
  fit <- stan(fit = mcmc_oos, data=datalist, init=list(inits[[1]]),
              pars=pars, chains=1, iter = 200, warmup = 100)
  waic_metrics <- waic(fit)
  lpd <- waic_metrics[["total"]]["elpd_loo"]
  return(lpd)
}

sfExport("sd_vec","cv.s2.grid","cv.fcn",
         "fold.idx.mat","yrs.vec","survD")
tmp.time=Sys.time()
score.list=sfClusterApplySR(1:n.grid,cv.fcn,perUpdate=1)
time.2=Sys.time()-tmp.time
sfStop()
time.2
score.cv.mat=matrix(unlist(score.list),n.beta,K)



####
####  Plot Shrinkage Trajectories and CV Score
####
score.cv.vec <- apply(score.cv.mat,1,sum)
sd.beta.opt <- sd_vec[which(score.cv.vec==max(score.cv.vec))]
png("cv_score_surv.png", width = 6, height=8, units = "in", res=100)
par(mfrow=c(2,1), mar=c(1,4.1,4.1,2.1))
plot(sd_vec^2,beta.mat[,1],type="l",lwd=3,ylab=bquote(beta[mean]),
     xlab=bquote(sigma[beta]^2), ylim=c(-1,1))
for(i in 2:ncol(beta.mat)){
  lines(sd_vec^2,beta.mat[,i],type="l",lwd=3,col=i)
}
abline(h=0, lty=2)
abline(v=sd.beta.opt^2,col=8,lwd=2)
legend("bottomright",legend = names(clim_covs), col = c(1:ncol(beta.mat)), 
       lty = 1, bty="n", ncol=2)
par(mar=c(5,4.1,2,2.1))
plot(sd_vec^2,score.cv.vec,type="l",lwd=3,ylab="Log Predictive Score (lppd)",
     xlab=bquote(sigma[beta]^2))
abline(v=sd.beta.opt^2,col=8,lwd=2)
dev.off()

