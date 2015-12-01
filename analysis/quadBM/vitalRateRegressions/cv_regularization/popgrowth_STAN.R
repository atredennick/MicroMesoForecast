## Script to estimate parameters for the quad-based model using STAN

library(rstan)
library(parallel)
library(ggmcmc)
library(snowfall)
library(rlecuyer)
source("../../../waic_fxns.R")

##  Read in data
#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)

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


growD <- subset(growD_all, Species=="PASM" & year != 33)
clim_covs <- growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
clim_covs$sizepptLag <- clim_covs$pptLag*log(growD$percLagCover)
clim_covs$sizeppt1 <- clim_covs$ppt1*log(growD$percLagCover)
clim_covs$sizeppt2 <- clim_covs$ppt2*log(growD$percLagCover)
clim_covs$sizetemp1 <- clim_covs$TmeanSpr1*log(growD$percLagCover)
clim_covs$sizetemp2 <- clim_covs$TmeanSpr2*log(growD$percLagCover)
groups <- as.numeric(as.factor(growD$group))
G <- length(unique(growD$group))
Yrs <- length(unique(growD$year))
yid <- as.numeric(as.factor(growD$year))

hold_data <- subset(growD_all, Species=="PASM" & year == 33)
clim_covs_out <- hold_data[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]

### Initialize Regularization MCMC
datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=growD$percCover, X=log(growD$percLagCover),
                 C=clim_covs, G=G, gid=groups, sd_clim=0.1)
pars=c("a_mu")
mcmc_reg <- stan(file="qbm_reg_cv.stan", data=datalist, pars=pars, chains=0)

### Initialize OOS-CV MCMC
datalist <- list(N=nrow(growD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=growD$percCover, X=log(growD$percLagCover),
                 C=clim_covs, G=G, gid=groups, sd_clim=0.1,
                 y_holdout=hold_data$percCover, X_out=log(hold_data$percLagCover),
                 C_out=clim_covs_out, npreds=nrow(hold_data))
pars=c("a_mu")
mcmc_oos <- stan(file="qbm_oos_cv.stan", data=datalist, pars=pars, chains=0)


####
####  Use Parallel Regularization Fit with Full Data Sets 
####
####  Issue this command in shell before starting R: export OMP_NUM_THREADS=1 
####
grow_pasm <- subset(growD_all, Species=="PASM")
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
  clim_covs <- grow_pasm[,clim_vars]
  clim_covs$sizepptLag <- clim_covs$pptLag*log(grow_pasm$percLagCover)
  clim_covs$sizeppt1 <- clim_covs$ppt1*log(grow_pasm$percLagCover)
  clim_covs$sizeppt2 <- clim_covs$ppt2*log(grow_pasm$percLagCover)
  clim_covs$sizetemp1 <- clim_covs$TmeanSpr1*log(grow_pasm$percLagCover)
  clim_covs$sizetemp2 <- clim_covs$TmeanSpr2*log(grow_pasm$percLagCover)
  groups <- as.numeric(as.factor(grow_pasm$group))
  G <- length(unique(grow_pasm$group))
  Yrs <- length(unique(grow_pasm$year))
  yid <- as.numeric(as.factor(grow_pasm$year))
  sd.now <- sd_vec[i]
  datalist <- list(N=nrow(grow_pasm), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=grow_pasm$percCover, 
                   X=log(grow_pasm$percLagCover),
                   C=clim_covs, G=G, gid=groups, sd_clim=sd.now)
  pars <- c("b2")
  fit <- stan(fit = mcmc_reg, data=datalist,
              pars=pars, chains=2, iter = 2000, warmup = 1000)
  long <- ggs(fit)
  longagg <- ddply(ggs(fit), .(Parameter), summarise,
                   avg_value = mean(value))
  return(longagg$avg_value)
}

sfExport("sd_vec", "grow_pasm")
tmp.time=Sys.time()
beta.list=sfClusterApplySR(1:n.beta,reg.fcn,perUpdate=1)
time.1=Sys.time()-tmp.time
sfStop()
time.1

beta.mat=matrix(unlist(beta.list),ncol=ncol(clim_covs),byrow=TRUE)


####
####  Set up cross validation and regularization loop
####
####  Issue this command in shell before starting R: export OMP_NUM_THREADS=1 
####
grow_pasm <- subset(growD_all, Species=="PASM")
n.beta <- 24
sd_vec <- seq(.1,1.5,length.out = n.beta)
yrs.vec <- unique(subset(growD_all, Species=="PASM")$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)

clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")

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
  df_train <- subset(grow_pasm, year!=yr.out)
  df_hold <- subset(grow_pasm, year==yr.out) 
  clim_covs <- df_train[,clim_vars]
  clim_covs$sizepptLag <- clim_covs$pptLag*log(df_train$percLagCover)
  clim_covs$sizeppt1 <- clim_covs$ppt1*log(df_train$percLagCover)
  clim_covs$sizeppt2 <- clim_covs$ppt2*log(df_train$percLagCover)
  clim_covs$sizetemp1 <- clim_covs$TmeanSpr1*log(df_train$percLagCover)
  clim_covs$sizetemp2 <- clim_covs$TmeanSpr2*log(df_train$percLagCover)
  groups <- as.numeric(as.factor(df_train$group))
  G <- length(unique(df_train$group))
  Yrs <- length(unique(df_train$year))
  yid <- as.numeric(as.factor(df_train$year))
  clim_covs_out <- df_hold[,clim_vars]
  clim_covs_out$sizepptLag <- clim_covs_out$pptLag*log(df_hold$percLagCover)
  clim_covs_out$sizeppt1 <- clim_covs_out$ppt1*log(df_hold$percLagCover)
  clim_covs_out$sizeppt2 <- clim_covs_out$ppt2*log(df_hold$percLagCover)
  clim_covs_out$sizetemp1 <- clim_covs_out$TmeanSpr1*log(df_hold$percLagCover)
  clim_covs_out$sizetemp2 <- clim_covs_out$TmeanSpr2*log(df_hold$percLagCover)
  
  datalist <- list(N=nrow(df_train), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=df_train$percCover, 
                   X=log(df_train$percLagCover),
                   C=clim_covs, G=G, gid=groups, sd_clim=sd.now,
                   y_holdout=df_hold$percCover, X_out=log(df_hold$percLagCover),
                   C_out=clim_covs_out, npreds=nrow(df_hold))
  pars <- c("log_lik")
  fit <- stan(fit = mcmc_oos, data=datalist,
              pars=pars, chains=2, iter = 2000, warmup = 1000)
  waic_metrics <- waic(fit)
  lpd <- waic_metrics[["total"]]["elpd_loo"]
  return(lpd)
}

sfExport("sd_vec","cv.s2.grid","cv.fcn",
         "fold.idx.mat","yrs.vec","grow_pasm")
tmp.time=Sys.time()
score.list=sfClusterApplySR(1:n.grid,cv.fcn,perUpdate=1)
time.2=Sys.time()-tmp.time
sfStop()
time.2
score.cv.mat=matrix(unlist(score.list),n.beta,K)


####
####  Shrinkage Trajectories and CV Score
####
score.cv.vec=apply(score.cv.mat,1,sum)
sd.beta.opt <- sd_vec[which(score.cv.vec==max(score.cv.vec))]
png("cv_score.png", width = 6, height=8, units = "in", res=100)
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

