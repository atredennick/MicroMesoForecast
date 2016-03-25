##
##  R Script for Recruitment Ricker Model (for IPM)
##
##  Includes Model Selection via Bayesian Regularization
##
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 12-15-2015
##
##  NOTICE: 
##  As is, this script will take ~20 minutes to run.
##  Production runs (more MCMC iterations) required use of 
##  Utah State University's High Performance Computing facility.
##  See hpc_scripts/ directory for fully parallelized versions.
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
sppList=sort(c("BOGR","HECO","PASM","POSE"))
do_species <- 4

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../../speciesData/", doSpp, "/edited/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../../speciesData/", doSpp, "/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }
  
  sppD$Group <- as.factor(substr(sppD$quad,1,1)) #add by Chengjin
  sppD <- sppD[,c("quad","year","NRquad","totParea","Group")]
  names(sppD)[3] <- paste("R.",sppList[spp],sep="")
  names(sppD)[4] <- paste("cov.",sppList[spp],sep="")
  if(spp==1){
    D <- sppD
  }else{
    D <- merge(D,sppD,all=T)
  }
}
D[is.na(D)]=0  # replace missing values 

# calculate mean cover by group and year
tmpD <- D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
tmpD <- aggregate(tmpD[,4:NCOL(tmpD)],
                  by=list("year"=tmpD$year,"Group"=tmpD$Group),FUN=mean)
names(tmpD)[3:NCOL(tmpD)] <- paste("Gcov.",sppList,sep="")
D <- merge(D,tmpD,all.x=T)

parents1 <- as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100 ##convert from absolute cover to [1,100] range
parents2 <- as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100

##for species 1
tmp1L <- which(parents1[,1]==0) ##lcoal
tmp1G <- which(parents2[,1]==0) ##Group
tmp1 <- intersect(tmp1L,tmp1G)
##for species 2
tmp2L <- which(parents1[,2]==0)
tmp2G <- which(parents2[,2]==0)
tmp2 <- intersect(tmp2L,tmp2G)
##for species 3
tmp3L <- which(parents1[,3]==0)
tmp3G <- which(parents2[,3]==0)
tmp3 <- intersect(tmp3L,tmp3G)
##for species 4
tmp4L <- which(parents1[,4]==0)
tmp4G <- which(parents2[,4]==0)
tmp4 <- intersect(tmp4L,tmp4G)

tmp <- unique(c(tmp1,tmp2,tmp3,tmp4))

if(length(tmp)>0){
  parents1 <- parents1[-tmp,] ##remove them
  parents2 <-parents2[-tmp,] ##remove them
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])[-tmp,] ##remove them  
  year <- D$year[-tmp] ##remove them
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]-length(tmp) ##reduce
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group))[-tmp] ##remove them ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
} else {
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])
  year <- D$year
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group)) ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
}

tmpY <- melt(y)
tmpP1 <- melt(parents1)
tmpP2 <- melt(parents2)
allD <- data.frame(species=tmpY$Var2,
                   year=rep(year,4),
                   group=rep(Group,4),
                   recruits=tmpY$value,
                   parents1=tmpP1$value,
                   parents2=tmpP2$value)

##  Add in climate data
climD <- read.csv("../../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
allD <- merge(allD,climD)



####
####  Collate Data for Stan Model Initialization
####
## Initialize shrinkage model with full data set
recD <- subset(allD, species==paste("R.",sppList[do_species],sep=""))
clim_covs <- recD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))
yid <- as.numeric(as.factor(recD$year))
datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=recD$recruits, C=clim_covs, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups, tau=1)
pars <- c("b2")
mcmc_reg <- stan(file="recruitment_reg_cv.stan", data=datalist, 
                 pars=pars, chains=0)

##  Initialize out-of-sample cross validation model
holdyear <- unique(allD$year)[1]
recD <- subset(allD, species==paste0("R.",sppList[do_species]) & year!=holdyear)
clim_covs <- recD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))
yid <- as.numeric(as.factor(recD$year))

holdD <- subset(allD, species==paste0("R.",sppList[do_species]) & year==holdyear)
clim_holds <- holdD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups_out <- as.numeric(holdD$group)
Gout <- length(unique(holdD$group))
npreds <- nrow(holdD)
y_holdout <- holdD$recruits
parents1_out <- holdD$parents1
parents2_out <- holdD$parents2

datalist <- list(N=nrow(recD), Yrs=Yrs, yid=(recD$year-31),
                 Covs=length(clim_covs), Y=recD$recruits, C=clim_covs, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups, tau=1,
                 npreds=npreds, y_holdout=y_holdout, C_out=clim_holds,
                 parents1_out=parents1_out, parents2_out=parents2_out,
                 gid_out=groups_out, Gout=Gout)

pars=c("b2")
mcmc_oos <- stan(file = "recruitment_oos_cv.stan", data=datalist,
                     pars=pars, chains=0)



####
####  Use Parallel Regularization Fit with Full Data Sets 
####
####  Issue this command in shell before starting R: export OMP_NUM_THREADS=1 
####
recD <- subset(allD, species==paste("R.",sppList[do_species],sep=""))
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
  clim_covs <- recD[,clim_vars]
  groups <- as.numeric(recD$group)
  G <- length(unique(recD$group))
  Yrs <- length(unique(recD$year))
  yid <- as.numeric(as.factor(recD$year))
  datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=recD$recruits, C=clim_covs, 
                   parents1=recD$parents1, parents2=recD$parents2,
                   G=G, gid=groups, tau=sd_vec[i])
  pars <- c("b2")
  inits <- list()
  inits=list()
  inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                  gint=rep(0,G), sig_G=1, u=0.4,
                  dd=-1,theta=1,b2=rep(0,length(clim_covs)))
  fit <- stan(fit = mcmc_reg, data=datalist, init=list(inits[[1]]),
              pars=pars, chains=1, iter = 1000, warmup = 500)
  long <- ggs(fit)
  longagg <- ddply(ggs(fit), .(Parameter), summarise,
                   avg_value = mean(value))
  return(longagg$avg_value)
}

clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
sfExport("sd_vec", "recD", "clim_vars")
tmp.time <- Sys.time()
beta.list <- sfClusterApplySR(1:n.beta,reg.fcn,perUpdate=1)
time.1 <- Sys.time()-tmp.time
sfStop()
time.1

beta.mat <- matrix(unlist(beta.list),ncol=ncol(clim_covs),byrow=TRUE)



####
####  Fit Each Recruitment Model using Parallelized C-V 
####
####  Issue this command in shell before starting R: export OMP_NUM_THREADS=1 
####
recD <- subset(allD, species==paste("R.",sppList[do_species],sep=""))
n.beta <- 24
sd_vec <- seq(.1,1.5,length.out = n.beta)
yrs.vec <- unique(recD$year)
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
  df_train <- subset(recD, year!=yr.out)
  df_hold <- subset(recD, year==yr.out) 
  clim_covs <- df_train[,clim_vars]
  groups <- as.numeric(df_train$group)
  G <- length(unique(df_train$group))
  Yrs <- length(unique(df_train$year))
  yid <- as.numeric(as.factor(df_train$year))
  sd.now <- sd_vec[cv.s2.grid[i,1]]
  
  groups_out <- as.numeric(holdD$group)
  Gout <- length(unique(holdD$group))
  npreds <- nrow(holdD)
  y_holdout <- holdD$recruits
  parents1_out <- holdD$parents1
  parents2_out <- holdD$parents2
  
  datalist <- list(N=nrow(df_train), Yrs=Yrs, yid=yid,
                   Covs=length(clim_covs), Y=df_train$recruits, C=clim_covs, 
                   parents1=df_train$parents1, parents2=df_train$parents2,
                   G=G, gid=groups, tau=sd.now,
                   npreds=npreds, y_holdout=y_holdout, C_out=clim_holds,
                   parents1_out=parents1_out, parents2_out=parents2_out,
                   gid_out=groups_out, Gout=Gout)
  pars <- c("log_lik")
  inits <- list()
  inits=list()
  inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                  gint=rep(0,G), sig_G=1, u=0.4,
                  dd=-1,theta=1,b2=rep(0,length(clim_covs)))
  fit <- stan(fit = mcmc_oos, data=datalist, init=list(inits[[1]]),
              pars=pars, chains=1, iter = 1000, warmup = 500)
  waic_metrics <- waic(fit)
  lpd <- waic_metrics[["total"]]["elpd_loo"]
  return(lpd)
}

sfExport("sd_vec","cv.s2.grid","cv.fcn",
         "fold.idx.mat","yrs.vec","recD")
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
# png("cv_score_surv.png", width = 6, height=8, units = "in", res=100)
par(mfrow=c(2,1), mar=c(1,4.1,4.1,2.1))
plot(sd_vec^2,beta.mat[,1],type="l",lwd=3,ylab=bquote(beta[mean]),
     xlab=bquote(sigma[beta]^2), ylim=c(-1,1))
for(i in 2:ncol(beta.mat)){
  lines(sd_vec^2,beta.mat[,i],type="l",lwd=3,col=i)
}
abline(h=0, lty=2)
abline(v=sd.beta.opt^2,col=8,lwd=2)
legend("bottomright",legend = names(clim_covs), col = c(1:ncol(beta.mat)), 
       lty = 1, bty="n", ncol=2, cex = 0.5)
par(mar=c(5,4.1,2,2.1))
plot(sd_vec^2,score.cv.vec,type="l",lwd=3,ylab="Log Predictive Score (lppd)",xlab=bquote(sigma[beta]^2))
abline(v=sd.beta.opt^2,col=8,lwd=2)
# dev.off()



####
####  Fit Optimal Predictive Model Using Full Dataset
####
recD <- subset(allD, species==paste("R.",sppList[do_species],sep=""))
clim_covs <- recD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))
yid <- as.numeric(as.factor(recD$year))
datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                 Covs=length(clim_covs), Y=recD$recruits, C=clim_covs, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups, tau=sd.beta.opt)
pars <- c("b2")
inits <- list()
inits=list()
inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                gint=rep(0,G), sig_G=1, u=0.4,
                dd=-1,theta=1,b2=rep(0,length(clim_covs)))
fit <- stan(fit = mcmc_reg, data=datalist, init=list(inits[[1]]),
            pars=pars, chains=1, iter = 2000, warmup = 1000)


