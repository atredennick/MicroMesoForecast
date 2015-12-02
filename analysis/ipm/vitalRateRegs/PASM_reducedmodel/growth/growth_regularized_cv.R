##  Script for Bayesian regularization by cross-validation
##  for growth regressions

##  Regularization implemented on the prior for the climate effects

rm(list=ls())

####
####  Load libraries -----------------------------------------------------------
####
library(rstan)
library(plyr)
library(reshape2)
library(matrixStats)
source("../waic_fxns.R")



####
####  Read in data -------------------------------------------------------------
####
sppList=sort(c("BOGR", "HECO", "PASM", "POSE"))
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
    sppD <- read.csv(paste("../../../../speciesData/", doSpp, "/edited/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../../speciesData/", doSpp, "/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

climD <- read.csv("../../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
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
growD_all <- merge(growD, crowd, by=c("species", "X"))

growD <- subset(growD_all, species=="BOGR")
clim_invest <- ddply(growD, .(year), summarise,
                     pptLag = mean(pptLag),
                     ppt1 = mean(ppt1),
                     ppt2 = mean(ppt2),
                     TmeanSpr1 = mean(TmeanSpr1),
                     TmeanSpr2 = mean(TmeanSpr2))
cor(clim_invest[,-1])
det(cor(clim_invest[,-1]))
X.cor.test <- scale(clim_invest[,-1])
CN <- kappa(X.cor.test)


####
####  Create datalist with train and holdout data ------------------------------
####
set.seed(12345)
"%w/o%" <- function(x, y) x[!x %in% y] # x without y
growD <- subset(growD_all, species=="BOGR")
##  Train on 75% of years, validate on 25%
yrs <- unique(growD$year)
train_years <- sample(c(1:length(yrs)), size = 0.99*length(yrs))
grow_train <- growD[growD$year %in% yrs[train_years], ]
grow_hold <- growD[growD$year %w/o% yrs[train_years], ]

##  Climate covariates
clim_train <- as.data.frame(grow_train[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")])
# clim_train$sizepptLag <- clim_train$pptLag*log(grow_train$area.t0)
# clim_train$sizeppt1 <- clim_train$ppt1*log(grow_train$area.t0)
# clim_train$sizeppt2 <- clim_train$ppt2*log(grow_train$area.t0)
# clim_train$sizeTmeanSpr1 <- clim_train$TmeanSpr1*log(grow_train$area.t0)
# clim_train$sizeTmenaSpr2 <- clim_train$TmeanSpr2*log(grow_train$area.t0)

clim_hold <- as.data.frame(grow_hold[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")])
# clim_hold$sizepptLag <- clim_hold$pptLag*log(grow_hold$area.t0)
# clim_hold$sizeppt1 <- clim_hold$ppt1*log(grow_hold$area.t0)
# clim_hold$sizeppt2 <- clim_hold$ppt2*log(grow_hold$area.t0)
# clim_hold$sizeTmeanSpr1 <- clim_hold$TmeanSpr1*log(grow_hold$area.t0)
# clim_hold$sizeTmeanSpr2 <- clim_hold$TmeanSpr2*log(grow_hold$area.t0)

##  Group and year vectors; crowding
group_train <- as.numeric(grow_train$Group)
ngrp_train <- length(unique(grow_train$Group))
gid_train <- as.numeric(grow_train$Group)
nyrs_train <- length(unique(grow_train$year))
yid_train <- as.numeric(as.factor(grow_train$year))
W_train <- cbind(grow_train$W, grow_train$W*log(grow_train$area.t0))

group_hold <- as.numeric(grow_hold$Group)
ngrp_hold <- length(unique(grow_hold$Group))
gid_hold <- as.numeric(grow_hold$Group)
nyrs_hold <- length(unique(grow_hold$year))
yid_hold <- as.numeric(as.factor(grow_hold$year))
W_hold <- cbind(grow_hold$W, grow_hold$W*log(grow_hold$area.t0))
Wish <- diag(length(clim_train))

datalist <- list(N=nrow(grow_train), Yrs=nyrs_train, yid=yid_train,
                 Covs=length(clim_train), Y=log(grow_train$area.t1), 
                 X=log(grow_train$area.t0),
                 C=clim_train, W=W_train, G=ngrp_train, gid=gid_train,
                 npreds=nrow(grow_hold), y_holdout=log(grow_hold$area.t1),
                 X_out=log(grow_hold$area.t0), C_out=clim_hold, W_out=W_hold,
                 gid_out=gid_hold, tau_betas=0.01)


####
####  Compile and fit via STAN -------------------------------------------------
####
## Set reasonable initial values for three chains
inits <- list()
inits[[1]] <- list(a_mu=0, a=rep(0,nyrs_train), b1_mu=0.01,
                   b1=rep(0.01,nyrs_train), gint=rep(0,ngrp_train), 
                   w=c(0,0), sig_b1=0.5, sig_a=0.5, 
                   tau=0.5, tauSize=0.5, sig_G=0.5, 
                   b2=rep(0,length(clim_train)))
pars <- c("b2")
fitted <- stan(file = "growth.stan", data=datalist, init = list(inits[[1]]),
               pars=pars, chains=1, iter=200, warmup = 100)
plot(fitted)
# waic_metrics <- waic(fitted)

