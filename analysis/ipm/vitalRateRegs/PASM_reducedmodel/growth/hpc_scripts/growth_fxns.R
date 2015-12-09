##
##  R Functions for Parallelized Stan Fits for survival logistic model (for IPM)
##
##
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 12-6-2015
##

####
####  Regularization Without Cross Validation (collects betas)
####
reg.fcn <- function(i){
  library(rstan)
  library(ggmcmc)
  library(plyr)
  clim_covs <- growD[,clim_vars]
  groups <- as.numeric(growD$Group)
  G <- length(unique(growD$Group))
  gid <- as.numeric(growD$Group)
  nyrs <- length(unique(growD$year))
  W <- cbind(growD$W, growD$W*log(growD$area.t0))
  yid <- as.numeric(as.factor(growD$year))
  sd_now <- sd_vec[i]
  datalist <- list(N=nrow(growD), Yrs=nyrs, yid=yid,
                   Covs=ncol(clim_covs), Y=log(growD$area.t1), 
                   X=log(growD$area.t0), C=clim_covs, W=W, G=G, gid=gid,
                   tau_beta=sd_now)
  
  pars <- c("b2")
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,nyrs), b1_mu=0.01,
                     b1=rep(0.01,nyrs), gint=rep(0,G), 
                     w=c(0,0), sig_b1=0.5, sig_a=0.5, 
                     tau=0.5, tauSize=0.5, sig_G=0.5, 
                     b2=rep(0,ncol(clim_covs)))
  inits[[2]] <- list(a_mu=0.5, a=rep(0.5,nyrs), b1_mu=0.1,
                     b1=rep(0.1,nyrs), gint=rep(0.4,G), 
                     w=c(0.1,0.1), sig_b1=0.2, sig_a=0.2, 
                     tau=0.2, tauSize=0.2, sig_G=0.2, 
                     b2=rep(0.5,ncol(clim_covs)))
  inits[[3]] <- list(a_mu=1, a=rep(1,nyrs), b1_mu=0.2,
                     b1=rep(0.2,nyrs), gint=rep(0.2,G), 
                     w=c(0.5,0.5), sig_b1=0.01, sig_a=0.01, 
                     tau=0.05, tauSize=0.05, sig_G=0.05, 
                     b2=rep(-0.5,ncol(clim_covs)))
  fit <- stan(fit = mcmc_reg, data=datalist, init=inits,
              pars=pars, chains=3, iter = 2000, warmup = 1000)
  long <- ggs(fit)
  longagg <- ddply(ggs(fit), .(Parameter), summarise,
                   avg_value = mean(value))
  return(longagg$avg_value)
}




####
####  Regularization and Cross-Validation (collects lppd)
####
cv.fcn <- function(i){
  library(rstan)
  library(matrixStats)
  k <- cv.s2.grid[i,2]
  fold.idx <- fold.idx.mat[,k] 
  yr.out <- yrs.vec[fold.idx]
  sd.now <- sd_vec[cv.s2.grid[i,1]]
  df_train <- subset(growD, year!=yr.out)
  df_hold <- subset(growD, year==yr.out) 
  clim_covs <- df_train[,clim_vars]
  groups <- as.numeric(df_train$Group)
  G <- length(unique(df_train$Group))
  nyrs <- length(unique(df_train$year))
  W <- cbind(df_train$W, df_train$W*log(df_train$area.t0))
  yid <- as.numeric(as.factor(df_train$year))
  
  clim_covs_oos <- df_hold[,clim_vars]
  W_oos <- cbind(df_hold$W, df_hold$W*log(df_hold$area.t0))
  
  datalist <- list(N=nrow(df_train), Yrs=nyrs, yid=yid,
                   Covs=ncol(clim_covs), Y=log(df_train$area.t1), X=log(df_train$area.t0),
                   C=clim_covs, W=W, G=G, gid=groups, tau_beta=sd.now,
                   npreds=nrow(df_hold), y_holdout=log(df_hold$area.t1), Xhold=log(df_hold$area.t0),
                   Chold=clim_covs_oos, Whold=W_oos)
  pars <- c("log_lik")
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,nyrs), b1_mu=0.01,
                     b1=rep(0.01,nyrs), gint=rep(0,G), 
                     w=c(0,0), sig_b1=0.5, sig_a=0.5, 
                     tau=0.5, tauSize=0.5, sig_G=0.5, 
                     b2=rep(0,ncol(clim_covs)))
  inits[[2]] <- list(a_mu=0.5, a=rep(0.5,nyrs), b1_mu=0.1,
                     b1=rep(0.1,nyrs), gint=rep(0.4,G), 
                     w=c(0.1,0.1), sig_b1=0.2, sig_a=0.2, 
                     tau=0.2, tauSize=0.2, sig_G=0.2, 
                     b2=rep(0.5,ncol(clim_covs)))
  inits[[3]] <- list(a_mu=1, a=rep(1,nyrs), b1_mu=0.2,
                     b1=rep(0.2,nyrs), gint=rep(0.2,G), 
                     w=c(0.5,0.5), sig_b1=0.01, sig_a=0.01, 
                     tau=0.05, tauSize=0.05, sig_G=0.05, 
                     b2=rep(-0.5,ncol(clim_covs)))
  fit <- stan(fit = mcmc_oos, data=datalist, init=inits,
              pars=pars, chains=3, iter = 2000, warmup = 1000)
  waic_metrics <- waic(fit)
  lpd <- waic_metrics[["total"]]["elpd_loo"]
  return(lpd)
}

