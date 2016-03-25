##
##  R Functions for Parallelized Stan Fits for survival logistic model (for IPM)
##
##
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 1-20-2016
##

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
  df_train <- subset(recD, year!=yr.out)
  df_hold <- subset(recD, year==yr.out) 
  ##  Create and scale interaction covariates
  df_train$ppt1TmeanSpr1 <- df_train$ppt1*df_train$TmeanSpr1
  df_train$ppt2TmeanSpr2 <- df_train$ppt2*df_train$TmeanSpr2
  clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
  clim_covs <- df_train[,clim_vars_all]
  # Get scalers for climate covariates from training data
  clim_means <- colMeans(clim_covs)
  clim_sds <- apply(clim_covs, 2, FUN = sd)
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  groups <- as.numeric(df_train$group)
  G <- length(unique(df_train$group))
  Yrs <- length(unique(df_train$year))
  yid <- as.numeric(as.factor(df_train$year))
  sd.now <- sd_vec[cv.s2.grid[i,1]]
  
  df_hold$ppt1TmeanSpr1 <- df_hold$ppt1*df_hold$TmeanSpr1
  df_hold$ppt2TmeanSpr2 <- df_hold$ppt2*df_hold$TmeanSpr2
  clim_covs_oos <- df_hold[,clim_vars_all]
  for(j in 1:ncol(clim_covs_oos)){
    clim_covs_oos[,j] <- (clim_covs_oos[,j] - clim_means[j])/clim_sds[j]
  }
  groups_out <- as.numeric(df_hold$group)
  Gout <- length(unique(df_hold$group))
  npreds <- nrow(df_hold)
  y_holdout <- df_hold$recruits
  parents1_out <- df_hold$parents1
  parents2_out <- df_hold$parents2
  
  datalist <- list(N=nrow(df_train), Yrs=Yrs, yid=yid,
                   Covs=ncol(clim_covs), Y=df_train$recruits, C=clim_covs, 
                   parents1=df_train$parents1, parents2=df_train$parents2,
                   G=G, gid=groups, tau=sd.now,
                   npreds=npreds, y_holdout=y_holdout, C_out=clim_covs_oos,
                   parents1_out=parents1_out, parents2_out=parents2_out,
                   gid_out=groups_out, Gout=Gout)
  pars <- c("log_lik")
  inits <- list()
  inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                  gint=rep(0,G), sig_G=1, u=0.4,
                  dd=-1,theta=1,b2=rep(0,ncol(clim_covs)))
  inits[[2]]=list(a=rep(0.5,Yrs), a_mu=0.2, sig_a=10,
                  gint=rep(0,G), sig_G=0.1,  u=0.7,
                  dd=-0.05,theta=1.5, b2=rep(0.5,ncol(clim_covs))) 
  inits[[3]]=list(a=rep(1,Yrs), a_mu=0.5, sig_a=5,
                  gint=rep(-0.1,G), sig_G=0.5,  u=0.5,
                  dd=-0.1,theta=1.2, b2=rep(-0.5,ncol(clim_covs))) 
  fit <- stan(fit = mcmc_oos, data=datalist, init=inits,
              pars=pars, chains=3, iter = 2000, warmup = 1000)
  waic_metrics <- waic(fit)
  lpd <- waic_metrics[["total"]]["lpd"]
  return(lpd)
}
