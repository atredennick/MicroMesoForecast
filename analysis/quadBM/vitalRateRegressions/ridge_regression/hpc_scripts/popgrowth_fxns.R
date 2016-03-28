## Function for Parallel OOS-CV for Population Growth Model

cv.fcn <- function(i){
  library(rstan)
  library(matrixStats)
  k <- cv.s2.grid[i,2]
  fold.idx <- fold.idx.mat[,k] 
  yr.out <- yrs.vec[fold.idx]
  sd.now <- sd_vec[cv.s2.grid[i,1]]
  df_train <- subset(grow_data, year!=yr.out)
  df_hold <- subset(grow_data, year==yr.out) 
  ##  Create and scale interaction covariates
  df_train$ppt1TmeanSpr1 <- df_train$ppt1*df_train$TmeanSpr1
  df_train$ppt2TmeanSpr2 <- df_train$ppt2*df_train$TmeanSpr2
  df_train$sizepptLag <- df_train$pptLag*log(df_train$percLagCover)
  df_train$sizeppt1 <- df_train$ppt1*log(df_train$percLagCover)
  df_train$sizeppt2 <- df_train$ppt2*log(df_train$percLagCover)
  df_train$sizeTmeanSpr1 <- df_train$TmeanSpr1*log(df_train$percLagCover)
  df_train$sizeTmeanSpr2 <- df_train$TmeanSpr2*log(df_train$percLagCover)
  clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                     "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
  clim_covs <- df_train[,clim_vars_all]
  # Get scalers for climate covariates from training data
  clim_means <- colMeans(clim_covs)
  clim_sds <- apply(clim_covs, 2, FUN = sd)
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  groups <- as.numeric(as.factor(df_train$group))
  G <- length(unique(df_train$group))
  Yrs <- length(unique(df_train$year))
  yid <- as.numeric(as.factor(df_train$year))
  
  df_hold$ppt1TmeanSpr1 <- df_hold$ppt1*df_hold$TmeanSpr1
  df_hold$ppt2TmeanSpr2 <- df_hold$ppt2*df_hold$TmeanSpr2
  df_hold$sizepptLag <- df_hold$pptLag*log(df_hold$percLagCover)
  df_hold$sizeppt1 <- df_hold$ppt1*log(df_hold$percLagCover)
  df_hold$sizeppt2 <- df_hold$ppt2*log(df_hold$percLagCover)
  df_hold$sizeTmeanSpr1 <- df_hold$TmeanSpr1*log(df_hold$percLagCover)
  df_hold$sizeTmeanSpr2 <- df_hold$TmeanSpr2*log(df_hold$percLagCover)
  clim_covs_oos <- df_hold[,clim_vars_all]
  for(j in 1:ncol(clim_covs_oos)){
    clim_covs_oos[,j] <- (clim_covs_oos[,j] - clim_means[j])/clim_sds[j]
  }
  
  datalist <- list(N=nrow(df_train), Yrs=Yrs, yid=yid,
                   Covs=ncol(clim_covs), Y=df_train$percCover, 
                   X=log(df_train$percLagCover),
                   C=clim_covs, G=G, gid=groups, sd_clim=sd.now,
                   y_holdout=df_hold$percCover, X_out=log(df_hold$percLagCover),
                   C_out=clim_covs_oos, npreds=nrow(df_hold))
  pars <- c("log_lik")
  fit <- stan(fit = mcmc_oos, data=datalist,
              pars=pars, chains=3, iter = 2000, warmup = 1000)
  waic_metrics <- waic(fit)
  lpd <- waic_metrics[["total"]]["elpd_loo"]
  return(lpd)
}

