##  R Script to Aggregate Simulation Results from One-Step-Ahead
##  Forecasts and to Make Table of Results for Manuscript
##
##  Author: Andrew Tredennick
##  Last update: 4-13-2016



####
####  LIBRARIES
####
library(plyr)
library(reshape2)
library(xtable)



####
####  PRELIMINARIES
####
path2ipm <- "../analysis/ipm/simulations/results/one_step_validation/"
path2qbm <- "../analysis/quadBM/simulations/results/one_step_validation/"
path2save <- "../manuscript/components/"



####
####  IPM 1-STEP FORECAST RESULTS
####
# Read in data
ipm_onestep <- readRDS(paste0(path2ipm,"ipm_loyo_forecasts_combined.RDS"))
# Remove predictions >1 for fair comparison to [0,1] truncated quad model
# resets <- which(ipm_onestep[,"cover.t1"]>1)
# ipm_onestep[resets, "cover.t1"] <- 1
# Average predictions over quad-year reps
avg_ipm <- ddply(ipm_onestep, .(species, quad, t1), summarise,
                 avg_prediction = median(cover.t1),
                 observation = mean(cover.t0),
                 quant_dist = quantile(cover.t1*100, probs = 0.95) - quantile(cover.t1, probs = 0.05))
# Calculate error
avg_ipm$error <- avg_ipm$observation*100 - avg_ipm$avg_prediction*100
# Aggregate over species
ipmstats <- ddply(avg_ipm, .(species), summarise,
                  mae = mean(abs(error)),
                  mean_cover = mean(observation*100),
                  quant_dist = mean(quant_dist))
ipmstats$model <- "IPM"



####
####  QBM 1-STEP FORECAST RESULTS
####
# Read in data
qbm_onestep <- readRDS(paste0(path2qbm, "qbm_one-step_forecasts_combined.RDS"))
# Average predictions over quad-year reps
avg_qbm <- ddply(qbm_onestep, .(species, quad, year), summarise,
                 avg_prediction = median(cover.t1),
                 observation = mean(cover.t0),
                 quant_dist = quantile(cover.t1*100, probs = 0.95) - quantile(cover.t1, probs = 0.05))
# Calculate error
avg_qbm$error <- avg_qbm$observation*100 - avg_qbm$avg_prediction*100
# Aggregate over species
qbmstats <- ddply(avg_qbm, .(species), summarise,
                  mae = mean(abs(error)),
                  mean_cover = mean(observation*100),
                  quant_dist = mean(quant_dist))
qbmstats$model <- "QBM"




####
####  CALCULATE P-VALUES FROM T-TESTS ON MODEL ERRORS
####  [SEE YE et al. 2015 (PNAS) FOR EXAMPLE USING MODEL ACCURACIES]
####
mae_list <- list()
mae_pvals <- list()
species <- unique(qbmstats$species)
for(spp in species){
  tmp1 <- subset(avg_ipm, species==spp)
  tmp2 <- subset(avg_qbm, species==spp)
  err1 <- as.data.frame(abs(tmp1$observation - tmp1$avg_prediction))
  err2 <- as.data.frame(abs(tmp2$observation - tmp2$avg_prediction))
  colnames(err1) <- "resid"
  err1$model <- "ipm"
  colnames(err2) <- "resid"
  err2$model <- "qbm"
  errd <-  rbind(err1, err2)
  mae_ttest <- with(errd, t.test(resid~model, alternative = "less"))
  mae_list[[spp]] <- mae_ttest
  mae_pvals[[spp]] <- mae_ttest[["p.value"]]
}

mae_ps <- melt(mae_pvals)



####
####  MAKE LATEX TABLE
####
stats_table <- rbind(ipmstats, qbmstats)
stats_table <- stats_table[ order(stats_table[,"species"]), ]
stats_table <- stats_table[,c("species", "model", "mae",
                              "quant_dist","mean_cover")]

table_caption <- "Accuracy (mean absolute error, MAE) and precision (90\\% Distance) 
                  of out of sample predictions. Forecasts were made without random 
                  year effects; only climate covariates could explain year-to-year 
                  variation. 90\\% Distance refers to the average distance between the 
                  upper and lower 90th percentiles of the 100 predicted values for 
                  each quadrat-year combination."

colnames(stats_table) <- c("Species", "Model", "MAE", "90% Distance", "Mean Obs. Cover")
print.xtable(xtable(stats_table, caption = table_caption), type="latex", comment=FALSE,
             include.rownames=FALSE, caption.placement="top", file = paste0(path2save,"error_table.tex"))
