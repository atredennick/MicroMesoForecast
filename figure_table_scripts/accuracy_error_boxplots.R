##  R Script to Aggregate Simulation Results from One-Step-Ahead
##  Forecasts and to Make Figure 2 for Manuscript
##  Also Outputs a Table of Significance Tests for rho Comparisons
##
##  Author: Andrew Tredennick
##  Last update: 4-13-2016



####
####  LIBRARIES
####
library(plyr)
library(reshape2)
library(xtable)
library(ggplot2)



####
####  PRELIMINARIES
####
path2ipm <- "../analysis/ipm/simulations/results/one_step_validation/"
path2qbm <- "../analysis/quadBM/simulations/results/one_step_validation/"
path2save <- "../manuscript/components/"

##  Source statistical functions for comparing correlations
source("rho_comp_functions.R")


####
####  IPM 1-STEP FORECAST RESULTS
####
# Read in data
ipm_onestep <- readRDS(paste0(path2ipm,"ipm_loyo_forecasts_combined.RDS"))
torms <- which(is.na(ipm_onestep$cover.t1)==TRUE) # get rid of NAs; only a couple
ipm_onestep <- ipm_onestep[-torms,]

# Average predictions over quad-year reps
avg_ipm <- ddply(ipm_onestep, .(species, quad, t1), summarise,
                 avg_prediction = median(cover.t1),
                 observation = mean(obs.cover.t1),
                 quant_dist = quantile(cover.t1*100, probs = 0.95) - quantile(cover.t1, probs = 0.05))
# Calculate error
avg_ipm$error <- avg_ipm$observation*100 - avg_ipm$avg_prediction*100

# Aggregate by species and t1
ipmstats <- ddply(avg_ipm, .(species, t1), summarise,
                  mae = mean(abs(error)),
                  rho = cor(observation, avg_prediction, method = "p"),
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
                 avg_prediction = median(pred_cover.t1),
                 observation = mean(obs_cover.t1),
                 quant_dist = quantile(pred_cover.t1*100, probs = 0.95) - quantile(pred_cover.t1, probs = 0.05))
# Calculate error
avg_qbm$error <- avg_qbm$observation*100 - avg_qbm$avg_prediction*100

# Aggregate by species
qbmstats <- ddply(avg_qbm, .(species, year), summarise,
                  mae = mean(abs(error)),
                  rho = cor(observation, avg_prediction, method = "p"),
                  mean_cover = mean(observation*100),
                  quant_dist = mean(quant_dist))
qbmstats$model <- "QBM"

# REPLACE mean QBM cover with cover from IPM because it is the correct value
qbmstats$mean_cover <- ipmstats$mean_cover


testq <- subset(qbm_onestep, quad=="B1"&sim==1)
testi <- subset(ipm_onestep, quad=="B1"&rep==1)
testi$year <- testi$t0+1900
test_df <- merge(testq, testi)
if(mean(with(test_df, obs_cover.t1-obs.cover.t1))!=0) { stop("observed values not matching up") }



####
####  CALCULATE STATISTICS FROM NO-WEATHER SIMULATIONS (DENS. DEP. ONLY)
####
# IPM
path2ipm2 <- "../analysis/ipm/simulations/results/one_step_validation_densdep_only/"
# Read in data
ipm_onestep_dd <- readRDS(paste0(path2ipm2,"ipm_loyo_forecasts_noweather_combined.RDS"))
torms <- which(is.na(ipm_onestep_dd$cover.t1)==TRUE) # get rid of NAs; only a couple
ipm_onestep_dd <- ipm_onestep_dd[-torms,]
# Average predictions over quad-year reps
avg_ipm_dd <- ddply(ipm_onestep_dd, .(species, quad, t1), summarise,
                    avg_prediction = median(cover.t1),
                    observation = mean(obs.cover.t1),
                    quant_dist = quantile(cover.t1*100, probs = 0.95) - quantile(cover.t1, probs = 0.05))
# Calculate error
avg_ipm_dd$error <- avg_ipm_dd$observation*100 - avg_ipm_dd$avg_prediction*100
# Aggregate over species
ipmstats_dd <- ddply(avg_ipm_dd, .(species, t1), summarise,
                     mae = mean(abs(error)),
                     rho = cor(observation, avg_prediction, method = "p"),
                     mean_cover = mean(observation*100),
                     quant_dist = mean(quant_dist))
ipmstats_dd$model <- "IPM"

#QBM
path2qbm2 <- "../analysis/quadBM/simulations/results/one_step_validation_densdep_only/"
qbm_onestep_dd <- readRDS(paste0(path2qbm2, "qbm_one-step_forecasts_combined_dendeponly.RDS"))
# Average predictions over quad-year reps
avg_qbm_dd <- ddply(qbm_onestep_dd, .(species, quad, year), summarise,
                    avg_prediction = median(pred_cover.t1),
                    observation = mean(obs_cover.t1),
                    quant_dist = quantile(pred_cover.t1*100, probs = 0.95) - quantile(pred_cover.t1, probs = 0.05))
# Calculate error
avg_qbm_dd$error <- avg_qbm_dd$observation*100 - avg_qbm_dd$avg_prediction*100
# Aggregate over species
qbmstats_dd <- ddply(avg_qbm_dd, .(species, year), summarise,
                     mae = mean(abs(error)),
                     rho = cor(observation, avg_prediction, method = "p"),
                     mean_cover = mean(observation*100),
                     quant_dist = mean(quant_dist))
qbmstats_dd$model <- "QBM"

# REPLACE mean QBM cover with cover from IPM because it is the correct value
qbmstats_dd$mean_cover <- ipmstats_dd$mean_cover



####
####  MAKE FIGURE 2 -- ACCURACY BOXPLOTS
####
ipmstats$weather <- "yes"
ipmstats_dd$weather <- "no"
qbmstats$weather <- "yes"
qbmstats_dd$weather <- "no"
colnames(qbmstats_dd)[2] <- "t1"
colnames(qbmstats)[2] <- "t1"
plot_dat <- rbind(ipmstats,ipmstats_dd, qbmstats,qbmstats_dd)
plot_dat$modelweather <- paste0(plot_dat$model,plot_dat$weather)

library(ggthemes)
library(gridExtra)

ggplot(plot_dat, aes(x=species, y=rho, fill=modelweather))+
  stat_boxplot(geom ='errorbar', aes(color=modelweather)) +
  geom_boxplot(outlier.size = 0, color="white", coef = 0)+
  # geom_bar(stat="identity", position="dodge", width=0.8, color="black")+
  scale_fill_manual(values=c("coral", "darkred", "dodgerblue", "darkblue"),
                    labels = c("IPM - no climate", "IPM - climate", "QBM - no climate", "QBM - climate"))+
  scale_color_manual(values=c("coral", "darkred", "dodgerblue", "darkblue"),
                    labels = c("IPM - no climate", "IPM - climate", "QBM - no climate", "QBM - climate"))+
  ylab(expression(paste("Forecast Accuracy (", rho, ")")))+
  xlab("Species")+
  scale_y_continuous(limits=c(0,1))+
  theme_few()+
  theme(legend.position=c(0.87,0.15),
        legend.background = element_rect(fill=NA, size=0.5, linetype=1, color="white"),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18))
ggsave(paste0(path2save,"forecast_accuracy_boxplot.png"), height = 5, width=6.5, units="in", dpi=72)



####
####  MAKE FIGURE S1
####
ggplot(plot_dat, aes(x=species, y=mae, fill=modelweather))+
  geom_bar(stat="identity", position="dodge", width=0.8, color="black")+
  scale_fill_manual(values=c("coral", "darkred", "dodgerblue", "darkblue"),
                    labels = c("IPM - no climate", "IPM - climate", "QBM - no climate", "QBM - climate"))+
  ylab(paste("Forecast Error (MAE)"))+
  xlab("Species")+
  theme_few()+
  theme(legend.position=c(0.87,0.88),
        legend.background = element_rect(fill=NA, size=0.5, linetype=1, color="grey"),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))
ggsave(paste0(path2save,"forecast_errors.png"), height = 5, width=6.5, units="in", dpi=72)


####
####  SIGNIFICANCE TESTS FOR ACCURACY MEASURES
####
ipm_quadyear_weather <- ddply(ipm_onestep, .(quad, t1, species), summarise,
                              median_pred = median(cover.t1),
                              obs = mean(obs.cover.t1))
ipm_quadyear_noweather <- ddply(ipm_onestep_dd, .(quad, t1, species), summarise,
                                median_pred = median(cover.t1),
                                obs = mean(obs.cover.t1))
test1 <- subset(ipm_quadyear_noweather, species=="POSE")
test2 <- subset(ipm_quadyear_weather, species=="POSE")
rho_comp(x1=test2$median_pred, x2=test1$median_pred, y=test1$obs)

qbm_quadyear_weather <- ddply(qbm_onestep, .(quad, year, species), summarise,
                              median_pred = median(pred_cover.t1),
                              obs = mean(obs_cover.t1))
qbm_quadyear_noweather <- ddply(qbm_onestep_dd, .(quad, year, species), summarise,
                                median_pred = median(pred_cover.t1),
                                obs = mean(obs_cover.t1))
test1 <- subset(ipm_quadyear_noweather, species=="POSE")
test2 <- subset(ipm_quadyear_weather, species=="POSE")
rho_comp(x1=test2$median_pred, x2=test1$median_pred, y=test1$obs)

ipm4merge_noweather <- ipm_quadyear_noweather
ipm4merge_noweather$t1 <- ipm4merge_noweather$t1+1900-1
qbmVipm_df_noweather <- merge(ipm4merge_noweather, qbm_quadyear_noweather, 
                              by.x=c("quad", "t1", "species"),
                              by.y=c("quad", "year", "species"))
test1 <- subset(qbmVipm_df_noweather, species=="POSE")
rho_comp(x1=test1$median_pred.x, x2=test1$median_pred.y, y=test1$obs.x)

ipm4merge_weather <- ipm_quadyear_weather
ipm4merge_weather$t1 <- ipm4merge_weather$t1+1900-1
qbmVipm_df_weather <- merge(ipm4merge_weather, qbm_quadyear_weather, 
                            by.x=c("quad", "t1", "species"),
                            by.y=c("quad", "year", "species"))
test1 <- subset(qbmVipm_df_weather, species=="POSE")
rho_comp(x1=test1$median_pred.x, x2=test1$median_pred.y, y=test1$obs.x)



####
####  SIGNIFICANCE TESTS FOR ERROR COMPARISONS (SI MATERIAL)
####
##  Function from Ye et al. ---
print_comparison_table <- function()
{
  compute_p_values <- function(x1, x2, y)
  {
    index <- is.finite(x1) & is.finite(x2) & is.finite(y)
    x1 <- x1[index]
    x2 <- x2[index]
    y <- y[index]
    err1 <- abs(y - x1)
    err2 <- abs(y - x2)
    mae_ttest <- t.test(err1, err2, paired = TRUE, alternative = "less")
    mae_df <- mae_ttest$parameter
    mae_statistic <- mae_ttest$statistic
    mae_p <- mae_ttest$p.value
    rho_ttest <- rho_comp(x1, x2, y)
    rho_df <- rho_ttest$df
    rho_statistic <- rho_ttest$statistic
    rho_p <- rho_ttest$p.value
    return(data.frame(mae_df, mae_statistic, mae_p, rho_df, rho_statistic, rho_p))
  }
  
  preds <- readRDS("all_model_predictions.RDS")
  
  ## ONLY NORMALIZE IF NOT COMPARING BY SPECIES...
  #   preds <- list(BOGR=subset(pred_df, species=="BOGR"),
  #                 HECO=subset(pred_df, species=="HECO"),
  #                 PASM=subset(pred_df, species=="PASM"),
  #                 POSE=subset(pred_df, species=="POSE"))
  #   # normalize by mean obs value
  #   preds_n <- lapply(names(preds), function(spp) {
  #     df <- preds[[spp]]
  #     sigma <- sd(df$obs, na.rm = TRUE)
  #     mu <- mean(df$obs, na.rm = TRUE)
  #     df$obs <- (df$obs - sigma) / mu
  #     df$ipm_dd <- (df$ipm_dd - sigma) / mu
  #     df$ipm_clim <- (df$ipm_clim - sigma) / mu
  #     df$qbm_dd <- (df$qbm_dd - sigma) / mu
  #     df$qbm_clim <- (df$qbm_clim - sigma) / mu
  #     df$species <- spp
  #     return(df)
  #   })
  #   preds_n <- do.call(rbind, preds_n)
  #   preds_n$species <- factor(preds_n$species)
  
  temp_table_all <- list()
  for(do_spp in unique(preds$species)){
    preds_n <- subset(preds, species==do_spp)
    compare_from <- list(preds_n$ipm_clim, 
                         preds_n$qbm_clim, 
                         preds_n$ipm_dd, 
                         preds_n$ipm_clim)
    compare_to <- list(preds_n$ipm_dd, 
                       preds_n$qbm_dd, 
                       preds_n$qbm_dd, 
                       preds_n$qbm_clim)
    comparison_names <- list("climate IPM vs. simple IPM", 
                             "climate QBM vs. simple QBM", 
                             "simple IPM vs. simple QBM", 
                             "climate IPM vs. climate QBM")
    temp_table <- do.call(rbind, lapply(1:4, function(i) {
      temp <- compute_p_values(compare_from[[i]], compare_to[[i]], preds_n$obs)
      return(data.frame(comparison = comparison_names[[i]], 
                        performance_measure = c("rho", "MAE"), 
                        test_type = "t",
                        test_statistic = c(temp$rho_statistic, temp$mae_statistic), 
                        df = c(temp$rho_df, temp$mae_df), 
                        p_value = c(temp$rho_p, temp$mae_p)))
    }))
    temp_table$species = do_spp
    temp_table_all <- rbind(temp_table_all, temp_table)
  }
  
  my_table <- xtable(temp_table_all, digits = 3)
  print(my_table, type = "html", file = "stats_table.html", include.rownames = FALSE)
  print(my_table, type = "latex", file = "stats_table.latex", include.rownames = FALSE)
  return()
}

print_comparison_table()
