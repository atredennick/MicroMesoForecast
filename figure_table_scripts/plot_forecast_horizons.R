##  R Script to plot forecast horizons for both model types and all species.
##
##  Author: Andrew Tredennick
##  Last update: 4-20-2016

# Clear the workspace
rm(list=ls())


####
####  PRELIMINARIES
####
path2ipms <- "../analysis/ipm/simulations/results/forecast_horizon/"
path2qbms <- "../analysis/quadBM/simulations/results/forecast_horizon/"
path2figs <- "../manuscript/components/"
species_list <- c("BOGR", "HECO", "PASM", "POSE")



####
####  LOAD LIBRARIES
####
library(ggplot2)  # plotting
library(ggthemes) # plotting themes
library(reshape2) # data wrangling
library(plyr)     # data wrangling/summarizing
library(xtable)   # writing latex/html tables



####
####  READ IN SIMULATION RESULTS AND COLLATE/FORMAT
####
# Saved simulation results are the final year (t) of cover predicted by the 
# model from some starting year t-x, where x ranges from 1 to 12 years. Each
# forecast simulation was run 50 times with different parameters so that we can
# plot uncertainty around the forecast horizon trend.

##  IPM results
out_ipm <- list() # empty storage for data frames
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2ipms,tmp_file))
  tmp_sims$error <- with(tmp_sims, pred.cover.t1-obs.cover.t1)
  tmp_agg <- ddply(tmp_sims, .(horizon,species), summarise,
                   num_preds = length(error),
                   mean_error = mean(abs(error)),
                   sd_error = sd(abs(error)),
                   rmse = sqrt(mean(error^2)),
                   rho = cor(pred.cover.t1,obs.cover.t1))
  out_ipm <- rbind(out_ipm, tmp_agg)
}


##  QBM results
out_qbm <- list()
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2qbms,tmp_file))
  # First aggregate over the reps since QBM is run over many reps, even with mean parameter values
  # NOTE: this is not necessary for the IPM because it is deterministic for a given set of parameter values
  tmp_agg1 <- ddply(tmp_sims, .(horizon,species,quad,year.start,year.pred), summarise,
                    mean_pred = mean(pred.cover.t1),
                    mean_obs = mean(obs.cover.t1))
  tmp_agg1$error <- with(tmp_agg1, mean_pred - mean_obs)
  tmp_agg <- ddply(tmp_agg1, .(horizon,species), summarise,
                   num_preds = length(error),
                   mean_error = mean(abs(error)),
                   sd_error = sd(abs(error)),
                   rmse = sqrt(mean(error^2)),
                   rho = cor(mean_pred,mean_obs))
  out_qbm <- rbind(out_qbm, tmp_agg)
}

##  Combine the results
max_ipm_year <- max(out_ipm$horizon)
max_qbm_year <- max(out_qbm$horizon)
out_ipm$model <- "AIPM"
out_qbm$model <- "BQBM"
out_all <- rbind(out_ipm, out_qbm)



####
####  MAKE THE PLOT
####
out_plot <- subset(out_all, num_preds>49 & horizon<9)
mycols <- c("#969C43", "#B46FA1")
mycols <- c("grey50", "black")
ggplot(data=out_plot, aes(x=horizon, y=rho))+
  # geom_ribbon(aes(x=yearsbefore, ymax=(mean_error+sd_error)*100, 
                  # ymin=(mean_error-sd_error)*100), alpha=0.3, color=NA)+
  geom_hline(aes(yintercept=0.5), linetype="dashed", color="dodgerblue")+
  geom_line(aes(linetype=model))+
  geom_point(aes(shape=model),size=2)+
  facet_wrap("species", ncol=1)+
  xlab("Forecasted Years")+
  ylab(expression(paste("Forecast Accuracy (", rho, ")")))+
  # ylab("Mean Absolute Error (% Cover)")+
  scale_x_continuous(breaks=seq(1,13,by=1))+
  scale_y_continuous(limits=c(0,1))+
  scale_shape_discrete(name="", labels=c("IPM", "QBM"))+
  scale_linetype_discrete(name="", labels=c("IPM", "QBM"))+
  theme_few()+
  theme(legend.position=c(0.8,0.96),
        legend.background=element_rect(NA))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))
ggsave(paste0(path2figs,"forecast_horizon.png"), width = 7, height = 16, units = "cm", dpi=300)

ggsave(paste0(path2figs,"forecast_horizon.pdf"), width = 3, height = 8, units = "in", dpi=600)



####
####  FOR PRESENTATIONS
####
mycols <- c("#969C43", "#B46FA1")
ggplot(data=out_plot, aes(x=horizon, y=rho))+
  geom_hline(aes(yintercept=0.5), linetype="dashed", color="grey50")+
  geom_line(aes(color=model))+
  # geom_point(size=2)+
  facet_wrap("species", ncol=2)+
  xlab("Forecasted Years")+
  ylab(expression(paste("Forecast Accuracy (", rho, ")")))+
  scale_x_continuous(breaks=seq(1,13,by=1))+
  scale_y_continuous(limits=c(0,1))+
  # scale_shape_discrete(name="", labels=c("IPM", "QBM"))+
  scale_color_manual(values = c("dodgerblue", "#E78C33"), labels=c("IPM", "QBM"), name = NULL)+
  theme_few()+
  theme(legend.background=element_rect(NA))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))





####
####  t-TESTS FOR EACH TIME HORIZON (TABLE S24)
####
##  Source statistical functions for comparing correlations
source("rho_comp_functions.R")

##  Read in IPM results
out_ipm <- list() # empty storage for data frames
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2ipms,tmp_file))
  out_ipm <- rbind(out_ipm, tmp_sims)
}
# get rid of sim.num column since not necessary
out_ipm <- subset(out_ipm, select = -c(sim.num))
# Rename IPM columns
colnames(out_ipm)[which(colnames(out_ipm)=="pred.cover.t1")] <- "pred.cover.t1.ipm"
colnames(out_ipm)[which(colnames(out_ipm)=="obs.cover.t1")] <- "obs.cover.t1.ipm"

##  QBM results
out_qbm <- list()
for(ispp in 1:length(species_list)){
  tmp_file <- paste0(species_list[ispp],"_final_year_cover.RDS")
  tmp_sims <- readRDS(paste0(path2qbms,tmp_file))
  # First aggregate over the reps since QBM is run over many reps, even with mean parameter values
  # NOTE: this is not necessary for the IPM because it is deterministic for a given set of parameter values
  tmp_agg1 <- ddply(tmp_sims, .(horizon,species,quad,year.start,year.pred), summarise,
                    mean_pred = mean(pred.cover.t1),
                    mean_obs = mean(obs.cover.t1))
  out_qbm <- rbind(out_qbm, tmp_agg1)
}
# Rename QBM columns
colnames(out_qbm)[which(colnames(out_qbm)=="mean_pred")] <- "pred.cover.t1.qbm"
colnames(out_qbm)[which(colnames(out_qbm)=="mean_obs")] <- "obs.cover.t1.qbm"

out_ipm$year.start <- out_ipm$year.start+1900-1 # make a calendar year for merging; minus 1 to match QBM
out_ipm$year.pred <- out_ipm$year.pred+1900-1 # make a calendar year for merging; minus 1 to match QBM
all_predictions <- merge(out_ipm, out_qbm)
saveRDS(all_predictions, "all_model_predictions_horizons.RDS")



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
#     mae_ttest <- t.test(err1, err2, paired = TRUE, alternative = "less")
#     mae_df <- mae_ttest$parameter
#     mae_statistic <- mae_ttest$statistic
#     mae_p <- mae_ttest$p.value
    rho_ttest <- as.data.frame(rho_comp(x1, x2, y))
    rho_df <- rho_ttest$df
    rho_statistic <- rho_ttest$statistic
    rho_p <- rho_ttest$p.value
    return(data.frame(rho_df, rho_statistic, rho_p))
  }
  
  preds <- readRDS("all_model_predictions_horizons.RDS")
  # Exclude long horizons where we don't have enough predictions
  preds <- subset(preds, horizon < 9)
  
  temp_table_all <- list()
  for(do_horizon in unique(preds$horizon)){
    preds_h <- subset(preds, horizon==do_horizon)
    for(do_spp in unique(preds$species)){
      preds_n <- subset(preds_h, species==do_spp)
      compare_from <- list(preds_n$pred.cover.t1.ipm)
      compare_to <- list(preds_n$pred.cover.t1.qbm)
      comparison_names <- list("climate IPM vs. climate QBM")
      temp_table <- do.call(rbind, lapply(1:1, function(i) {
        temp <- compute_p_values(compare_from[[i]], compare_to[[i]], preds_n$obs.cover.t1.qbm)
        return(data.frame(comparison = comparison_names[[i]], 
                          performance_measure = c("rho"), 
                          test_type = "t",
                          test_statistic = c(temp$rho_statistic), 
                          df = c(temp$rho_df), 
                          p_value = c(temp$rho_p)))
      }))
      temp_table$species = do_spp
      temp_table$horizon = do_horizon
      temp_table_all <- rbind(temp_table_all, temp_table)
    } # end spp loop
  } # end horizon loop
  
  
  my_table <- xtable(temp_table_all, digits = 3)
  print(my_table, type = "html", file = "horizon_stats_table.html", include.rownames = FALSE)
  print(my_table, type = "latex", file = "horizon_stats_table.latex", include.rownames = FALSE)
  return()
}

print_comparison_table()


##  Function from Ye et al. ---
print_comparison_table_opposite <- function()
{
  compute_p_values <- function(x1, x2, y)
  {
    index <- is.finite(x1) & is.finite(x2) & is.finite(y)
    x1 <- x1[index]
    x2 <- x2[index]
    y <- y[index]
    err1 <- abs(y - x1)
    err2 <- abs(y - x2)
    #     mae_ttest <- t.test(err1, err2, paired = TRUE, alternative = "less")
    #     mae_df <- mae_ttest$parameter
    #     mae_statistic <- mae_ttest$statistic
    #     mae_p <- mae_ttest$p.value
    rho_ttest <- as.data.frame(rho_comp(x1, x2, y))
    rho_df <- rho_ttest$df
    rho_statistic <- rho_ttest$statistic
    rho_p <- rho_ttest$p.value
    return(data.frame(rho_df, rho_statistic, rho_p))
  }
  
  preds <- readRDS("all_model_predictions_horizons.RDS")
  # Exclude long horizons where we don't have enough predictions
  preds <- subset(preds, horizon < 9)
  
  temp_table_all <- list()
  for(do_horizon in unique(preds$horizon)){
    preds_h <- subset(preds, horizon==do_horizon)
    for(do_spp in unique(preds$species)){
      preds_n <- subset(preds_h, species==do_spp)
      compare_from <- list(preds_n$pred.cover.t1.qbm)
      compare_to <- list(preds_n$pred.cover.t1.ipm)
      comparison_names <- list("climate QBM vs. climate IPM")
      temp_table <- do.call(rbind, lapply(1:1, function(i) {
        temp <- compute_p_values(compare_from[[i]], compare_to[[i]], preds_n$obs.cover.t1.qbm)
        return(data.frame(comparison = comparison_names[[i]], 
                          performance_measure = c("rho"), 
                          test_type = "t",
                          test_statistic = c(temp$rho_statistic), 
                          df = c(temp$rho_df), 
                          p_value = c(temp$rho_p)))
      }))
      temp_table$species = do_spp
      temp_table$horizon = do_horizon
      temp_table_all <- rbind(temp_table_all, temp_table)
    } # end spp loop
  } # end horizon loop
  
  
  my_table <- xtable(temp_table_all, digits = 3)
  print(my_table, type = "html", file = "horizon_stats_table_opposite.html", include.rownames = FALSE)
  print(my_table, type = "latex", file = "horizon_stats_table_opposite.latex", include.rownames = FALSE)
  return()
}

print_comparison_table_opposite()

