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
                 observation = mean(obs.cover.t1),
                 quant_dist = quantile(cover.t1*100, probs = 0.95) - quantile(cover.t1, probs = 0.05))
# Calculate error
avg_ipm$error <- avg_ipm$observation*100 - avg_ipm$avg_prediction*100
# Aggregate over species
ipmstats <- ddply(avg_ipm, .(species), summarise,
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
# Aggregate over species
qbmstats <- ddply(avg_qbm, .(species), summarise,
                  mae = mean(abs(error)),
                  rho = cor(observation, avg_prediction, method = "p"),
                  mean_cover = mean(observation*100),
                  quant_dist = mean(quant_dist))
qbmstats$model <- "QBM"

# REPLACE mean QBM cover with cover from IPM because it is the correct value
qbmstats$mean_cover <- ipmstats$mean_cover


testq <- subset(qbm_onestep, quad=="B1"&sim==1)
testi <- subset(ipm_onestep, quad=="B1"&rep==1)


####
####  CALCULATE P-VALUES FROM T-TESTS ON MODEL ERRORS
####  [SEE YE et al. 2015 (PNAS) FOR EXAMPLE USING MODEL ACCURACIES]
####
# mae_list <- list()
# mae_pvals <- list()
# species <- unique(qbmstats$species)
# for(spp in species){
#   tmp1 <- subset(avg_ipm, species==spp)
#   tmp2 <- subset(avg_qbm, species==spp)
#   err1 <- as.data.frame(abs(tmp1$observation - tmp1$avg_prediction))
#   err2 <- as.data.frame(abs(tmp2$observation - tmp2$avg_prediction))
#   colnames(err1) <- "resid"
#   err1$model <- "ipm"
#   colnames(err2) <- "resid"
#   err2$model <- "qbm"
#   errd <-  rbind(err1, err2)
#   mae_ttest <- with(errd, t.test(resid~model, alternative = "less"))
#   mae_list[[spp]] <- mae_ttest
#   mae_pvals[[spp]] <- mae_ttest[["p.value"]]
# }
# 
# mae_ps <- melt(mae_pvals)



####
####  CALCULATE STATISTICS FROM NO-WEATHER SIMULATIONS (DENS. DEP. ONLY)
####
# IPM
path2ipm2 <- "../analysis/ipm/simulations/results/one_step_validation_densdep_only/"
# Read in data
ipm_onestep_dd <- readRDS(paste0(path2ipm2,"ipm_loyo_forecasts_noweather_combined.RDS"))
# Average predictions over quad-year reps
avg_ipm_dd <- ddply(ipm_onestep_dd, .(species, quad, t1), summarise,
                 avg_prediction = median(cover.t1),
                 observation = mean(obs.cover.t1),
                 quant_dist = quantile(cover.t1*100, probs = 0.95) - quantile(cover.t1, probs = 0.05))
# Calculate error
avg_ipm_dd$error <- avg_ipm_dd$observation*100 - avg_ipm_dd$avg_prediction*100
# Aggregate over species
ipmstats_dd <- ddply(avg_ipm_dd, .(species), summarise,
                  mae = mean(abs(error)),
                  rho = cor(observation, avg_prediction, method = "p"),
                  mean_cover = mean(observation*100),
                  quant_dist = mean(quant_dist))
ipmstats_dd$model <- "IPM"

#QBM
path2qbm2 <- "../analysis/quadBM/simulations/results/one_step_validation_densdep_only/"
qbm_onestep_dd <- readRDS(paste0(path2qbm2, "qbm_one-step_forecasts_combined.RDS"))
# Average predictions over quad-year reps
avg_qbm_dd <- ddply(qbm_onestep_dd, .(species, quad, year), summarise,
                 avg_prediction = median(pred_cover.t1),
                 observation = mean(obs_cover.t1),
                 quant_dist = quantile(pred_cover.t1*100, probs = 0.95) - quantile(pred_cover.t1, probs = 0.05))
# Calculate error
avg_qbm_dd$error <- avg_qbm_dd$observation*100 - avg_qbm_dd$avg_prediction*100
# Aggregate over species
qbmstats_dd <- ddply(avg_qbm_dd, .(species), summarise,
                  mae = mean(abs(error)),
                  rho = cor(observation, avg_prediction, method = "p"),
                  mean_cover = mean(observation*100),
                  quant_dist = mean(quant_dist))
qbmstats_dd$model <- "QBM"

# REPLACE mean QBM cover with cover from IPM because it is the correct value
qbmstats_dd$mean_cover <- ipmstats_dd$mean_cover



####
####  MAKE A FIGURE
####
ipmstats$weather <- "yes"
ipmstats_dd$weather <- "no"
qbmstats$weather <- "yes"
qbmstats_dd$weather <- "no"
plot_dat <- rbind(ipmstats,ipmstats_dd, qbmstats,qbmstats_dd)
plot_dat$modelweather <- paste0(plot_dat$model,plot_dat$weather)

library(ggthemes)
library(gridExtra)
ggplot(plot_dat, aes(x=species, y=rho, fill=modelweather))+
  geom_bar(stat="identity", position="dodge", width=0.8, color="black")+
  scale_fill_manual(values=c("coral", "darkred", "dodgerblue", "darkblue"),
                    labels = c("IPM - no climate", "IPM - climate", "QBM - no climate", "QBM - climate"))+
  ylab(expression(paste("Forecast Accuracy (", rho, ")")))+
  xlab("Species")+
  theme_few()+
  theme(legend.position=c(0.87,0.88),
        legend.background = element_rect(fill=NA, size=0.5, linetype=1, color="grey"),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))
ggsave(paste0(path2save,"forecast_accuracy.png"), height = 5, width=6.5, units="in", dpi=72)



####
####  SIGNIFICANCE TESTS FOR ACCURACY MEASURES
####
# from Wilcox' Robust Statistics R package
elimna<-function(m)
{
  #
  # remove any rows of data having missing values
  #
  if(is.null(dim(m)))m<-as.matrix(m)
  ikeep<-c(1:nrow(m))
  for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
  elimna<-m[ikeep[ikeep>=1],]
  elimna
}

pcorhc4_err <- function(x, y, CN = FALSE)
{
  z1 <- (x - mean(x)) / sd(x)
  z2 <- (y - mean(y)) / sd(y)
  ans <- olshc4(z1, z2, alpha = 0.05, CN = CN)
  return(ans$ci[2, 6])
}

# from Wilcox' Robust Statistics R package
pcorhc4<-function(x,y,alpha=.05,CN=F)
{
  #
  #   Compute a .95 confidence interval for Pearson's correlation coefficient.
  #   using the HC4 method
  #
  # CN=F, degrees of freedom are n-p; seems better for general use.
  # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
  #
  xy<-elimna(cbind(x,y))
  x<-xy[,1]
  y<-xy[,2]
  z1=(x-mean(x))/sqrt(var(x))
  z2=(y-mean(y))/sqrt(var(y))
  ans=olshc4(z1,z2,alpha=alpha,CN=CN)
  list(r=ans$r,ci=ans$ci[2,3:4],p.value=ans$ci[2,5])
}
TWOpov_err <- function(x,y,CN=F)
{
  #
  # Comparing two dependent correlations: Overlapping case
  #
  # x is assumed to be a matrix with 2 columns
  #
  #  Compare correlation of x[,1] with y to x[,2] with y
  #
  # returns p-value
  if(!is.matrix(x))stop("x should be a matrix")
  if(ncol(x)!=2)stop("x should be a matrix with two columns")
  xy=elimna(cbind(x,y))
  x1=xy[,1]
  x2=xy[,2]
  y=xy[,3]
  r12=cor(x1,y)
  r13=cor(x2,y)
  r23=cor(x1,x2)
  err12 <- pcorhc4_err(x1,y,CN=CN)
  err13 <- pcorhc4_err(x2,y,CN=CN)
  corhat=((r23-.5*r12*r13)*(1-r12^2-r13^2-r23^2)+r23^3)/((1-r12^2)*(1-r13^2))
  err_correction_term = 2*corhat*(err12)*(err13)
  err_diff <- sqrt(err12^2 + err13^2 - err_correction_term)
  return(err_diff)
}
# from Wilcox' Robust Statistics R package
olshc4<-function(x,y,alpha=.05,CN=FALSE,xout=FALSE,outfun=outpro,HC3=FALSE,...)
{
  #
  # Compute confidence for least squares
  # regression using heteroscedastic method
  # recommended by Cribari-Neto (2004).
  # CN=F, degrees of freedom are n-p
  # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
  # All indications are that CN=F is best for general use.
  #
  #  HC3=TRUE, will replace the HC4 estimator with the HC3 estimator.
  #
  x<-as.matrix(x)
  if(nrow(x) != length(y))stop("Length of y does not match number of x values")
  m<-cbind(x,y)
  m<-elimna(m)
  y<-m[,ncol(x)+1]
  x=m[,1:ncol(x)]
  n=length(y)
  nrem=n
  n.keep=length(y)
  x<-as.matrix(x)
  if(xout){
    flag<-outfun(x,...)$keep
    x<-as.matrix(x)
    x<-x[flag,]
    y<-y[flag]
    n.keep=length(y)
    x<-as.matrix(x)
  }
  temp<-lsfit(x,y)
  x<-cbind(rep(1,nrow(x)),x)
  xtx<-solve(t(x)%*%x)
  h<-diag(x%*%xtx%*%t(x))
  n<-length(h)
  d<-(n*h)/sum(h)
  for(i in 1:length(d)){
    d[i]<-min(4, d[i])
  }
  if(HC3)d=2
  hc4<-xtx%*%t(x)%*%diag(temp$res^2/(1-h)^d)%*%x%*%xtx
  df<-nrow(x)-ncol(x)
  crit<-qt(1-alpha/2,df)
  if(CN)crit=qnorm(1-alpha/2)
  al<-ncol(x)
  p=al-1
  ci<-matrix(NA,nrow=al,ncol=6)
  lab.out=rep("Slope",p)
  dimnames(ci)<-list(c("(Intercept)",lab.out),c("Coef.","Estimates",
                                                "ci.lower","ci.upper","p-value","Std.Error"))
  for(j in 1:al){
    ci[j,1]<-j-1
    ci[j,2]<-temp$coef[j]
    ci[j,3]<-temp$coef[j]-crit*sqrt(hc4[j,j])
    ci[j,4]<-temp$coef[j]+crit*sqrt(hc4[j,j])
    test<-temp$coef[j]/sqrt(hc4[j,j])
    ci[j,5]<-2*(1-pt(abs(test),df))
    if(CN)ci[j,5]<-2*(1-pnorm(abs(test),df))
  }
  ci[,6]=sqrt(diag(hc4))
  list(n=nrem,n.keep=n.keep,ci=ci, cov=hc4)
}

rho_comp <- function(x1, x2, y)
  # computes p-value for cor(x1, y) > cor(x2, y) using 
  # t-test with df = length(y) - 2
{
  if(identical(x1, x2))
    return(0.5)
  n <- sum(is.finite(x1) & is.finite(x2) & is.finite(y))
  x1y <- cor(x1, y, use = "pairwise")
  x2y <- cor(x2, y, use = "pairwise")
  err <- TWOpov_err(x=as.matrix(cbind(x1, x2)), y)
  p_value <- 1 - pt((x1y - x2y) / err, df = n-2, lower.tail = TRUE)
  return(data.frame(df = n-2, statistic = (x1y - x2y) / err, p.value = p_value))
}

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
test1 <- subset(qbm_quadyear_noweather, species=="BOGR")
test2 <- subset(qbm_quadyear_weather, species=="BOGR")
rho_comp(x1=test2$median_pred, x2=test1$median_pred, y=test1$obs)


mae_list <- list()
mae_pvals <- list()
species <- unique(qbmstats$species)
for(spp in species){
  tmp1 <- subset(avg_ipm, species==spp)
  tmp2 <- subset(avg_ipm_dd, species==spp)
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
