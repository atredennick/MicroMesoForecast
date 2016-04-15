##  Figure script to look at correlations between IPM and QBM
##    climate effect estimates.

library(ggmcmc)
library(gridExtra)
library(ggplot2)
library(plyr)

allfiles <- list.files("../analysis/speciesData/")
removals <- grep("csv", allfiles)
species_list <- allfiles[-removals]

clim_covs <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2",
               "ppt1xTmeanSpr1", "ppt2xTmeanSpr2")
clim_mains <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2",
                "ppt1xTmeanSpr1", "ppt2xTmeanSpr2")


####
##  Get mean growth climate effects
####
growth_parms <- data.frame(species=NA, model=NA, parameter=NA, value=NA)
for(spp in species_list){
  tmp_data <- readRDS(paste("../analysis/ipm/vitalRateRegs/growth/growth_stanmcmc_", spp, ".RDS", sep=""))
  params <- as.character(tmp_data$Parameter)
  tmp_clim <- tmp_data[grep("b2", params), ]
  tmp_clim[,"Parameter"] <- rep(clim_covs, each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   value = median(value))
  tmp_agg2 <- tmp_agg[which(tmp_agg$Parameter %in% clim_mains), ]
  tmp_agg2$species <- spp
  tmp_agg2$model <- "ipm"
  colnames(tmp_agg2) <- tolower(colnames(tmp_agg2))
  growth_parms <- rbind(growth_parms,
                        tmp_agg2[,c("species", "model", "parameter", "value")])
}

####
##  Get mean survival climate effects
####
surv_parms <- data.frame(species=NA, model=NA, parameter=NA, value=NA)
params <- character(length = 12)
for(i in 1:12){
  params[i] <- paste("b2.",i,".", sep="")
}
for(spp in species_list){
  tmp_data <- ggs(readRDS(paste("../analysis/ipm/vitalRateRegs/survival/survival_stanmcmc_", spp, ".RDS", sep="")))
  keeps <- which(as.character(tmp_data$Parameter) %in% params)
  tmp_clim <- tmp_data[keeps, ]
  tmp_clim[,"Parameter"] <- rep(clim_covs, each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   value = median(value))
  tmp_agg2 <- tmp_agg[which(tmp_agg$Parameter %in% clim_mains), ]
  tmp_agg2$species <- spp
  tmp_agg2$model <- "ipm"
  colnames(tmp_agg2) <- tolower(colnames(tmp_agg2))
  surv_parms <- rbind(surv_parms,
                        tmp_agg2[,c("species", "model", "parameter", "value")])
}


####
##  Get mean recruitment climate effects
####
rec_parms <- data.frame(species=NA, model=NA, parameter=NA, value=NA)
params <- character(length = 12)
for(i in 1:12){
  params[i] <- paste("b2.",i,".", sep="")
}
for(spp in species_list){
  tmp_data <- ggs(readRDS(paste("../analysis/ipm/vitalRateRegs/recruitment/recruitment_stanmcmc_", spp, ".RDS", sep="")))
  keeps <- which(as.character(tmp_data$Parameter) %in% params)
  tmp_clim <- tmp_data[keeps, ]
  tmp_clim[,"Parameter"] <- rep(clim_mains, each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   value = median(value))
  tmp_agg2 <- tmp_agg[which(tmp_agg$Parameter %in% clim_mains), ]
  tmp_agg2$species <- spp
  tmp_agg2$model <- "ipm"
  colnames(tmp_agg2) <- tolower(colnames(tmp_agg2))
  rec_parms <- rbind(rec_parms,
                      tmp_agg2[,c("species", "model", "parameter", "value")])
}


####
##  Get QBM growth climate estimates
####
clim_covs2 <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2",
               "ppt1xTmeanSpr1", "ppt2xTmeanSpr2")

qbm_parms <- data.frame(species=NA, model=NA, parameter=NA, qbmvalue=NA)
params <- character(length = 12)
for(i in 1:12){
  params[i] <- paste("b2.",i,".", sep="")
}
for(spp in species_list){
  tmp_data <- readRDS(paste("../analysis/quadBM/vitalRateRegressions/truncNormModel/popgrowth_stanmcmc_", spp, ".RDS", sep=""))
  keeps <- which(tmp_data$Parameter %in% params)
  tmp_clim <- tmp_data[keeps, ]
  tmp_clim[,"Parameter"] <- rep(clim_covs2, each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   qbmvalue = median(value))
  tmp_agg2 <- tmp_agg[which(tmp_agg$Parameter %in% clim_mains), ]
  tmp_agg2$species <- spp
  tmp_agg2$model <- "qbm"
  colnames(tmp_agg2) <- tolower(colnames(tmp_agg2))
  qbm_parms <- rbind(qbm_parms,
                     tmp_agg2[,c("species", "model", "parameter", "qbmvalue")])
}


####
##  Format data, combine, and plot
####
growth_parms <- growth_parms[2:nrow(growth_parms), ]
growth_parms$vitalrate <- "growth"
surv_parms <- surv_parms[2:nrow(surv_parms), ]
surv_parms$vitalrate <- "surv"
rec_parms <- rec_parms[2:nrow(rec_parms), ]
rec_parms$vitalrate <- "rec"
ipm_parms <- rbind(growth_parms, surv_parms, rec_parms)

qbm_parms <- qbm_parms[2:nrow(qbm_parms), ]
qbm_parms$vitalrate <- "growth"
vitals <- c("surv", "rec")
qbmalone <- qbm_parms
for(dov in vitals){
  tmp <- qbmalone
  tmp$vitalrate <- dov
  qbm_parms <- rbind(qbm_parms, tmp)
}

ipm_parms <- ipm_parms[,c("species","parameter","value","vitalrate")]
qbm_parms <- qbm_parms[,c("species","parameter","qbmvalue","vitalrate")]
all_ests <- merge(ipm_parms, qbm_parms, by = c("species", "parameter", "vitalrate"))

# Calculate correlation for each group
cors <- ddply(all_ests, c("species", "vitalrate"), summarise, 
              cor = round(cor(value, qbmvalue), 2))

library(ggthemes)
ggplot(all_ests, aes(x=qbmvalue, y=value))+
  geom_point(size=3)+
  geom_smooth(method="lm", color="black", se=FALSE)+
  facet_grid(species~vitalrate)+
  geom_text(data=cors, aes(label=paste("r = ", cor, sep="")), x=-.05, y=0.12, size=6)+
  xlab("QBM Estimate")+
  ylab("IPM Estimate")+
  theme_few()

