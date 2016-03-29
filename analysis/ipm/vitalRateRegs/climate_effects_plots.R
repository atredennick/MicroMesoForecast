##  Script to plot climate effects from each stastical model

library(ggmcmc)
library(gridExtra)
library(plyr)
library(reshape2)

setwd("~/Repos/MicroMesoForecast/manuscript")

allfiles <- list.files("../analysis/speciesData/")
removals <- grep("csv", allfiles)
species_list <- allfiles[-removals]

clim_covs <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2",
               "ppt1xTmeanSpr1", "ppt2xTmeanSpr2",
               "sizeXpptLag", "sizeXppt1", "sizeXppt2",
               "sizeXTmeanSpr1", "sizeXTmeanSpr2")

## Plot function
make_plot <- function(data, col="black"){
  gtmp <- ggplot(data, aes(x = med, y = Parameter)) + 
    geom_vline(aes(xintercept=0), linetype=2)+
    geom_point(size = 3, color=col) + 
    geom_segment(aes(x = Low, xend = High, yend = Parameter), size = 1.5, color=col) + 
    geom_segment(aes(x = low, xend = high, yend = Parameter), size = 0.5, color=col) + 
    xlab("HPD") +
    ylab("Parameter")+
    ggtitle(spp)
}

#   gtmp <- ggplot(tmp_agg, aes(x = med, y = reorder(Parameter, med))) + 
#     geom_point(size = 3) + 
#     geom_segment(aes(x = Low, xend = High, yend = reorder(Parameter, med)), size = 1.5) + 
#     geom_segment(aes(x = low, xend = high, yend = reorder(Parameter, med)), size = 0.5) + 
#     geom_vline(aes(xintercept=0))+
#     xlab("HPD") +
#     ylab("Parameter")+
#     ggtitle(spp)

##  Growth
glist <- list()
for(spp in species_list){
  tmp_data <- readRDS(paste("../analysis/ipm/vitalRateRegs/growth/growth_stanmcmc_", spp, ".RDS", sep=""))
  tmp_clim <- tmp_data[grep("b2", tmp_data[,"Parameter"]), ]
  tmp_clim[,"Parameter"] <- rep(clim_covs, each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   med = median(value),
                   high = quantile(value, 0.975),
                   High = quantile(value, 0.875),
                   low = quantile(value, 0.025),
                   Low = quantile(value, 0.125))
  glist[[spp]] <- make_plot(tmp_agg, col="black")
}
# g_grow <- grid.arrange(glist[[1]], glist[[2]], glist[[3]], glist[[4]],
#                        ncol=4, nrow=1)


##  Survival
slist <- list()
params <- character(length = 12)
for(i in 1:12){
  params[i] <- paste("b2[",i,"]", sep="")
}
for(spp in species_list){
  tmp_data <- ggs(readRDS(paste("../analysis/ipm/vitalRateRegs/survival/survival_stanmcmc_", spp, ".RDS", sep="")))
  keeps <- which(tmp_data$Parameter %in% params)
  tmp_clim <- tmp_data[keeps, ]
  tmp_clim[,"Parameter"] <- rep(clim_covs, each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   med = median(value),
                   high = quantile(value, 0.975),
                   High = quantile(value, 0.875),
                   low = quantile(value, 0.025),
                   Low = quantile(value, 0.125))
  slist[[spp]] <- make_plot(tmp_agg, col="steelblue")
}
# g_surv <- grid.arrange(slist[[1]], slist[[2]], slist[[3]], slist[[4]],
#                        ncol=4, nrow=1)


##  Recruitment
rlist <- list()
params <- character(length = 12)
for(i in 1:7){
  params[i] <- paste("b2[",i,"]", sep="")
}
for(spp in species_list){
  tmp_data <- readRDS(paste("../analysis/ipm/vitalRateRegs/recruitment/recruitment_stanmcmc_", spp, ".RDS", sep=""))
  keeps <- which(tmp_data$Parameter %in% params)
  tmp_clim <- tmp_data[keeps, ]
  tmp_clim[,"Parameter"] <- rep(clim_covs[1:7], each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   med = median(value),
                   high = quantile(value, 0.975),
                   High = quantile(value, 0.875),
                   low = quantile(value, 0.025),
                   Low = quantile(value, 0.125))
  rlist[[spp]] <- make_plot(tmp_agg, col="coral")
}
# g_surv <- grid.arrange(rlist[[1]], rlist[[2]], rlist[[3]], rlist[[4]],
#                        ncol=4, nrow=1)


## Full plot
# all_plot <- grid.arrange(glist[[1]], glist[[2]], glist[[3]], glist[[4]],
#                          slist[[1]], slist[[2]], slist[[3]], slist[[4]],
#                          rlist[[1]], rlist[[2]], rlist[[3]], rlist[[4]],
#                          ncol=4, nrow=3)

grid.arrange(glist[[1]], slist[[1]], rlist[[1]], 
                         glist[[2]], slist[[2]], rlist[[2]],
                         glist[[3]], slist[[3]], rlist[[3]],
                         glist[[4]], slist[[4]], rlist[[4]],
                         ncol=3, nrow=4)

g <- arrangeGrob(glist[[1]], slist[[1]], rlist[[1]], 
                 glist[[2]], slist[[2]], rlist[[2]],
                 glist[[3]], slist[[3]], rlist[[3]],
                 glist[[4]], slist[[4]], rlist[[4]],
                 ncol=3, nrow=4) #generates g
# ggsave(file="components/figure/ipm_climeffs.pdf", g, width = 10, height = 10)
g
