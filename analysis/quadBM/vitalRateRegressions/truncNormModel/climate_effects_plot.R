##  Script to plot climate effects from the quad model

library(ggmcmc)
library(gridExtra)

setwd("~/Repos/MicroMesoForecast/manuscript")

allfiles <- list.files("../analysis/speciesData/")
removals <- grep("csv", allfiles)
species_list <- allfiles[-removals]

clim_covs <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2",
               "ppt1xTmeanSpr1", "ppt2xTmeanSpr2",
               "coverXpptLag", "coverXppt1", "coverXppt2",
               "coverXTmeanSpr1", "coverXTmeanSpr2")

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


glist <- list()
for(spp in species_list){
  tmp_data <- readRDS(paste("../analysis/quadBM/vitalRateRegressions/truncNormModel/popgrowth_stanmcmc_", spp, ".RDS", sep=""))
  tmp_clim <- tmp_data[grep("b2", tmp_data$Parameter), ]
  tmp_clim[,"Parameter"] <- rep(clim_covs, each=3000)
  tmp_agg <- ddply(tmp_clim, .(Parameter), summarise,
                   med = median(value),
                   high = quantile(value, 0.975),
                   High = quantile(value, 0.875),
                   low = quantile(value, 0.025),
                   Low = quantile(value, 0.125))
  glist[[spp]] <- make_plot(tmp_agg, col="black")
}
grid.arrange(glist[[1]], glist[[2]], glist[[3]], glist[[4]],
                       ncol=2, nrow=2)
g <- arrangeGrob(glist[[1]], glist[[2]], glist[[3]], glist[[4]],
                 ncol=2, nrow=2)
# ggsave(file="components/figure/qbm_climeffs.pdf", g, width = 7, height = 5)

