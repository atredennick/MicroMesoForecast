##  R Script to Plot Posterior Densities of Climate Effects
##
##  Author: Andrew Tredennick
##  Last update: 4-15-2016


# Clear the Workspace
rm(list=ls())



####
####  LOAD LIBRARIES
####
library(ggplot2)
library(ggthemes)
library(plyr)
library(reshape2)
library(ggmcmc)
library(gridExtra)



####
####  PRELIMINARIES
####
species <- c("BOGR", "HECO", "PASM", "POSE")
vitals <- c("growth", "survival", "recruitment", "cover")
path2ipms <- "../analysis/ipm/vitalRateRegs/"
path2qbms <- "../analysis/quadBM/vitalRateRegressions/truncNormModel/"
all_priors <- read.csv("../analysis/all_maxlppds.csv")


####
#### LOOP OVER MODELS AND SPECIES; PLOT THEM
####
pdf(file = "clim_posteriors.pdf", onefile = TRUE, height=3.5, width=3)
for(do_species in species){
  for(do_vital in vitals){
    if(do_vital=="cover") { path <- paste0(path2qbms, "popgrowth_stanmcmc_",do_species,".RDS") }
    if(do_vital!="cover") { path <- paste0(path2ipms,do_vital,"/",do_vital,"_stanmcmc_",do_species,".RDS") }
    if(do_vital=="survival"){ post_df <- ggs(readRDS(path)) }
    if(do_vital=="recruitment"){ post_df <- ggs(readRDS(path)) }
    if(do_vital=="cover"){ post_df <- readRDS(path) }
    if(do_vital=="growth"){ post_df <- readRDS(path) }
    
    post_clims <- post_df[grep("b2",post_df$Parameter),]
    post_clim_stats <- ddply(post_clims, .(Parameter), summarise,
                             average = mean(value),
                             upperci = quantile(value, 0.85),
                             lowerci = quantile(value, 0.15))
    post_clim_stats$col <- NA
    for(i in 1:nrow(post_clim_stats)){
      if(sign(post_clim_stats$lowerci[i]) == sign(post_clim_stats$upperci[i]) & post_clim_stats$average[i] > 0 ){
        post_clim_stats$col[i] <- "Cblue"
      }
      if(sign(post_clim_stats$lowerci[i]) == sign(post_clim_stats$upperci[i]) & post_clim_stats$average[i] < 0 ){
        post_clim_stats$col[i] <- "Ared"
      }
      if(sign(post_clim_stats$lowerci[i]) != sign(post_clim_stats$upperci[i]) ){
        post_clim_stats$col[i] <- "Bblack"
      }
    }
    
    post_clim_stats$param_name <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2",
                                    "TmeanSp1 \nX ppt1", "TmeanSpr2 \nX ppt2")
    post_clim_stats$rank <- c(1:nrow(post_clim_stats))
    post_clims <- merge(post_clims, post_clim_stats, by="Parameter")
    
    cols2use <- substr(sort(unique(post_clim_stats$col)),2,length(post_clim_stats$col))
    
      my_labeller <- function(variable, value){
        return(as.character(post_clim_stats$param_name))
      } 
    prior_sd <- as.numeric(subset(all_priors, species==do_species&vital==do_vital)["prior_stdev"])
    gout <- ggplot(post_clims, aes(x=value))+
      # geom_line(data=data.frame(x=rnorm(3000,0,prior_sd)), aes(x=x), stat = "density", adjust=3, color="darkorange", linetype=2)+
      geom_density(aes(fill=col),color=NA, adjust=3, alpha=0.5)+
      geom_vline(aes(xintercept=0), linetype=2, color="grey25")+
      facet_grid(rank~., labeller=my_labeller, scales="free")+
      xlab("Standardized Coefficient Value")+
      ylab("Density")+
      scale_fill_manual(values=cols2use)+
      guides(fill=FALSE)+
      ggtitle(paste(do_species,do_vital))+
      scale_y_continuous(breaks=NULL)+
      theme_few()+
      theme(strip.background=element_rect(fill="white"))+
      theme(strip.text.y = element_text(size = 0, angle = 0))+
      theme(axis.title=element_text(size=0), text=element_text(size=16))
    print(gout)
  }
}
dev.off()


####
####  ARRANGE GGPLOT OBJECT AND SAVE
####
# gplot <- grid.arrange(glist[[1]],glist[[2]],glist[[3]],glist[[4]],
#                       glist[[5]],glist[[6]],glist[[7]],glist[[8]],
#                       glist[[9]],glist[[10]],glist[[11]],glist[[12]],
#                       glist[[13]],glist[[14]],glist[[15]],glist[[16]],
#                       ncol=4,nrow=4)

