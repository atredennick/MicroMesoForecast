##  Quick script for checking model diagnostics and MCMC chains
##    for full quad-based model.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com

rm(list=ls(all=TRUE))

####
####  Load libraries ---------------------------------------------
####
library('ggplot2')
library('reshape2')


####
#### Read in diagnostics and MCMC chain -------------------------
####
mcmc_chains <- as.data.frame(readRDS("growthParamsMCMC.rds"))
chain_length <- nrow(mcmc_chains)/3
mcmc_chains$iter <- rep(c(1:chain_length),3)
mcmc_chains$chain <- as.character(rep(c(1:3),each=chain_length))
mcmc_melt <- melt(mcmc_chains, id.vars = c("iter", "chain"))
ggplot(mcmc_melt, aes(x=iter, y=value, color=chain, group=chain))+
  geom_line(alpha=0.5)+
  facet_wrap("variable", scales="free")+
  guides(color=FALSE)

## For saving (commented out for now)
# png("quad_model_mcmcs.png", width = 18, height = 15, units = "in", res = 100)
# ggplot(mcmc_melt, aes(x=iter, y=value, color=chain, group=chain))+
#   geom_line(alpha=0.5)+
#   facet_wrap("variable", scales="free")+
#   guides(color=FALSE)
# dev.off()

gelman_diags <- read.csv("growthGelman.csv")
ggplot(gelman_diags)+
  geom_histogram(aes(x=Point.est.), color="grey")
summary(gelman_diags$Point.est.)
