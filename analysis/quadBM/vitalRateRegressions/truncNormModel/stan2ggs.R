##  Script to unlist big stan fit objects and save as ggs files

library(rstan)
library(ggmcmc)

bigstan <- readRDS("popgrowth_stanfits.RDS")
spp_list <- c("BOGR", "HECO", "PASM", "POSE")

for(spp in spp_list){
  tmp <- bigstan[[spp]]
  long <- ggs(tmp)
  file <- paste("mcmc_popgrowth_", spp, ".RDS", sep="")
  saveRDS(long, file)
}