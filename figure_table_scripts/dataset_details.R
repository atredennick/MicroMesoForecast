##  Gets data information by species and vital rate

rm(list=ls())

library(plyr)

# Growth
grow_dat <- readRDS("../analysis/processed_data/grow_with_weather.RDS")
grow_obs <- ddply(grow_dat, .(species), summarise,
                  nobs = length(area.t1),
                  nquads = length(unique(quad)))

# Survival
surv_dat <- readRDS("../analysis/processed_data/surv_with_weather.RDS")
surv_obs <- ddply(surv_dat, .(species), summarise,
                  nobs = length(area),
                  nquads = length(unique(quad)))

# Recruitment
nquad <- numeric(4)
counter=1
for(spp in c("BOGR", "HECO", "PASM", "POSE")){
  if(spp=="BOGR") rec_tmp <- read.csv("../analysis/data_processing/speciesData/BOGR/edited/recArea.csv") 
  if(spp!="BOGR") rec_tmp <- read.csv(paste0("../analysis/data_processing/speciesData/",spp,"/recArea.csv"))
  nquad[counter] <- length(unique(rec_tmp$quad))
  counter=counter+1
}
rec_dat <- readRDS("../analysis/processed_data/rec_with_weather.RDS")
rec_obs <- ddply(rec_dat, .(species), summarise,
                  nobs = length(parents1))

# Percent cover
cover_dat <- readRDS("../analysis/processed_data/cover_with_weather.RDS")
cover_obs <- ddply(cover_dat, .(species), summarise,
                   nobs = length(propCover.t1),
                   nquads = length(unique(quad)))
