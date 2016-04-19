
library(reshape2)
library(plyr)

fitthin <- data.frame(Iteration=NA, Chain=NA, Parameter=NA,
                      value=NA, keep=NA, species=NA)
for(ispp in spp_list){
  fitlong <- readRDS(paste("../vitalRateRegs/growth/growth_stanmcmc_", ispp, ".RDS", sep=""))
  fitlong$keep <- "no"
  keepseq <- seq(from = 1, to = nrow(fitlong), by = 1)
  fitlong[keepseq,"keep"] <- "yes"
  tmp <- subset(fitlong, keep=="yes")
  tmp$species <- ispp
  fitthin <- rbind(fitthin,tmp)
}
fitthin <- fitthin[2:nrow(fitthin),]

##  Break up MCMC into regression components
# Climate effects
climeff_grow <- fitthin[grep("b2", fitthin$Parameter),]

# Yearly cover (size) effects
coveff_grow <- fitthin[grep("b1[.]", fitthin$Parameter),]
coveff_grow$yearid <- substr(coveff_grow$Parameter, 4, length(coveff_grow$Parameter))
coveff_grow$yearid <- unlist(strsplit(coveff_grow$yearid, split='[.]'))

# Mean cover effect
covermu_grow <- fitthin[grep("b1_mu", fitthin$Parameter),]

# Yearly intercepts
intercept_grow <- fitthin[grep("a", fitthin$Parameter),]
intercept_grow <- subset(intercept_grow, Parameter!="a_mu")
intercept_grow <- subset(intercept_grow, Parameter!="tau")
intercept_grow$yearid <- substr(intercept_grow$Parameter, 3, length(intercept_grow$Parameter))
intercept_grow$yearid <- unlist(strsplit(intercept_grow$yearid, split='[.]'))

# Mean intercept
interceptmu_grow <- fitthin[grep("a_mu", fitthin$Parameter),]

# Crowding effects
crowd_grow <- fitthin[grep("w", fitthin$Parameter),]

# Group effects
group_grow <- fitthin[grep("gint", fitthin$Parameter),]
group_grow$groupid <- substr(group_grow$Parameter, 6, length(group_grow$Parameter))
group_grow$groupid <- unlist(strsplit(group_grow$groupid, split='[.]'))

# Size-based variance parameters
tau_grow <- fitthin[grep("tau", fitthin$Parameter),]
tauSize_grow <- subset(tau_grow, Parameter=="tauSize")
tau_grow <- subset(tau_grow, Parameter=="tau")

## Get rid of big objects
rm(list = c("tmp","fitthin","fitlong"))

##  Define function to format survival coefficients
getGrowCoefs <- function(doYear, groupnum){
  # Get random chain and iteration for this timestep
  tmp4chain <- subset(climeff_grow, species=="BOGR")
  randchain <- sample(x = tmp4chain$Chain, size = 1)
  randiter <- sample(x = tmp4chain$Iteration, size = 1)
  
  # Get random effects if doYear!=NA
  if(is.na(doYear)==FALSE){
    tmp_intercept <- subset(intercept_grow, yearid==doYear &
                                       Iteration==randiter &
                                       Chain==randchain)
    tmp_size <- subset(coveff_grow, yearid==doYear &
                               Iteration==randiter &
                               Chain==randchain)
  }
  
  # Set mean intercept and slope if doYear==NA
  if(is.na(doYear)==TRUE){
    tmp_intercept <- subset(interceptmu_grow, Iteration==randiter &
                                         Chain==randchain)
    tmp_size <- subset(covermu_grow, Iteration==randiter &
                                Chain==randchain)
  }
  size_vec <- tmp_size$value
  intercept_vec <- tmp_intercept$value
  
  ##  Now do group, climate, and competition fixed effects
  # Group effects
  if(is.na(groupnum)==TRUE){
    tmp_group <- rep(0,length(spp_list))
    group_vec <- tmp_group
  }
  if(is.na(groupnum)==FALSE){
    tmp_group <- subset(group_grow, Iteration==randiter &
                               Chain==randchain &
                               groupid==groupnum)
    group_vec <- tmp_group$value
  }
  
  # Climate effects
  tmp_clim <- subset(climeff_grow, Iteration==randiter &
                              Chain==randchain)
  clim_mat <- matrix(tmp_clim$value, length(unique(tmp_clim$Parameter)), length(spp_list))
  
  # Crowding effect
  tmp_crowd <- subset(crowd_grow, Iteration==randiter &
                             Chain==randchain)
  crowd_mat <- matrix(tmp_crowd$value, length(unique(tmp_crowd$Parameter)), length(spp_list))
  
  ###TODO: subset the variance params...
  #Tau for size variance
  tmp_tau <- subset(tau_grow, Iteration==randiter &
                          Chain==randchain)
  tmp_tauSize <- subset(tauSize_grow, Iteration==randiter &
                                 Chain==randchain)
  tau_vec <- tmp_tau$value
  tauSize_vec <- tmp_tauSize$value
  
  ##  Collate all parameters for output
  Gpars=list(intcpt=intercept_vec, 
             intG=group_vec,
             slope=size_vec,
             nb=crowd_mat[1,],
             nbXsize=crowd_mat[2,],
             clim=clim_mat,
             sigma2.a=tau_vec,
             sigma2.b=tauSize_vec)
  return(Gpars)
}
