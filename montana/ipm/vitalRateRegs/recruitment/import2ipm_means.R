
library(reshape2)
library(plyr)
fitthin <- data.frame(Parameter=NA, value=NA, species=NA)
for(ispp in spp_list){
  fitlong <- readRDS(paste("../vitalRateRegs/recruitment/recruitment_stanmcmc_", ispp, ".RDS", sep=""))
  tmp <- ddply(fitlong, .(Parameter), summarise,
               value = mean(value))
  tmp$species <- ispp
  fitthin <- rbind(fitthin,tmp)
}
fitthin <- fitthin[2:nrow(fitthin),]

##  Break up MCMC into regression components
# Climate effects
climeff_rec <- fitthin[grep("b2", fitthin$Parameter),]

# Yearly intercepts
intercept_rec <- fitthin[grep("a", fitthin$Parameter),]
intercept_rec <- subset(intercept_rec, Parameter!="a_mu")
intercept_rec <- subset(intercept_rec, Parameter!="tau")
intercept_rec$yearid <- substr(intercept_rec$Parameter, 3, length(intercept_rec$Parameter))
intercept_rec$yearid <- unlist(strsplit(intercept_rec$yearid, split=']'))

# Mean intercept
interceptmu_rec <- fitthin[grep("a_mu", fitthin$Parameter),]

# Group effects
group_rec <- fitthin[grep("gint", fitthin$Parameter),]
group_rec$groupid <- substr(group_rec$Parameter, 6, length(group_rec$Parameter))
group_rec$groupid <- unlist(strsplit(group_rec$groupid, split=']'))

# Get parent density-dependent effects
densdep_rec <- fitthin[grep("dd", fitthin$Parameter),]

# Get mixing fraction u
mixfrac_rec <- fitthin[grep("u", fitthin$Parameter),]
mixfrac_rec <- subset(mixfrac_rec, Parameter!="a_mu")

# Get theta
theta_rec <- fitthin[grep("theta", fitthin$Parameter),]

# Get rid of big objects
rm(list = c("tmp","fitthin","fitlong"))

##  Define function to format survival coefficients
getRecCoefs <- function(doYear, groupnum){
  
  # Get random effects if doYear!=NA
  if(is.na(doYear)==FALSE){
    tmp_intercept <- subset(intercept_rec, yearid==doYear)
  }
  
  # Set mean intercept and slope if doYear==NA
  if(is.na(doYear)==TRUE){
    tmp_intercept <-interceptmu_rec
  }
  intercept_vec <- tmp_intercept$value
  
  ##  Now do group, climate, and competition fixed effects
  # Group effects
  if(is.na(groupnum)==TRUE){
    tmp_group <- rep(0,length(spp_list))
    group_vec <- tmp_group
  }
  if(is.na(groupnum)==FALSE){
    tmp_group <- subset(group_rec, groupid==groupnum)
    group_vec <- tmp_group$value
  }
  
  # Climate effects
  tmp_clim <- climeff_rec
  clim_mat <- matrix(tmp_clim$value, length(unique(tmp_clim$Parameter)), length(spp_list))
  
  # Density dependence effect
  tmp_dd <- densdep_rec
  dd_vec <- tmp_dd$value
  
  # Mixing fraction
  tmp_u <- mixfrac_rec
  u_vec <- tmp_u$value
  
  # Theta
  tmp_theta <- theta_rec
  theta_vec <- tmp_theta$value
  
  ##  Collate all parameters for output
  Rpars=list(intcpt=intercept_vec, 
             grpInt=group_vec,
             dd=dd_vec,
             clim=clim_mat,
             u=u_vec,
             theta=theta_vec)
  return(Rpars)
}

