##  R Script to Import and Format Fitted QBM Parameters

fetch_qbm_params <- function(do_species, dir_to_find, mean_params){
  if(mean_params==TRUE){
    ##  Load vital rate parameters
    fitlong <- readRDS(paste0(dir_to_find, "popgrowth_stanmcmc_",do_species, ".RDS"))
    fitthin <- ddply(fitlong, .(Parameter), summarise,
                     value = mean(value))
    ##  Break up MCMC into regression components
    # Climate effects
    climeff <- fitthin[grep("b2", fitthin$Parameter),]
    climeff$id <- substr(climeff$Parameter, 4, length(climeff$Parameter))
    climeff$id <- unlist(strsplit(climeff$id, split=']'))
    climeff <- climeff[with(climeff, order(as.numeric(id))),]
    # Mean cover (size) effects
    coveff <- fitthin[grep(glob2rx("b1_mu"), fitthin$Parameter),]
    # Mean intercept
    intercept <- fitthin[grep("a_mu", fitthin$Parameter),]
    # Lognormal sigma (called tau here)
    tau <- fitthin[grep("tau", fitthin$Parameter),]
  } # end IF for mean parameters
  
  if(mean_params==FALSE){
    fitthin <- fitlong
    ##  Break up MCMC into regression components
    # Climate effects
    climeff <- fitthin[grep("b2", fitthin$Parameter),]
    # Mean cover (size) effects
    coveff <- fitthin[grep(glob2rx("b1_mu"), fitthin$Parameter),]
    # Mean intercept
    intercept <- fitthin[grep("a_mu", fitthin$Parameter),]
    # Lognormal sigma (called tau here)
    tau <- fitthin[grep("tau", fitthin$Parameter),]
  } # end IF for full MCMC param output
  
  out_list <- list(climeff=climeff, coveff=coveff, intercept=intercept, tau=tau)
  return(out_list)
}
