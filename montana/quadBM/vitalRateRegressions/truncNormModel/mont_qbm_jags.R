mont_qbm_jags <- function(Y,size,X,groups,species,years,iters,adapt.iters,thins){
  ####
  ####  Setup variables
  ####
  
  X <- as.matrix(X)
  n_covs <- dim(X)[2]
  n <- length(Y) #get number of observations
  n_group <- length(unique(groups)) #get number of unique groups
  n_years <- length(unique(years)) #get number of unique years
  
  #Now write list for JAGS
  data_jags <- list(X=X, n_covs=n_covs, n=n, years=years, size=size,
                    n_group=n_group, n_years=n_years, Y=Y, groups=groups)
  
  ####
  ####  Write the model as a string for JAGS
  ####
  
  modelstring="
  model{
    #process and likelihood
    for(i in 1:n){
      mu[i] <- beta.0[years[i]] + beta.group[groups[i]] + beta.size[years[i]]*size[i] + inprod(beta[1:n_covs],X[i,1:n_covs])
      Y[i] ~ dlnorm(mu[i], tau)T(0,1)
    }
    
    #priors
    for(c in 1:n_covs) {
      beta[c] ~ dnorm(0, 1e-6)
    }
    for(y in 1:n_years){
      beta.0[y] ~ dnorm(beta.0.mu, beta.0.tau)
      beta.size[y] ~ dnorm(beta.size.mu, beta.size.tau)
    }
    for(g in 1:n_group){
      beta.group[g] ~ dnorm(0,beta.group.tau)
    } 
    beta.size.tau ~ dgamma(0.001, 0.001)
    beta.0.tau ~ dgamma(0.001, 0.001)
    beta.group.tau ~ dgamma(2, 0.5)
    beta.0.mu ~ dnorm(0, 1e-6)
    beta.size.mu ~ dnorm(0, 1e-6)
    tau ~ dgamma(0.001, 0.001)
  } #end model
  "
  
  
  ####
  ####  Call JAGS and fit model
  ####
  iterations <- iters
  adapt <- adapt.iters
  fit_model <- jags.model(textConnection(modelstring), data=data_jags, n.chains=3, n.adapt=adapt)
  update(fit_model, n.iter = (iterations*0.25))
  dic_model <- dic.samples(fit_model, n.iter = iterations, thin = thins, type = "pD")
  return(dic_model)
  
}#end of function
  