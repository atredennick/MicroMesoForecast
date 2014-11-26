model{
  #process model and likelihood
  for(i in 1:nObs){
    mu[i] <- intercept[spp[i]] + intG[spp[i],grp[i]] + betaSpp[spp[i]]*X[i] + nb[spp[i]]*crowd[i]

    Y[i] ~ dnorm(mu[i], error[spp[i]])
  }
  
  #priors
  for(s in 1:nSpp){
    betaSpp[s] ~ dnorm(0, 1e-6)
    nb[s] ~ dnorm(0, 1e-6)
    intercept[s] ~ dnorm(0, 1e-6)
    error[s] ~ dgamma(0.001, 0.001)
    intVarG[s] ~ dgamma(2, 0.5) 
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }
  
}