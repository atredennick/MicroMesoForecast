model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(p[i]) <- intercept[spp[i]] + intG[grp[i]] + betaSpp[spp[i]]*X[i]
    C[i] ~ dbern(p[i])
  }
  
  #priors
  for(s in 1:nSpp){
    betaSpp[s] ~ dnorm(0, 1e-6)
    intercept[s] ~ dnorm(0, 1e-6)
  }
  
  for(g in 1:nGrp){
    intG[g] ~ dnorm(0, intVarG)
  }
  
  intVarG ~ dgamma(2, 0.5) 
}