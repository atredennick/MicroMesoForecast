model{
  #process model and likelihood
  for(i in 1:nObs){
    logit(p[i]) <- intercept[spp[i]] + intG[spp[i],grp[i]] + betaSpp[spp[i]]*X[i]
    Y[i] ~ dbern(p[i])
  }
  
  #priors
  for(s in 1:nSpp){
    betaSpp[s] ~ dnorm(0, 1e-6)
    intercept[s] ~ dnorm(0, 1e-6)
    intVarG[s] ~ dgamma(2, 0.5) 
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }

}