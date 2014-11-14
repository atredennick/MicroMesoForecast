model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(gMu[i]) <- intercept[spp[i]] + intG[spp[i],grp[i]] + betaSpp[spp[i]]*X[i]
    p[i] <- gMu[i] * tau[spp[i]]
    q[i] <- (1 - gMu[i]) * tau[spp[i]]
    C[i] ~ dbeta(p[i], q[i])
  }
  
  #priors
  for(s in 1:nSpp){
    betaSpp[s] ~ dnorm(0, 1e-6)
    intercept[s] ~ dnorm(0, 1e-6)
    t0[s] ~ dnorm(0, .01) 
    tau[s] <- exp(t0[s])
    intVarG[s] ~ dgamma(2, 0.5) 
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }

}