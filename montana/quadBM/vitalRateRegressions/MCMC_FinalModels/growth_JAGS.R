model{
  #growth process and likelihood
  for(i in 1:nObs){
    logit(gMu[i]) <- intercept + group[grp[i]] + year[yrs[i]] + size[i]*(beta+yearB[yrs[i]])
    p[i] <- gMu[i] * tau
    q[i] <- (1 - gMu[i]) * tau
    Y[i] ~ dbeta(p[i], q[i])
  }
  
  #priors
  for(g in 1:nGrp){
    group[g] ~ dnorm(0, .001)
  }
  for(y in 1:nYrs){
    year[y] ~ dnorm(0, .001)
    yearB[y] ~ dnorm(0, .001)
  }
  intercept ~ dnorm(0, .001)
  beta ~ dnorm(0, .001)
  t0 ~ dnorm(0, .01) 
  tau <- exp(t0)
  
}