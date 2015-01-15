model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(gMu[i]) <- intYr[spp[i],yrs[i]] + beta[spp[i],yrs[i]]*log(X[i]) + temp1[spp[i]]*TmeanSpr1[i] + temp2[spp[i]]*TmeanSpr2[i] + rain1[spp[i]]*ppt1[i] + rain2[spp[i]]*ppt2[i]  
    log(tau[i]) <- intYrT[spp[i],yrs[i]] + betaT[spp[i],yrs[i]]*X[i]
    p[i] <- gMu[i] * tau[i]
    q[i] <- (1 - gMu[i]) * tau[i]
    C[i] ~ dbeta(p[i], q[i])
  }
  
  #priors
  for(s in 1:nSpp){
    betaSpp[s] ~ dnorm(0, 1e-6)
    temp1[s] ~ dnorm(temp1Mu, temp1Var)
    temp2[s] ~ dnorm(temp2Mu, temp2Var)
    rain1[s] ~ dnorm(rain1Mu, rain1Var)
    rain2[s] ~ dnorm(rain2Mu, rain2Var)
    intercept[s] ~ dnorm(0, 1e-6)
    intVaryY[s] ~ dgamma(0.001, 0.001)
    betaVar[s] ~ dgamma(0.001, 0.001)
    intVarG[s] ~ dgamma(2, 0.5) 
    
    betaSppT[s] ~ dnorm(0, 1e-6)
    interceptT[s] ~ dnorm(0, 1e-6)
    intVaryYT[s] ~ dgamma(0.001, 0.001)
    betaVarT[s] ~ dgamma(0.001, 0.001)
    for(y in 1:nYrs){
      intYr[s,y] ~ dnorm(intercept[s], intVaryY[s])
      beta[s,y] ~ dnorm(betaSpp[s], betaVar[s])
      intYrT[s,y] ~ dnorm(interceptT[s], intVaryYT[s])
      betaT[s,y] ~ dnorm(betaSppT[s], betaVarT[s])
    }
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }
  
  temp1Mu ~ dnorm(0,1e-6)
  temp2Mu ~ dnorm(0,1e-6)
  rain1Mu ~ dnorm(0,1e-6)
  rain2Mu ~ dnorm(0,1e-6)
  temp1Var ~ dgamma(0.001, 0.001)
  temp2Var ~ dgamma(0.001, 0.001)
  rain1Var ~ dgamma(0.001, 0.001)
  rain2Var ~ dgamma(0.001, 0.001)
}