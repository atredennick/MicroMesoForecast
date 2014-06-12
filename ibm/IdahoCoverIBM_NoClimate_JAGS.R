model{
  #Priors
  varObsf ~ dgamma(0.01,0.01)
  varPf ~ dgamma(0.01,0.01)
  varObs <- 1/varObsf
  varP <- 1/varPf
  
  for(i in 1:3){
    regVarf[i] ~ dgamma(.01,.01)
  }
  regVar <- 1/regVarf
  
  for(i in 1:2){
    alphaMean[i] ~ dnorm(0,100)
    betaMean[i] ~ dnorm(0,100)
    varAlphaf[i] ~ dgamma(.01,.01)
    varBetaf[i] ~ dgamma(.01,.01)
  }
  varAlpha <- 1/varAlphaf
  varBeta <- 1/varBetaf
  
  for(i in 1:1){
    gammaMean[i] ~ dnorm(0,100)
    varGammaf[i] ~ dgamma(.01,.01)
  }
  varGamma <- 1/varGammaf
  
  for(i in 1:nQuad){
    for(j in 1:2){
      alpha[i,j] ~ dnorm(alphaMean[j], varAlpha[j])
      beta[i,j] ~ dnorm(betaMean[j], varBeta[j])
    }
    
    for(j in 1:1){
      gamma[i,j] ~ dnorm(gammaMean[j], varGamma[j])
    }
    N0[i] ~ dunif(0,100) #initial value for N, it has to be between 0-100%
  }
  
  #Process models
  #initial conditions
  for(k in 1:nQuad){
    S[1,k] <- 1/(1+exp(alpha[k,1]+alpha[k,2]*N0[k]))
    G[1,k] <- beta[k,1]+beta[k,2]*N0[k]
    C[1,k] <- 1/(1+exp(gamma[k,1]))
    
    Sp[1,k] ~ dbern(S[1,k])
    Gp[1,k] ~ dnorm(G[1,k], regVar[2])
    Cp[1,k] ~ dbern(C[1,k])
    
    N.pm[1,k] <- S[1,k]*G[1,k]
    N.p[1,k] <- ifelse(N.pm[1,k]==0, C[1,k]*0.01, N.pm[1,k])
    
    N[1,k] ~ dlnorm(N.p[1,k], varP)
  }
  
  for(t in 2:nProc){
    for(k in 1:nQuad){
      S[t,k] <- 1/(1+exp(alpha[k,1]+alpha[k,2]*N[t-1,k]))
      G[t,k] <- beta[k,1]+beta[k,2]*N[t-1,k]
      C[t,k] <- 1/(1+exp(gamma[k,1]))
      
      Sp[t,k] ~ dbern(S[t,k])
      Gp[t,k] ~ dnorm(G[t,k], regVar[2])
      Cp[t,k] ~ dbern(C[t,k])
      
      N.pm[t,k] <- S[t,k]*G[t,k]
      N.p[t,k] <- ifelse(N.pm[t,k]==0, C[t,k]*0.01, N.pm[t,k])
      
      N[t,k] ~ dlnorm(N.p[t,k], varP)
    }
  }
  
  #Likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(N[timeN[i],quadN[i]], varObs)
  }
  
  #Derived quantities
  #calculate mean over quads for each year
  qCount <- 1
  for(i in 1:nProc){
    muN[i] <- mean(N[i,])
  }
}