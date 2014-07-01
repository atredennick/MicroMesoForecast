model{
  #Priors
  varObsf ~ dgamma(0.01,0.01)
  varPf ~ dgamma(0.01,0.01)
  varObs <- 1/varObsf
  varP <- 1/varPf
  
  regVarf ~ dgamma(.01,.01)
  regVar <- 1/regVarf
  
  #Random intercept and size effects
  for(i in 1:2){
    alphaMean[i] ~ dnorm(0,100)
    betaMean[i] ~ dnorm(0,100)
    varAlphaf[i] ~ dgamma(.01,.01)
    varBetaf[i] ~ dgamma(.01,.01)
  }
  varAlpha <- 1/varAlphaf
  varBeta <- 1/varBetaf
  
  #Random intercept for colonization
  gammaMean ~ dnorm(0,100)
  varGammaf ~ dgamma(.01,.01)
  varGamma <- 1/varGammaf
  
  #Climate effects
  for(i in 1:4){
    alphaClim[i] ~ dnorm(0,100)
    betaClim[i] ~ dnorm(0,100)
    gammaClim[i] ~ dnorm(0,100)
  }
  
  #Draw random quadrat coefficients from mean distributions
  for(i in 1:nQuad){
    for(j in 1:2){
      alpha[i,j] ~ dnorm(alphaMean[j], varAlpha[j])
      beta[i,j] ~ dnorm(betaMean[j], varBeta[j])
    }
    gamma[i] ~ dnorm(gammaMean, varGamma)
    N0[i] ~ dunif(0,100) #initial value for N, it has to be between 0-100%
  }
  
  #Process models
  #initial conditions
  for(k in 1:nQuad){
    logit(S[1,k]) <- alpha[k,1]+alpha[k,2]*N0[k]+alphaClim[1]*clim[1,1]+alphaClim[2]*clim[1,2]+alphaClim[3]*clim[1,3]+alphaClim[4]*clim[1,4]
    G[1,k] <- beta[k,1]+beta[k,2]*N0[k]+betaClim[1]*clim[1,1]+betaClim[2]*clim[1,2]+betaClim[3]*clim[1,3]+betaClim[4]*clim[1,4]
    logit(C[1,k]) <- gamma[k]+gammaClim[1]*clim[1,1]+gammaClim[2]*clim[1,2]+gammaClim[3]*clim[1,3]+gammaClim[4]*clim[1,4]
    
#     Sp[1,k] ~ dbern(S[1,k])
#     Gp[1,k] ~ dnorm(G[1,k], regVar)
#     Cp[1,k] ~ dbern(C[1,k])
#     
    N.pm[1,k] <- S[1,k]*G[1,k]
    N.p[1,k] <- ifelse(N.pm[1,k]==0, C[1,k]*0.01, N.pm[1,k])
    
    N[1,k] ~ dnorm(N.p[1,k], varP)
  }
  
  for(t in 2:nProc){
    for(k in 1:nQuad){
      logit(S[t,k]) <- alpha[k,1]+alpha[k,2]*N[t-1,k]+alphaClim[1]*clim[t,1]+alphaClim[2]*clim[t,2]+alphaClim[3]*clim[t,3]+alphaClim[4]*clim[t,4]
      G[t,k] <- beta[k,1]+beta[k,2]*N[t-1,k]+betaClim[1]*clim[1,1]+betaClim[2]*clim[t,2]+betaClim[3]*clim[t,3]+betaClim[4]*clim[t,4]
      logit(C[t,k]) <- gamma[k]+gammaClim[1]*clim[t,1]+gammaClim[2]*clim[t,2]+gammaClim[3]*clim[t,3]+gammaClim[4]*clim[t,4]
      
#       Sp[t,k] ~ dbern(S[t,k])
#       Gp[t,k] ~ dnorm(G[t,k], regVar)
#       Cp[t,k] ~ dbern(C[t,k])
      
      N.pm[t,k] <- S[t,k]*G[t,k]
      N.p[t,k] <- ifelse(N.pm[t,k]==0, C[t,k]*0.01, N.pm[t,k])
      
      N[t,k] ~ dnorm(N.p[t,k], varP)
    }
  }

  #Likelihoods
  for(i in 1:nObs){
    y[i] ~ dnorm(N[timeN[i],quadN[i]], varObs)
  }
  for(i in 1:nCol){
    yCol[i] ~ dbern(C[timeNcol[i], quadNcol[i]])
  }
  for(i in 1:nSurv){
    ySurv[i] ~ dbern(S[timeNsurv[i], quadNsurv[i]])
  }
  for(i in 1:nGrow){
    yGrow[i] ~ dnorm(G[timeNgrow[i], quadNgrow[i]], regVar)
  }

#Derived quantities
  #calculate mean over quads for each year
  for(i in 1:nProc){
    muN[i] <- mean(N[i,])
  }

  # Assess model fit using a sum-of-squares-type discrepancy
  for (i in 1:nObs) {
    predicted[i] <- N[timeN[i],quadN[i]]            # Predicted values
    residual[i] <- y[i]-predicted[i]                # Residuals for observed data                                     
    sq[i] <- pow(residual[i], 2)                    # Squared residuals
  
    # Generate replicate data and compute fit statistics for them
    y.new[i] ~ dnorm(N[timeN[i],quadN[i]], varObs)        # One new data set at each MCMC iteration
    sq.new[i] <- pow(y.new[i]-predicted[i], 2)            # Squared residuals for new data
  }

  fit <- sum(sq[])              # Sum of squared residuals for actual data set
  fit.new <- sum(sq.new[])      # Sum of squared residuals for new data set
  test <- step(fit.new-fit)     # Test whether new data set more extreme
  bpvalue <- mean(test) 

}
  
#   #Posterior predictive distribution
#   for(k in 1:nQuad){
#     Snew[1,k] <- 1/(1+exp(alpha[k,1]+alpha[k,2]*N0[k]+alphaClim[1]*clim[1,1]+alphaClim[2]*clim[1,2]+alphaClim[3]*clim[1,3]+alphaClim[4]*clim[1,4]))
#     Gnew[1,k] <- beta[k,1]+beta[k,2]*N0[k]+betaClim[1]*clim[1,1]+betaClim[2]*clim[1,2]+betaClim[3]*clim[1,3]+betaClim[4]*clim[1,4]
#     Cnew[1,k] <- 1/(1+exp(gamma[k]+gammaClim[1]*clim[1,1]+gammaClim[2]*clim[1,2]+gammaClim[3]*clim[1,3]+gammaClim[4]*clim[1,4]))
#     
#     Spnew[1,k] ~ dbern(Snew[1,k])
#     Gpnew[1,k] ~ dnorm(Gnew[1,k], regVar)
#     Cpnew[1,k] ~ dbern(Cnew[1,k])
#     
#     N.pmnew[1,k] <- Snew[1,k]*Gnew[1,k]
#     N.pnew[1,k] <- ifelse(N.pmnew[1,k]==0, Cnew[1,k]*0.01, N.pmnew[1,k])
#     
#     Nnew[1,k] ~ dlnorm(N.pnew[1,k], varP)
#   }
#   
#   for(t in 2:nProc){
#     for(k in 1:nQuad){
#       Snew[t,k] <- 1/(1+exp(alpha[k,1]+alpha[k,2]*N[t-1,k]+alphaClim[1]*clim[t,1]+alphaClim[2]*clim[t,2]+alphaClim[3]*clim[t,3]+alphaClim[4]*clim[t,4]))
#       Gnew[t,k] <- beta[k,1]+beta[k,2]*N[t-1,k]+betaClim[1]*clim[1,1]+betaClim[2]*clim[t,2]+betaClim[3]*clim[t,3]+betaClim[4]*clim[t,4]
#       Cnew[t,k] <- 1/(1+exp(gamma[k]+gammaClim[1]*clim[t,1]+gammaClim[2]*clim[t,2]+gammaClim[3]*clim[t,3]+gammaClim[4]*clim[t,4]))
#       
#       Spnew[t,k] ~ dbern(Snew[t,k])
#       Gpnew[t,k] ~ dnorm(Gnew[t,k], regVar)
#       Cpnew[t,k] ~ dbern(Cnew[t,k])
#       
#       N.pmnew[t,k] <- Snew[t,k]*Gnew[t,k]
#       N.pnew[t,k] <- ifelse(N.pmnew[t,k]==0, Cnew[t,k]*0.01, N.pmnew[t,k])
#       
#       Nnew[t,k] ~ dlnorm(N.pnew[t,k], varP)
#     }
#   }
#   
#   for(t in 1:nProc){
#     for(k in 1:nQuad){
#       y.new[t,k] ~ dnorm(Nnew[t,k], varObs)
#     }
#   }
