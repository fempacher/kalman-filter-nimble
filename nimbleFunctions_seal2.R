# Defining the 2 state Seal Model in nimbleCode

source("NIMBLE/helperFunctionsNimble.R")

nimbleMod_seal2 <- nimbleCode({
  # True model
  # Priors and constraints
  ###################################
  # Prior for initial population size
  # make init.pups/pups_survived 3rd coordinate
  #init.pups ~ dnorm(y0, tau/y0^2)
  x[1,3] ~ dnorm(y0, tau/y0^2)
  x[1,1] <- abs(round(x[1,3]))
  # make init.adults.nokids/adults_survived 4th coordinate
  #init.adults.nokids ~ dnegbin(alpha, x[1,1])
  x[1,4] ~ dnegbin(alpha, x[1,1])
  x[1,2] <- x[1,4] + x[1,1]
  ###################################
  # Prior for max pup survival
  phi_p ~ dbeta(2.87, 1.78)
  # Prior for adult survival
  # phi_a0 ~ dbeta(1.6,  1.2)
  # phi_a <- 0.8 + 0.17*phi_a0
  phi_a ~ dbetashift(1.6, 1.2, 0.8, 0.97)
  # Prior for fecundity
  # alpha0 ~ dbeta(2, 1.5)
  # alpha <- 0.6 + 0.4*alpha0
  alpha ~ dbetashift(2, 1.5, 0.6, 1)
  # Prior for CC parameter
  # chi ~ dgamma(4, 1250)
  chi ~ dgamma(4, 1/800)
  # prior for CC shape parameter
  rho ~ dgamma(4, 1/2.5)
  # convert chi to beta
  beta <- 1/chi*((alpha*0.5*phi_p)/(1-phi_a)-1)^(1/rho)
  # Prior for precision
  tau ~ dgamma(2.1, 1/66.67)
  # Likelihood
  ###################################
  # State process
  for (t in 1:(totalT)){
    real_phi_p[t] <- phi_p/(1+(beta*x[t,1])^rho)
    # survived pups is 3rd coordinate, survived adults is 4th coordinate
    # survived_pups[t] ~ dbin(0.5*real_phi_p[t], x[t,1])
    # survived_adults[t] ~ dbin(phi_a, x[t,2])
    x[t+1,3] ~ dbin(0.5*real_phi_p[t], x[t,1])
    x[t+1,4] ~ dbin(phi_a, x[t,2])
    x[t+1,2] <- x[t+1,3] + x[t+1,4]
    x[t+1,1] ~ dbin(alpha, x[t+1,2])
  }
  ###################################
  # Observation process
  for (t in 2:(totalT+1)) {
    y[t-1] ~ dnorm(x[t,1], tau/x[t,1]^2)
  }
})

nimbleMod_seal2_bkf <- nimbleCode({
  # Kalman filter of the NDLM approximation
  # Priors and constraints
  ###################################
  #######
  # Prior for max pup survival
  phi_p ~ dbeta(2.87, 1.78)
  # Prior for adult survival
  phi_a ~ dbetashift(1.6, 1.2, 0.8, 0.97)
  # Prior for fecundity
  alpha ~ dbetashift(2, 1.5, 0.6, 1)
  # Prior for CC parameter
  chi ~ dgamma(4, 1/800)
  # prior for CC shape parameter
  rho ~ dgamma(4, 1/2.5)
  # convert chi to beta
  beta <- 1/chi*((alpha*0.5*phi_p)/(1-phi_a)-1)^(1/rho)
  # Prior for precision
  tau ~ dgamma(2.1, 1/66.67)
  
  sigma2 <- y0^2/tau
  P_upd11[1] <- sigma2
  P_upd22[1] <- 1/(alpha^2)*((1-alpha)*y0+sigma2)
  P_upd12[1] <- 1/alpha*sigma2
  
  # Set matrices as far as possible
  T12 <- alpha*phi_a
  # Initialisation
  x_upd[1, 1:2] <-  y0*c(1, 1/alpha)
  
  for(t in 1:totalT){
    # Calculate matrices
    real_phi_p[t] <- phi_p/(1+(beta*x_upd[t, 1])^rho)
    T21[t] <- 0.5*real_phi_p[t]
    T11[t] <-  alpha*0.5*real_phi_p[t]
    
    # Prediction
    x_pred[t, 1] <- T11[t]*x_upd[t, 1] + T12*x_upd[t, 2]
    x_pred[t, 2] <- T21[t]*x_upd[t, 1] + phi_a*x_upd[t, 2]
    
    P_pred11[t] <- (T11[t]*P_upd11[t]+T12*P_upd12[t])*T11[t]+(T11[t]*P_upd12[t]+T12*P_upd22[t])*T12 +
      alpha*(1-alpha)*(0.5*real_phi_p[t]*x_upd[t, 1]+phi_a*x_upd[t, 2])
    P_pred12[t] <- (T11[t]*P_upd11[t]+T12*P_upd12[t])*T21[t]+(T11[t]*P_upd12[t]+T12*P_upd22[t])*phi_a +
      alpha*(0.5*real_phi_p[t]*(1-0.5*real_phi_p[t])*x_upd[t, 1] + phi_a*(1-phi_a)*x_upd[t, 2])
    P_pred22[t] <- (T21[t]*P_upd11[t]+phi_a*P_upd12[t])*T21[t]+(T21[t]*P_upd12[t]+phi_a*P_upd22[t])*phi_a +
      0.5*real_phi_p[t]*(1-0.5*real_phi_p[t])*x_upd[t, 1] + phi_a*(1-phi_a)*x_upd[t, 2]
        
    # Aid variable
    v_pred[t] <- y[t] - x_pred[t, 1]
    F_matrix[t] <- P_pred11[t] + (x_pred[t, 1]^2 + P_pred11[t])/tau
    
    # Distribution of y
    y[t] ~ dnorm(x_pred[t, 1], var = F_matrix[t])
    
    F_inv[t] <- F_matrix[t]^(-1)
    
    # Update
    x_upd[t+1, 1] <- x_pred[t, 1] + P_pred11[t]*F_inv[t]*v_pred[t]
    x_upd[t+1, 2] <- x_pred[t, 2] + P_pred12[t]*F_inv[t]*v_pred[t]
    P_upd11[t+1] <- P_pred11[t] - F_inv[t]*P_pred11[t]^2
    P_upd12[t+1] <- P_pred12[t] - F_inv[t]*P_pred11[t]*P_pred12[t]
    P_upd22[t+1] <- P_pred22[t] - F_inv[t]*P_pred12[t]^2
  }
})

