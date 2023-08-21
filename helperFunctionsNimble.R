# Defining the 2 state Seal Model in nimbleCode
# x as state, not pups and adults separately

library(nimble)

# Define new nimble distribution: Shifted beta distribution
dbetashift <- nimbleFunction(
  run = function(x = double(0),
                 shape1 = double(0, default = 1),
                 shape2 = double(0, default = 1),
                 min = double(0, default = 0),
                 max = double(0, default = 1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    x_scaled <- (x-min)/(max-min)
    logProb <- dbeta(x_scaled, shape1 = shape1, shape2=shape2, log = TRUE) - log(max-min)
    if(log) return(logProb)
    else return(exp(logProb)) 
  })

rbetashift <- nimbleFunction(
  run = function(n = integer(0),
                 shape1 = double(0, default = 1),
                 shape2 = double(0, default = 1),
                 min = double(0, default = 0),
                 max = double(0, default = 1)) {
    returnType(double(0))
    if(n != 1) print("rbetashift only allows n = 1; using n = 1.")
    dev_scaled <- rbeta(n=1, shape1=shape1,shape2 = shape2)
    dev <- dev_scaled*(max-min)+min
    return(dev)
  })

# Define new nimble distribution: Shifted gamma distribution
dgammashift <- nimbleFunction(
  run = function(x = double(0),
                 shape = double(0, default = 1),
                 rate = double(0, default = 1),
                 shift = double(0, default = 0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    x_scaled <- x-shift
    logProb <- dgamma(x_scaled, shape = shape, rate = rate, log = TRUE)
    if(log) return(logProb)
    else return(exp(logProb)) 
  })

rgammashift <- nimbleFunction(
  run = function(n = integer(0),
                 shape = double(0, default = 1),
                 rate = double(0, default = 1),
                 shift = double(0, default = 0)) {
    returnType(double(0))
    if(n != 1) print("rgammashift only allows n = 1; using n = 1.")
    dev_scaled <- rgamma(n=1, shape=shape, rate = rate)
    dev <- dev_scaled + shift
    return(dev)
  })


# Function to do one step in the Kalman filter

oneStepKalman <- nimbleFunction(
  run = function(x_upd = double(1), P_upd = double(2),
                 y = double(1),
                 phi_p = double(0),  phi_a = double(0), alpha = double(0),
                 betaNS = double(0), betaIH = double(0), betaOH = double(0),
                 betaOrk = double(0), rho = double(0), tau = double(0),
                 omega = double(0)) { # type declarations
    
    real_phi_pNS <- phi_p/(1+(betaNS*x_upd[1])^rho)
    real_phi_pIH <- phi_p/(1+(betaIH*x_upd[8])^rho)
    real_phi_pOH <- phi_p/(1+(betaOH*x_upd[15])^rho)
    real_phi_pOrk <- phi_p/(1+(betaOrk*x_upd[22])^rho)
    
    QNS[1,1:7] <- c(phi_a*alpha*(1-phi_a*alpha)*(x_upd[6]+x_upd[7]),0,0,0,
                      0,0,alpha*(x_upd[6]+x_upd[7])*phi_a*(1-phi_a))
    QNS[2,1:7] <- c(0,x_upd[1]*0.5*real_phi_pNS*(1-0.5*real_phi_pNS),0,0,0,0,0)
    QNS[3,1:7] <- c(0,0,x_upd[2]*phi_a*(1-phi_a),0,0,0,0)
    QNS[4,1:7] <- c(0,0,0,x_upd[3]*phi_a*(1-phi_a),0,0,0)
    QNS[5,1:7] <- c(0,0,0,0,x_upd[4]*phi_a*(1-phi_a),0,0)
    QNS[6,1:7] <- c(0,0,0,0,0,x_upd[5]*phi_a*(1-phi_a),0)
    QNS[7,1:7] <- c(alpha*(x_upd[6]+x_upd[7])*phi_a*(1-phi_a),0,0,0,0,
                      0,(x_upd[6]+x_upd[7])*phi_a*(1-phi_a))
    QIH[1,1:7] <- c(phi_a*alpha*(1-phi_a*alpha)*(x_upd[13]+x_upd[14]),0,0,0,
                      0,0,alpha*(x_upd[13]+x_upd[14])*phi_a*(1-phi_a))
    QIH[2,1:7] <- c(0,x_upd[8]*0.5*real_phi_pIH*(1-0.5*real_phi_pIH),0,0,0,0,0)
    QIH[3,1:7] <- c(0,0,x_upd[9]*phi_a*(1-phi_a),0,0,0,0)
    QIH[4,1:7] <- c(0,0,0,x_upd[10]*phi_a*(1-phi_a),0,0,0)
    QIH[5,1:7] <- c(0,0,0,0,x_upd[11]*phi_a*(1-phi_a),0,0)
    QIH[6,1:7] <- c(0,0,0,0,0,x_upd[12]*phi_a*(1-phi_a),0)
    QIH[7,1:7] <- c(alpha*(x_upd[13]+x_upd[14])*phi_a*(1-phi_a),0,0,0,0,
                      0,(x_upd[13]+x_upd[14])*phi_a*(1-phi_a))
    QOH[1,1:7] <- c(phi_a*alpha*(1-phi_a*alpha)*(x_upd[20]+x_upd[21]),0,0,0,
                      0,0,alpha*(x_upd[20]+x_upd[21])*phi_a*(1-phi_a))
    QOH[2,1:7] <- c(0,x_upd[15]*0.5*real_phi_pOH*(1-0.5*real_phi_pOH),0,0,0,0,0)
    QOH[3,1:7] <- c(0,0,x_upd[16]*phi_a*(1-phi_a),0,0,0,0)
    QOH[4,1:7] <- c(0,0,0,x_upd[17]*phi_a*(1-phi_a),0,0,0)
    QOH[5,1:7] <- c(0,0,0,0,x_upd[18]*phi_a*(1-phi_a),0,0)
    QOH[6,1:7] <- c(0,0,0,0,0,x_upd[19]*phi_a*(1-phi_a),0)
    QOH[7,1:7] <- c(alpha*(x_upd[20]+x_upd[21])*phi_a*(1-phi_a),0,0,0,0,
                      0,(x_upd[20]+x_upd[21])*phi_a*(1-phi_a))
    QOrk[1,1:7] <- c(phi_a*alpha*(1-phi_a*alpha)*(x_upd[27]+x_upd[28]),0,0,0,
                       0,0,alpha*(x_upd[27]+x_upd[28])*phi_a*(1-phi_a))
    QOrk[2,1:7] <- c(0,x_upd[22]*0.5*real_phi_pOrk*(1-0.5*real_phi_pOrk),0,0,0,0,0)
    QOrk[3,1:7] <- c(0,0,x_upd[23]*phi_a*(1-phi_a),0,0,0,0)
    QOrk[4,1:7] <- c(0,0,0,x_upd[24]*phi_a*(1-phi_a),0,0,0)
    QOrk[5,1:7] <- c(0,0,0,0,x_upd[25]*phi_a*(1-phi_a),0,0)
    QOrk[6,1:7] <- c(0,0,0,0,0,x_upd[26]*phi_a*(1-phi_a),0)
    QOrk[7,1:7] <- c(alpha*(x_upd[27]+x_upd[28])*phi_a*(1-phi_a),0,0,0,0,
                       0,(x_upd[27]+x_upd[28])*phi_a*(1-phi_a))
    
    Q[1:7,1:7] <- QNS[1:7,1:7]
    Q[8:14,8:14] <- QIH[1:7,1:7]
    Q[15:21,15:21] <- QOH[1:7,1:7]
    Q[22:28,22:28] <- QOrk[1:7,1:7]
    for(i in 1:7){
      Q[8:28,i] <- rep(0, times = 21)
      Q[1:7,i+7] <- rep(0, times = 7)
      Q[15:28,i+7] <- rep(0, times = 14)
      Q[1:14,i+14] <- rep(0, times = 14)
      Q[22:28,i+14] <- rep(0, times = 7)
      Q[1:21,i+21] <- rep(0, times = 21)
    }
    
    
    T_matrix[1,1:7] <- c(0,0,0,0,0,phi_a*alpha, phi_a*alpha)
    T_matrix[2,1:7] <- c(0.5*real_phi_pNS,0,0,0,0,0,0)
    T_matrix[3,1:7] <- c(0,phi_a,0,0,0,0,0)
    T_matrix[4,1:7] <- c(0,0,phi_a,0,0,0,0)
    T_matrix[5,1:7] <- c(0,0,0,phi_a,0,0,0)
    T_matrix[6,1:7] <- c(0,0,0,0,phi_a,0,0)
    T_matrix[7,1:7] <- c(0,0,0,0,0,phi_a,phi_a)
    T_matrix[8,8:14] <- c(0,0,0,0,0,phi_a*alpha, phi_a*alpha)
    T_matrix[9,8:14] <- c(0.5*real_phi_pIH,0,0,0,0,0,0)
    T_matrix[10,8:14] <- c(0,phi_a,0,0,0,0,0)
    T_matrix[11,8:14] <- c(0,0,phi_a,0,0,0,0)
    T_matrix[12,8:14] <- c(0,0,0,phi_a,0,0,0)
    T_matrix[13,8:14] <- c(0,0,0,0,phi_a,0,0)
    T_matrix[14,8:14] <- c(0,0,0,0,0,phi_a,phi_a)
    T_matrix[15,15:21] <- c(0,0,0,0,0,phi_a*alpha, phi_a*alpha)
    T_matrix[16,15:21] <- c(0.5*real_phi_pOH,0,0,0,0,0,0)
    T_matrix[17,15:21] <- c(0,phi_a,0,0,0,0,0)
    T_matrix[18,15:21] <- c(0,0,phi_a,0,0,0,0)
    T_matrix[19,15:21] <- c(0,0,0,phi_a,0,0,0)
    T_matrix[20,15:21] <- c(0,0,0,0,phi_a,0,0)
    T_matrix[21,15:21] <- c(0,0,0,0,0,phi_a,phi_a)
    T_matrix[22,22:28] <- c(0,0,0,0,0,phi_a*alpha, phi_a*alpha)
    T_matrix[23,22:28] <- c(0.5*real_phi_pOrk,0,0,0,0,0,0)
    T_matrix[24,22:28] <- c(0,phi_a,0,0,0,0,0)
    T_matrix[25,22:28] <- c(0,0,phi_a,0,0,0,0)
    T_matrix[26,22:28] <- c(0,0,0,phi_a,0,0,0)
    T_matrix[27,22:28] <- c(0,0,0,0,phi_a,0,0)
    T_matrix[28,22:28] <- c(0,0,0,0,0,phi_a,phi_a)
    for(i in 1:7){
      T_matrix[8:28,i] <- rep(0, times = 21)
      T_matrix[1:7,i+7] <- rep(0, times = 7)
      T_matrix[15:28,i+7] <- rep(0, times = 14)
      T_matrix[1:14,i+14] <- rep(0, times = 14)
      T_matrix[22:28,i+14] <- rep(0, times = 7)
      T_matrix[1:21,i+21] <- rep(0, times = 21)
    }
    
    R[1,1:4] <- c(x_upd[1]^2/tau,0,0,0)
    R[2,1:4] <- c(0,x_upd[8]^2/tau,0,0)
    R[3,1:4] <- c(0,0,x_upd[15]^2/tau,0)
    R[4,1:4] <- c(0,0,0,x_upd[22]^2/tau)
    
    # Prediction
    x_pred[1:28] <- T_matrix[1:28, 1:28, t]%*%x_upd[1:28]
    P_pred[1:28,1:28, t] <- T_matrix[1:28,1:28]%*%P_upd[1:28,1:28]%*%t(T_matrix[1:28,1:28]) + Q[1:28,1:28]
    
    y_pred[1:4] <-Z[1:4,1:28]%*%x_pred[1:28]
    F_matrix[1:4, 1:4] <- Z[1:4,1:28]%*%P_pred[1:28,1:28]%*%t(Z[1:4,1:28]) + R[1:4,1:4]
    F_inv[1:4,1:4] <- inverse(F_matrix[1:4, 1:4])
    
    v_pred[1:4] <- y[1:4] - y_pred[1:4]
    
    # Update
    G[1:28,1:4] <- P_pred[1:28,1:28]%*%t(Z[1:4,1:28])
    K[1:28,1:4] <- G[1:28,1:4]%*%F_inv[1:4,1:4]
    
    x[1:28] <- x_pred[1:28] + K[1:28,1:4]%*%v_pred[1:4]
    P[1:28, 1:28] <- P_pred[1:28,1:28] - G[1:28,1:4]%*%F_inv[1:4,1:4]%*%t(G[1:28,1:4])
    
    exampleNimListDef <- nimbleList(x = double(1), P = double(2),
                                    y_pred = double(1), F_matrix = double(2))
    exampleNimList <- exampleNimListDef$new(x = x, P = P, y_pred  = y_pred,
                                            F_matrix = F_matrix)
    return(exampleNimList)
    returnType(exampleNimListDef())  # return type declaration
  } )












