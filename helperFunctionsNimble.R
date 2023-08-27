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










