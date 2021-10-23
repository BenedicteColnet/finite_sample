# Inputs are 10 independent variables uniformly distributed on the interval [0,1], 
# only 5 out of these 10 are actually used. Outputs are created according to the formula
# y = 10 sin(π x1 x2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5 + e
generate_simulation_friedman_constant_cate <- function(n_obs = 1000, p = 12, ATE = 3, independent_covariate = TRUE){
  
  # generate multivariate gaussian vector
  
  if(independent_covariate){
    X = rmvnorm(n = n_obs, mean = rep(0, p), sigma = diag(p))
  } else {
    cov_mat = toeplitz(0.7^(0:(p - 1)))
    X = rmvnorm(n = n_obs, mean = rep(0, p), sigma = cov_mat)
  }
  
  b = 10 * sin(pi*X[,2]*X[,4]) + 20* (X[,3] - 0.5)**2 + 10*X[,5] + 5*X[,6]
  e = 1/(1 + exp(-X[,1]) + exp(-X[,2]) + exp(-X[,3]))
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, e = e)
  simulation$Y_0 <- simulation$b + rnorm(n_obs)
  simulation$Y_1 <- simulation$b + ATE + rnorm(n_obs)
  simulation$A <- rbinom(n_obs, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  return(simulation)
}


# simulation set up according to Wager & Nie
generate_simulation_wager_nie <- function(n = 1000, p = 12, setup = "D"){
  
  # set-ups
  if (setup == "A"){
    X = matrix(runif(n*p, min=0, max=1), n, p)
    b = sin(pi * X[,2] * X[,3]) + 2 * (X[,4] - 0.5)^2 + X[,5] + 0.5 * X[,6]
    eta = 0.1
    e = pmax(eta, pmin(sin(pi * X[,1] * X[,2] * X[,3]), 1-eta))
    tau = (X[,3] + X[,4]) / 2
  } else if (setup == "B"){
    X = matrix(rnorm(n * p), n, p)
    b = pmax(0, X[,2] + X[,3], X[,4]) + pmax(0, X[,5] + X[,6])
    e = 0.5
    tau = X[,2] + log(1 + exp(X[,3]))
  } else if (setup == "C") {
    X = matrix(rnorm(n * p), n, p)
    b = 2 * log(1 + exp(X[,2] + X[,3] + X[,4]+ X[,5]+ X[,6]))
    e = 1/(1 + exp(X[,1] + X[,2] + X[,3]))
    tau = rep(3, n)
  } else if (setup == "D") {
    X = matrix(rnorm(n*p), n, p)
    b = (pmax(X[,2] + X[,3] + X[,4], 0) + pmax(X[,5] + X[,6], 0)) / 2
    e = 1/(1 + exp(-X[,1]) + exp(-X[,2] + exp(-X[,3])))
    tau = pmax(X[,2] + X[,3] + X[,4], 0) - pmax(X[,5] + X[,6], 0)
  } else if (setup == "simple"){
    X = matrix(rnorm(n*p), n, p)
    b = X[,2] 
    e = 1/(1 + exp(-X[,1]) + exp(-X[,2] + exp(-X[,3])))
    tau = X[,3] + X[,4] + X[,5] + X[,6]
  } else {
    print("error in setup")
    break
  }
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, tau = tau, e = e)
  simulation$Y_0 <- simulation$b - 0.5*simulation$tau + rnorm(n, mean = 0, sd = 0.1)
  simulation$Y_1 <- simulation$b + 0.5*simulation$tau + rnorm(n, mean = 0, sd = 0.1)
  simulation$A <- rbinom(n, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  return(simulation)
}

# # Empirically compute true ATEs
# a_big_simulation <- generate_simulation(n = 50000, setup = "A")
# ATE_A <- mean(a_big_simulation$Y_1) - mean(a_big_simulation$Y_0)
# 
# a_big_simulation <- generate_simulation(n = 50000, setup = "B")
# ATE_B <- mean(a_big_simulation$Y_1) - mean(a_big_simulation$Y_0)
# 
# a_big_simulation <- generate_simulation(n = 50000, setup = "C")
# ATE_C <- mean(a_big_simulation$Y_1) - mean(a_big_simulation$Y_0)
# 
# a_big_simulation <- generate_simulation(n = 50000, setup = "D")
# ATE_D <- mean(a_big_simulation$Y_1) - mean(a_big_simulation$Y_0)
# 
# ATE <- c("A" = ATE_A, 
#          "B" = ATE_B,
#          "C" = ATE_C,
#          "D" = ATE_D)


# if ATE = 3 the naive ATE is overestimated
generate_simulation_linear_constant_cate <- function(n_obs = 1000, p = 12, ATE = 3, independent_covariate = TRUE){
  
  # generate multivariate gaussian vector
  
  if(independent_covariate){
    X = rmvnorm(n = n_obs, mean = rep(0, p), sigma = diag(p))
  } else {
    cov_mat = toeplitz(0.7^(0:(p - 1)))
    X = rmvnorm(n = n_obs, mean = rep(0, p), sigma = cov_mat)
  }
  
  b = X[,4] + X[,5] + X[,6]+ X[,7]+ X[,8] +  X[,9] +  X[,10]
  e = 1/(1 + exp(-X[,1]) + exp(-X[,2]) + exp(-X[,3])+ exp(-X[,4])+ exp(-X[,5])+ exp(-X[,6]) + exp(-X[,7]))
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, e = e)
  simulation$Y_0 <- simulation$b + rnorm(n_obs)
  simulation$Y_1 <- simulation$b + ATE + rnorm(n_obs)
  simulation$A <- rbinom(n_obs, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  return(simulation)
}