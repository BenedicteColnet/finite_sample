generate_simulation <- function(n_obs = 1000, p = 12, all_covariates_output = FALSE, independent_covariate = FALSE){
  
  # generate multivariate gaussian vector
  if(independent_covariate){
    cov_mat = diag(p)
  } else {
    cov_mat = toeplitz(0.6^(0:(p - 1)))
  }
  
  X = matrix(rnorm(n_obs*p), n_obs, p)
  b = 2 * log(1 + exp(X[,1] + X[,2] + X[,3] + X[,4])) + X[,6:p]%*%rep(3,p-6+1) + X[,p]*X[,p]
  
  e = ifelse(X[,4] > 0, 1/(1 + exp( 1 -X[,3] - X[,2])), 1/(1 + exp(-X[,1] -X[,2])))
  
  tau = pmax(X[,3] + X[,5], 0)
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, tau = tau, e = e)
  simulation$mu_0 <- simulation$b - 0.5*simulation$tau 
  simulation$mu_1 <- simulation$b + 0.5*simulation$tau 
  simulation$Y_0 <- simulation$b - 0.5*simulation$tau + rnorm(n_obs, mean = 0, sd = 0.1)
  simulation$Y_1 <- simulation$b + 0.5*simulation$tau + rnorm(n_obs, mean = 0, sd = 0.1)
  simulation$A <- rbinom(n_obs, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }
  return(simulation)
}



# simulation set up according to Wager & Nie
generate_simulation_wager_nie <- function(n = 1000, p = 12, setup = "D", all_covariates_output = FALSE){
  
  # set-ups
  if (setup == "A"){
    X = matrix(runif(n*p, min=0, max=1), n, p)
    b = 3*sin(pi * X[,1] * X[,2]) + 4 * (X[,3] - 0.5)^2 -5*cos(pi * X[,4] * X[,5]) + 2 * (X[,6] - 0.5)^2
    + ifelse(X[,7] > 0, 3, -3) 
    + ifelse(X[,8] > 0 & X[,9] < 0, 10, 20)
    + 3*pmax(X[,10], 0.4)
    + 0.4*log(1+X[,11]+X[,12])
    eta = 0.1
    e = pmax(0.5*X[,6], pmin(sin(pi * X[,1] * X[,2] * X[,3] * X[,4] * X[,5]), 1-eta))
    tau = (X[,1] + X[,2]) / 2 + (1+1.4*X[,7])^2
  
  } else if (setup == "B"){
    
    X = matrix(rnorm(n * p), n, p)
    b = pmax(0, X[,1] + X[,2]) + pmax(0, X[,3] + X[,4]) + 0.5 * X[,5] + 0.4 * X[,6]^3 + 4*sin(pi * X[,7] * X[,8]) + X[,9] - X[,10]
    e = 1/(1 + exp(-2 - 0.8*(-X[,1] - X[,2] - X[,3] - X[,4] - X[,5] - X[,6])))
    tau = X[,11] + log(1 + exp(X[,12]))
    
  } else if (setup == "C") {
    
    X = matrix(rnorm(n * p), n, p)
    b = 2 * log(1 + exp(X[,1] + X[,2])) + X[,3] + X[,4] + 0.5 * X[,5] + 2 * X[,6] + X[,7]*X[,8] + log(pmax(1, X[,9]) + pmax(4, X[,10] + pmax(1, X[,11])))
    e = 1/(1 + exp(X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]))
    tau = rep(10, n) + X[,12]
    
  } else if (setup == "D") {
    
    X = matrix(rnorm(n*p), n, p)
    b = (pmax(X[,1] + X[,2] + X[,3], 0) + pmax(X[,4] + X[,5], 0)) / 2 + X[,4] + 0.5 * X[,5] + 2 * X[,6] + 4*(X[,7] +  X[,8] +  X[,9] +  X[,10] + X[,11] +  X[,12])
    e = pmin(sin(pi*X[,6]), 1/(1 + exp(-X[,1]) + exp(-X[,2]) + exp(-X[,3]) + exp(-X[,4]) + exp(-X[,5])))
    tau = pmax(X[,1] + X[,2], 0) - pmax(X[,3] + X[,4] + X[,5], 0)
    
  } else {
    print("error in setup")
    break
  }
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, tau = tau, e = e)
  simulation$mu_0 <- simulation$b - 0.5*simulation$tau
  simulation$mu_1 <- simulation$b + 0.5*simulation$tau
  simulation$Y_0 <- simulation$mu_0 + rnorm(n, mean = 0, sd = 0.1)
  simulation$Y_1 <- simulation$mu_1 + rnorm(n, mean = 0, sd = 0.1)
  simulation$A <- rbinom(n, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }
  return(simulation)
}

generate_simulation_linear <- function(n_obs = 1000, independent_covariate = TRUE, constant_cate = FALSE, all_covariates_output = FALSE){
  
  p = 50
  
  # generate multivariate gaussian vector
  if(independent_covariate){
    cov_mat = diag(p)
  } else {
    cov_mat = toeplitz(0.6^(0:(p - 1)))
  }
  
  
  X = rmvnorm(n = n_obs, mean = rep(1, p), sigma = cov_mat)
  
  # generate baseline and propensity scores
  b = X[,1:30]%*%rep(1,30)
  e = 1/(1 + exp(-2 - 0.8*(-X[,1] - X[,2] - X[,3] - X[,4] - X[,5] - X[,6])))
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, e = e)
  simulation$mu_0 <- simulation$b
  simulation$Y_0 <- simulation$mu_0 + rnorm(n_obs)
  
  if(constant_cate){
    ATE = 3
    simulation$mu_1 <- simulation$b + ATE
  } else {
    simulation$mu_1 <- simulation$b + simulation$X.2 + simulation$X.3 + simulation$X.4
  }
  
  simulation$Y_1 <- simulation$mu_1 + 2*rnorm(n_obs)
  
  simulation$A <- rbinom(n_obs, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }
}
