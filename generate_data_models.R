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
    b = sin(pi * X[,2] * X[,3]) + 2 * (X[,4] - 0.5)^2 + X[,5]^2 + 3 * X[,6] + 3*X[,4]
    eta = 0.1
    e = pmax(eta, pmin(sin(pi * X[,2] * X[,3]), 1-eta))
    tau = (X[,3]*X[,4]) / 2
  } else if (setup == "B"){
    X = matrix(rnorm(n * p), n, p)
    b = pmax(0, X[,2] + X[,3], X[,4]) + pmax(0, X[,5] + X[,6])
    e = 0.5
    tau = X[,2] + log(1 + exp(X[,3]))
  } else if (setup == "C") {
    X = matrix(rnorm(n * p), n, p)
    b = 2 * log(1 + exp(X[,1] + X[,2] + X[,3])) + X[,4]^2 + log(1 + abs(X[,5])) + 3*X[,6]
    e = 1/(1 + exp(X[,1] + X[,2]))
    tau = rep(1, n)
  } else if (setup == "D") {
    X = matrix(rnorm(n * p), n, p)
    b = (pmax(X[,1] + X[,2] + X[,3], 0) + pmax(X[,4] + X[,5], 0)) / 2
    e = 1/(1 + exp(-X[,1]) + exp(-X[,2]))
    tau = pmax(X[,1] + X[,2] + X[,3], 0) - pmax(X[,4] + X[,5], 0)
  } else if (setup == "D.bis"){
    X = matrix(rnorm(n * p), n, p)
    b = 10*X[,1] + (pmax(X[,1] + X[,2] + X[,3], 0) + 3*pmax(X[,4] + X[,5], 0)) + 3*X[,6]*X[,6] + 0.4*log(1+exp(X[7,]))
    e = 1/(1 + exp(-X[,1]) + exp(-X[,2]))
    tau = -9*X[,1] + X[,2] - pmax(X[,4] + X[,5], 0)
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
