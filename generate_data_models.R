# simulation set up according to Wager & Nie
generate_simulation_wager_nie <- function(n = 1000, p = 12, setup = "D", all_covariates_output = FALSE){
  
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
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }
  return(simulation)
}




generate_simulation_linear <- function(n_obs = 1000, independent_covariate = FALSE, constant_cate = TRUE, all_covariates_output = FALSE){
  
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
  e = 1/(1 + exp(-4 - 0.8*(-X[,1] - X[,2] - X[,3] - X[,4] - X[,5] - X[,6])))
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, e = e)
  simulation$Y_0 <- simulation$b + rnorm(n_obs)
  
  if(constant_cate){
    ATE = 3
    simulation$Y_1 <- simulation$b + ATE + 2*rnorm(n_obs)
  } else {
    simulation$Y_1 <- simulation$b + simulation$X.2 + simulation$X.3 + simulation$X.4 + 2*rnorm(n_obs)
  }
  
  simulation$A <- rbinom(n_obs, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }
  
}

# original simulation set-up
complex_model <- function(n_obs = 1000,  all_covariates_output = FALSE){
  
  # generate multivariate gaussian vector for 5 covariates
  cov_mat = toeplitz(0.6^(0:(5 - 1)))
  X = rmvnorm(n = n_obs, mean = rep(1, 5), sigma = cov_mat)
  
  # generate another vector with 5 uniform covariates
  X_bis = matrix(runif(n_obs*5, min=0, max=1), n_obs, 5)
  
  # a binomal covariates
  X.1 <- rbinom(n_obs, 1, prob = 0.2)
  
  # a covariates depending on one of the uniform
  X.2 <- rbinom(n_obs, 1, prob = (0.75 * X_bis[,1] + (0.25 * (1 - X_bis[,1]))))
  
  
  X <- data.frame( X.1, X.2, X, X_bis)
  names(X) <- paste("X.", 1:12)
  
  eta = 0.1
  e = pmax(eta,as.numeric(X[,1]==0 & X[,2]==0 & X[,3]<0.6)*0.8, pmin(sin(pi * X[,1] * X[,2]), 1-eta), 1/(1 + exp(-X[,3]) + exp(-X[,4] + exp(-X[,5]))))
  
  b = (pmax(X[,2] + X[,3] + X[,4], 0) + pmax(X[,5] + X[,6], 0)) / 2 + X[,10]*X[,11] + X[,2]*4
  tau = (X[,1] + X[,7] + X[,8] + X[,9])*3
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, tau = tau, e = e)
  simulation$Y_0 <- simulation$b - 0.5*simulation$tau + rnorm(n_obs, mean = 0, sd = 0.1)
  simulation$Y_1 <- simulation$b + 0.5*simulation$tau + rnorm(n_obs, mean = 0, sd = 0.1)
  simulation$A <- rbinom(n_obs, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:12), "A", "Y")]
    return(simulation)
  }
  return(simulation)
}


generate_simulation_logit_binary <- function(n_obs = 1000, independent_covariate = FALSE, all_covariates_output = FALSE){
  
  p = 12
  
  # generate multivariate gaussian vector
  if(independent_covariate){
    cov_mat = diag(p)
  } else {
    cov_mat = toeplitz(0.6^(0:(p - 1)))
  }
  
  
  X = rmvnorm(n = n_obs, mean = rep(1, p), sigma = cov_mat)
  
  # generate baseline and propensity scores
  b = X[,1:12]%*%rep(3,12)
  e = 1/(1 + exp(-4 - 0.8*(-X[,1] - X[,2] - X[,3] - X[,4])))
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, e = e)
  simulation$prob_Y_0 <- 1/(1 + exp(simulation$b))
  simulation$Y_0 <- rbinom(n_obs, size = 1, prob = simulation$prob_Y_0)
  simulation$prob_Y_1 <- 1/(1 + exp(simulation$b + 4))
  simulation$Y_1 <- rbinom(n_obs, size = 1, prob = simulation$prob_Y_1)
  
  simulation$A <- rbinom(n_obs, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }
  
}