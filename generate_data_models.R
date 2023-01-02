# simulation set up according to Wager & Nie
generate_simulation_wager_nie <- function(n = 1000, p = 6, setup = "D", all_covariates_output = FALSE){
  
  # set-ups
  if (setup == "A"){
    
    X = matrix(runif(n*p, min=0, max=1), n, p)
    b = sin(pi * X[,1] * X[,2]) + 2 * (X[,3] - 0.5)^2 + X[,4] + 0.5 * X[,5]
    eta = 0.1
    e = pmax(eta, pmin(sin(pi * X[,1] * X[,2]), 1-eta))
    tau = (X[,1] + X[,2]) / 2
  
  } else if (setup == "B"){
    
    X = matrix(rnorm(n * p), n, p)
    b = pmax(0, X[,1] + X[,2], X[,3]) + pmax(0, X[,4] + X[,5])
    e = 1/(1 + exp( 1 -X[,1] - X[,2]))  # this part is changed
    tau = X[,1] + log(1 + exp(X[,2]))
    
  } else if (setup == "C") {
    
    X = matrix(rnorm(n * p), n, p)
    b = 2 * log(1 + exp(X[,1] + X[,2] + X[,3]))
    e = 1/(1 + exp(X[,2] + X[,3]))
    tau = rep(1, n)
    
  } else if (setup == "D") {
    
    X = matrix(rnorm(n*p), n, p)
    b = (pmax(X[,1] + X[,2] + X[,3], 0) + pmax(X[,4] + X[,5], 0)) / 2
    e = 1/(1 + exp(-X[,1]) + exp(-X[,2]))
    tau = pmax(X[,1] + X[,2] + X[,3], 0) - pmax(X[,4] + X[,5], 0)
    
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
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y", "Y_1", "mu_1", "Y_0", "mu_0", "e")]
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }

}

generate_simulation_linear <- function(n = 1000, independent_covariate = TRUE, constant_cate = FALSE, all_covariates_output = FALSE){
  
  p = 50
  
  # generate multivariate gaussian vector
  if(independent_covariate){
    cov_mat = diag(p)
  } else {
    cov_mat = toeplitz(0.6^(0:(p - 1)))
  }
  
  
  X = rmvnorm(n = n, mean = rep(1, p), sigma = cov_mat)
  
  # generate baseline and propensity scores
  b = X[,1:30]%*%rep(1,30)
  e = 1/(1 + exp(-2 - 0.8*(-X[,1] - X[,2] - X[,3] - X[,4] - X[,5] - X[,6])))
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, e = e)
  simulation$mu_0 <- simulation$b
  simulation$Y_0 <- simulation$mu_0 + rnorm(n)
  
  if(constant_cate){
    ATE = 3
    simulation$mu_1 <- simulation$b + ATE
  } else {
    simulation$mu_1 <- simulation$b + simulation$X.2 + simulation$X.3 + simulation$X.4
  }
  
  simulation$Y_1 <- simulation$mu_1 + 2*rnorm(n)
  
  simulation$A <- rbinom(n, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  
  if(all_covariates_output){
    return(simulation)
  } else {
    simulation <- simulation[, c(paste0("X.", 1:p), "A", "Y")]
    return(simulation)
  }
}


generate_simulation_leborgne_foucher <- function(n = 1000, setup = "simple", all_covariates_output = FALSE){
  
  
  if (setup == "simple"){
    
    # Generate covariates - independent
    
    bfort <- log(3) 
    bmodere <- log(1.5)
    btreat <- log(1.75)
    
    .x1 <- rbinom(n, 1, prob=0.5)
    .x2 <- rbinom(n, 1, prob=0.5)
    .x3 <- rnorm(n, 0, 1)
    
    .x4 <- rbinom(n, 1, prob=0.5)
    .x5 <- rbinom(n, 1, prob=0.5)
    .x6 <- rnorm(n, 0, 1)
    
    .x7 <- rbinom(n, 1, prob=0.5)
    .x8 <- rbinom(n, 1, prob=0.5)
    .x9 <- rnorm(n, 0, 1)
    
    
    # Generate propensity score
    bx <- -0.8 + bfort*.x1 + bmodere*.x2 - bfort*.x4 - bmodere*.x5 + bfort*.x7 + bmodere*.x8
    e <- (exp(bx) / (1 + exp(bx)))
    rm(bx)
    
    A <- rbinom(n, 1, prob = e)
    
    
    # Outcome
    bx_1 <- -0.8 + btreat + bfort*.x1 - bfort*.x2 - bfort*.x3 + bmodere*.x4 - bmodere*.x5 + bmodere*.x6
    bx_0 <- -0.8 + bfort*.x1 - bfort*.x2 - bfort*.x3 + bmodere*.x4 - bmodere*.x5 + bmodere*.x6
    
    pr.Y_1 <- (exp(bx_1) / (1 + exp(bx_1)))
    pr.Y_0 <- (exp(bx_0) / (1 + exp(bx_0)))
    
    Y_1 <- rbinom(n, 1, prob = pr.Y_1) 
    Y_0 <- rbinom(n, 1, prob = pr.Y_0) 
    
    # Store covariates in dataframe
    simulation <- data.frame(X.1=.x1, X.2 = .x2, X.3 = .x3, X.4 = .x4, X.5 = .x5, X.6 = .x6, X.7 = .x7, X.8 = .x8, X.9 = .x9)
    
    covariate.set = paste0("X.", 1:9)
    
  } else if (setup == "complex") {
    
    # Generate covariates
    
    b0.t <- (-0.4)
    b0.o <- (-1.1)
    b1.o <- log(2)
    b1.t <- log(2)
    b1.c <- log(2)
    
    .x1 <- rnorm(n, 0, 1)
    .x2 <- rnorm(n, b0.t + b1.t * .x1, 1)
    .x3 <- rnorm(n, b0.t - b1.t * .x1 - b1.t * .x2, 1)
    .x4 <- rnorm(n, b0.t + b1.t * .x3, 1)
    .x5 <- rnorm(n, 0, 1)
    .x6 <- 1 * (rnorm(n, 0, 1) > 0.66) # prevalence ~25%
    .x7 <- 1 * (rnorm(n, b0.t - b1.t * .x5, 1) > (-0.40)) # prevalence ~50%
    .x8 <- rnorm(n, b0.t - b1.t * .x6, 1)
    .x9 <- 1 * (rnorm(n, b0.t + b1.t * .x7, 1) > (-0.80)) # prevalence ~75%
    .x10 <- rnorm(n, b0.t + b1.t * .x8, 1)
    .x11 <- rnorm(n, 0, 1)
    .x12 <- 1 * (rnorm(n, b0.t + b1.t * .x9, 1) > (0.84)) # prevalence ~25%
    .x13 <- 1 * (rnorm(n, b0.t + b1.t * .x10, 1) > (-0.09)) # prevalence ~50%
    .x14 <- rnorm(n, b0.t - b1.t * .x12 - b1.t * .x11, 1)
    .x15 <- rnorm(n, b0.t - b1.t * .x12, 1)
    .x16 <- 1 * (rnorm(n, 0, 1) > (-0.66)) # prevalence ~75%
    .x17 <- 1 * (rnorm(n, b0.t - b1.t * .x16, 1) > (-0.92)) # prevalence ~50%
    .x18 <- rnorm(n, 0, 1)
    .x19 <- 1 * (rnorm(n, 0, 1) > (0.66))  # prevalence ~25%
    .x20 <- 1 * (rnorm(n, 0, 1) > (0.66))  # prevalence ~25%
    .x21 <- rnorm(n, 0, 1)
    .x22 <- 1 * (rnorm(n, 0, 1) > (0.66)) # prevalence ~25%
    
    
    
    # Generate propensity score
    bx <- b0.t + b1.t * .x1 -
      b1.t * .x3 +
      b1.t * .x5 -
      b1.t * .x7  +
      b1.t * .x9  -
      b1.t * .x11 + 
      b1.t * .x13 -
      b1.t * .x15 -
      b1.t * .x17 +
      b1.t * .x19 -
      b1.t * .x21
    e <- (exp(bx) / (1 + exp(bx)))
    rm(bx)
    
    A <- rbinom(n, 1, prob = e)
    
    
    # Outcome
    
    # median(data.obs$x2) # -0.40
    med.x2 <- -0.40
    # median(data.obs$x15) # -0.57
    med.x15 <- -0.57
    
    
    bx_1 <- b0.o + b1.o * (.x2 > med.x2) -
      b1.o * .x3 + (b1.o / 2) * (.x3^2) + b1.o * .x6 +
      b1.o * .x7 + b1.o * .x10 + b1.o * 0.5 * (.x11^2) -
      b1.o * .x14 - b1.o * (.x15 > med.x15) +  b1.o * .x18 +
      b1.o * .x19 + b1.c  + b1.o * 0.5 * .x18
     
    
    bx_0 <- b0.o + b1.o * (.x2 > med.x2) -
      b1.o * .x3 + (b1.o / 2) * (.x3^2) + b1.o * .x6 +
      b1.o * .x7 + b1.o * .x10 + b1.o * 0.5 * (.x11^2) -
      b1.o * .x14 - b1.o * (.x15 > med.x15) +  b1.o * .x18 +
      b1.o * .x19 + b1.c * 0 + b1.o * 0.5 * 0 * .x18
    
    pr.Y_1 <- (exp(bx_1) / (1 + exp(bx_1)))
    pr.Y_0 <- (exp(bx_0) / (1 + exp(bx_0)))
    
    Y_1 <- rbinom(n, 1, prob = pr.Y_1) 
    Y_0 <- rbinom(n, 1, prob = pr.Y_0) 
    
    # Store covariates in dataframe
    simulation <- data.frame(X.1 = .x1, X.2 = .x2, X.3 = .x3, X.4 = .x4, X.5 = .x5, X.6 = .x6, X.7 = .x7, X.8 = .x8, X.9 = .x9,
                             X.10 = .x10, X.11 = .x11, X.12 = .x12, X.13 = .x13, X.14 = .x14, X.15 = .x15, X.16 = .x16, X.17 = .x17,
                             X.18 = .x18, X.19 = .x1, X.20 = .x20, X.21 = .x21, X.22 = .x22)
    
    
    covariate.set = paste0("X.", 1:22)
    
  } else{
    
    print("Error in setup")
    break
    
  }
    
  simulation$e <- e
  simulation$A <- A
  simulation$Y_1 <- Y_1
  simulation$Y_0 <- Y_0
  
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  

  # return output
  if(all_covariates_output){
    simulation <- simulation[, c(covariate.set, "A", "Y", "Y_1", "Y_0", "e")]
    return(simulation)
  } else {
    simulation <- simulation[, c(covariate.set, "A", "Y")]
    return(simulation)
  }
}
