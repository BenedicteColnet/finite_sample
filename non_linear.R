# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(grf)
library(glmnet)
library(ranger) # efficient forest
library(splines) # function bs() for splines
library(mvtnorm) # rmvnorm
library(tmle)
library(SuperLearner)
library(AIPW)

source("estimators.R")
source("generate_data_models.R")

generate_simulation <- function(n = 1000,  return_oracles = FALSE) {
  
  # same as lunceford
  tau_0 <- c(-1, -1, 1, 1)
  tau_1 <- tau_0 * -1
  Sigma_X3 <- matrix(c(1, 0.5, -0.5, -0.5, 
                       0.5, 1, -0.5, -0.5,
                       -0.5, -0.5, 1, 0.5, 
                       -0.5, -0.5, 0.5, 1), ncol = 4, byrow = TRUE)
  
  
  # beginning of simulations
  X.3 <- rbinom(n, 1, prob = 0.2)
  V.3 <- rbinom(n, 1, prob = (0.75 * X.3 + (0.25 * (1 - X.3))))
  hold <- rmvnorm(n,  mean = rep(0, 4), Sigma_X3)
  colnames(hold) <- c("X.1", "V.1", "X.2", "V.2")
  hold <- cbind(hold, X.3, V.3)
  hold <- apply(hold, 1, function(x){
    x[1:4] <- x[1:4] + tau_1^(x[5])*tau_0^(1 - x[5])
    x})
  X <- t(hold)[, c("X.1", "X.2", "X.3", "V.1", "V.2", "V.3")]
  
  # model for the propensity scores
  beta = c(0.6, -0.6, 0.6)
  e <- ifelse(X[,2] < 1, plogis(X[, 1:3] %*% beta), abs(sin(pi * 0.5 *X[,2])))
  mu_0 <- 5*X[,3]*(1 / (1 + exp(-X[,1]))) + 5* sin(X[,2]) + ifelse(X[,5] > 0 & X[,6], 9, -9) + 3*X[,4]*X[,4]
  mu_1 <- mu_0 - 5*X[,3]*(1 / (1 + exp(-X[,1]))) + 3*(1-X[,3])*(1 / (1 + exp(+X[,1] + X[,2]))) 
  
  simulation <- data.frame(X, e = e, mu_0 = mu_0, mu_1 = mu_1)
  simulation$A <- rbinom(n, size = 1, prob = simulation$e)
  simulation$Y_0 <- simulation$mu_0 + rnorm(n, mean = 0, sd = 0.1)
  simulation$Y_1 <- simulation$mu_1 + rnorm(n, mean = 0, sd = 0.1)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  return(simulation)
  
}


results.linear <- data.frame("sample.size" = c(),
                             "estimate" = c(),
                             "estimator" = c(),
                             "subset" = c(),
                             "simulation" = c(),
                             "cross-fitting" = c(),
                             "nuisance" = c(),
                             "term.A" = c(), 
                             "term.B" = c(), 
                             "term.C" = c(),
                             "term.D" = c(), 
                             "term.E" = c(), 
                             "term.F" = c())

different_subset_tested <- c("extended",
                             "smart",
                             "minimal")

for (sample.size in c(1000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:100){
    
    # generate a simulation
    a_simulation <- generate_simulation(n = sample.size, return_oracles = TRUE)
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- c("X.1", "X.2", "X.3", "V.1", "V.2", "V.3")
        X_outcome <- c("X.1", "X.2", "X.3", "V.1", "V.2", "V.3")
      } else if (method == "smart"){
        X_treatment <- c("X.1", "X.2", "X.3")
        X_outcome <- c("X.1", "X.2", "X.3", "V.1", "V.2", "V.3")
      } else if (method == "minimal"){
        X_treatment <- c("X.1", "X.2", "X.3")
        X_outcome <- c("X.1", "X.2", "X.3")
      } else {
        stop("error in subset.")
      }
      
      for (number_of_folds in c(2)){
        
        custom_aipw <- aipw_forest(X_treatment, 
                                   X_outcome, 
                                   dataframe = a_simulation,
                                   min.node.size.if.forest = 1,
                                   n.folds = number_of_folds,
                                   return.decomposition = TRUE)
        
        custom_ipw <- ipw_forest(covariates_names = X_treatment, 
                                dataframe = a_simulation,
                                min.node.size.if.forest = 1,
                                n.folds = number_of_folds,
                                return.decomposition = TRUE)
        
        
        custom_gformula <- t_learner_forest(covariates_names = X_outcome, 
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 1,
                                 n.folds = number_of_folds,
                                 return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = rep(sample.size, 3),
                              "estimate" = c(custom_aipw["aipw"],
                                             custom_gformula,
                                             custom_ipw),
                              "estimator" = rep(c("aipw",
                                                  "t-learner",
                                                  "ipw"),1),
                              "subset" = rep(method, 3),
                              "simulation" = rep("non-linear", 3),
                              "cross-fitting" = c(number_of_folds, NA, NA),
                              "nuisance" = rep("forest",3),
                              "term.A" = rep(custom_aipw["term.A"], 3), 
                              "term.B" = rep(custom_aipw["term.B"], 3), 
                              "term.C" = rep(custom_aipw["term.C"], 3),
                              "term.D" = rep(custom_aipw["term.D"], 3), 
                              "term.E" = rep(custom_aipw["term.E"], 3), 
                              "term.F" = rep(custom_aipw["term.F"], 3))
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}


write.csv(x=results.linear, file="./data/non_linear.csv")
