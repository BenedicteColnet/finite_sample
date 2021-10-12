# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(grf)
library(glmnet)
library(ranger) # efficient forest
library(splines) # function bs() for splines
library(SuperLearner)
library(tmle) # tmle


# simulation set up according to Wager & Nie
generate_simulation <- function(n = 1000, p = 12, setup = "D"){
  
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
    tau = rep(1, n)
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

# Estimators

## Custom AIPW with forest
aipw_forest <- function(covariates_names_vector_treatment, 
                        covariates_names_vector_outcome,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        n.folds = 2,
                        min.node.size.if.forest = 5) {
  
  n <- nrow(dataframe)
  indices <- split(seq(n), sort(seq(n) %% n.folds))
  
  t0 = rep(0, n)
  t1 = rep(1, n)
  
  X_t <- dataframe[, covariates_names_vector_treatment]
  X_o <- dataframe[, covariates_names_vector_outcome]
  W <- dataframe[, treatment_name]
  Y <- dataframe[, outcome_name]
  
  xt <- cbind(X_o, W)
  xt0 <- cbind(X_o, t0)
  xt1 <- cbind(X_o, t1)
  
  mu.hat.1 <- rep(NA, n)
  mu.hat.0 <- rep(NA, n)
  e.hat <- rep(NA, n)
  
  # cross-fitting of nuisance parameters
  for (idx in indices) {
    # Estimation
    outcome.model.treated <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 1, covariates_names_vector_treatment], dataframe[-idx & dataframe[,treatment_name] == 1, outcome_name], num.trees = 100, min.node.size = min.node.size.if.forest)
    outcome.model.control <- regression_forest(dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_treatment], dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name], num.trees = 100, min.node.size = min.node.size.if.forest)
    
    #propensity.model <- probability_forest(dataframe[-idx, covariates_names_vector_outcome], as.factor(dataframe[,treatment_name][-idx]), num.trees=100, min.node.size=min.node.size.if.forest)
    propensity.model <- cv.glmnet(x=dataframe[-idx, covariates_names_vector_outcome], y=W[-idx], family="binomial")
    
    # Prediction
    mu.hat.1[idx] <- predict(outcome.model.treated, data = xt1[idx,])$predictions
    mu.hat.0[idx] <- predict(outcome.model.control, data = xt0[idx,])$predictions
    e.hat[idx] <- predict(propensity.model, data = X_t[idx,])$predictions[,2]
  }
  
  # compute estimates
  aipw = mean(mu.hat.1 - mu.hat.0
              + W / e.hat * (Y -  mu.hat.1)
              - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  
  g_formula = mean(mu.hat.1) - mean(mu.hat.0)
  
  res = c("aipw" = aipw,
          "ipw" = ipw,
          "t-learner" = g_formula)
  
  return(res)
}

## Wrapper for causal forest

causal_forest_wrapper <- function(covariates_names_vector, 
                                  dataframe,
                                  outcome_name = "Y",
                                  treatment_name = "A"){
  
  forest <- causal_forest(
    X=dataframe[, covariates_names_vector],  
    W=dataframe[, treatment_name],
    Y=dataframe[, outcome_name],
    num.trees = 100)
  forest.ate <- average_treatment_effect(forest)
  return(forest.ate[[1]])
}

## Wrapper for TMLE

tmle_wrapper <- function(covariates_names_vector, 
                         dataframe,
                         outcome_name = "Y",
                         treatment_name = "A"){
  
  SL.library<- c("SL.glm", "SL.glmnet")
  # ranger takes way more time
  #SL.library<- c("SL.glm", "SL.glm.interaction", "SL.glmnet", "SL.ranger")
  
  TMLE <- tmle(Y = dataframe[,outcome_name],
               A = dataframe[,treatment_name],
               W = dataframe[,covariates_names_vector],
               family = "gaussian",
               Q.SL.library = SL.library,
               g.SL.library = SL.library)
  
  return(TMLE$estimates$ATE$psi)
}



###########

# Simulation pipeline

results <- data.frame("sample.size" = c(),
                      "estimate" = c(),
                      "estimator" = c(),
                      "subset" = c(),
                      "simulation" = c())

different_subset_tested <- c("all.covariates",
                             "smart",
                             "minimal.set")

for (sample.size in c(100, 300, 1000, 3000)){
  print(paste0("sample size: ", sample.size))
  for (i in 1:30){
    print(paste0("repetitions: ", i))
    for (simulation_setup in c("A", "B", "C", "D")){
      
      # generate a simulation
      a_simulation <- generate_simulation(n= sample.size, setup = simulation_setup)
      
      
      # compute estimator
      for (method in different_subset_tested){
        if (method == "all.covariates"){
          X_treatment <- paste0("X.", 1:12)
          X_outcome <- paste0("X.", 1:12)
          X_naive <- paste0("X.", 1:12)
        } else
          if (method == "smart"){
          X_treatment <- paste0("X.", 1:3)
          X_outcome <- paste0("X.", 2:6)
          X_naive <- paste0("X.", 2:6)
        } else if (method == "minimal.set"){
          X_treatment <- paste0("X.", 2:3)
          X_outcome <- paste0("X.", 2:3)
          X_naive <- paste0("X.", 2:3)
        }
        #causal_forest_estimate <- causal_forest_wrapper(X_naive, dataframe = a_simulation)
        custom_aipw_forest <- aipw_forest(X_treatment, X_outcome, dataframe = a_simulation)
        #tmle_estimate <- tmle_wrapper(X_naive, dataframe = a_simulation)
        
        new.row <- data.frame("sample.size" = rep(sample.size, 3),
                              "estimate" = c(custom_aipw_forest["ipw"], 
                                             custom_aipw_forest["t-learner"], 
                                             custom_aipw_forest["aipw"]),
                              "estimator" = c("ipw", 
                                              "t-learner",
                                              "aipw"),
                              "subset" = rep(method, 3),
                              "simulation" = rep(simulation_setup, 3))
        results <- rbind(results, new.row)
      }
    }
  }
}

results$sample.size <- as.factor(results$sample.size)


# Save the results

write.csv(x=results, file="./data/2021-10-12-aipw-forest-outcome-glmnet-propensity.csv")