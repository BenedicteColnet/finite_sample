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

results.linear <- data.frame("sample.size" = c(),
                             "estimate" = c(),
                             "estimator" = c(),
                             "subset" = c(),
                             "simulation" = c(),
                             "cross-fitting" = c(),
                             "independence" = c(),
                             "nuisance" = c())

different_subset_tested <- c("extended",
                             "smart",
                             "minimal")

for (sample.size in c(50, 100, 200, 300, 500, 1000, 3000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    
    # generate a simulation
    a_simulation_for_mu <- generate_simulation_DML(n = sample.size)
    a_simulation_for_e <- generate_simulation_DML(n = sample.size)
    a_simulation_for_estimate <- generate_simulation_DML(n = sample.size)
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", 1:30)
        X_outcome <- paste0("X.", 1:30)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 1:3)
        X_outcome <- paste0("X.", 1:30)
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", 1:3)
        X_outcome <- paste0("X.", 1:3)
      } else {
        stop("error in subset.")
      }
      
      for (number_of_folds in c(5)){
        
        # fit models
        outcome.model.treated <-  ranger(Y ~ .,  
                                         num.trees = 5000, 
                                         max.depth = 0,
                                         min.node.size = 1, 
                                         data = a_simulation_for_mu[a_simulation_for_mu$A == 1, c("Y", X_outcome)])
        outcome.model.control <-  ranger(Y ~ .,  
                                         num.trees = 5000, 
                                         max.depth = 0,
                                         min.node.size = 1, 
                                         data = a_simulation_for_mu[a_simulation_for_mu$A == 0, c("Y", X_outcome)])
        
        propensity.model.on.same.fold <- probability_forest(a_simulation_for_mu[, X_treatment], 
                                               as.factor(a_simulation_for_mu[, "A"]), 
                                               num.trees = 5000, 
                                               min.node.size=1)
        
        propensity.model.on.other.fold <- probability_forest(a_simulation_for_e[, X_treatment], 
                                                            as.factor(a_simulation_for_e[, "A"]), 
                                                            num.trees = 5000, 
                                                            min.node.size=1)
        
        
        # prediction and estimation
        mu.hat.1 <- predict(outcome.model.treated, a_simulation_for_estimate[, X_outcome])$predictions
        mu.hat.0 <- predict(outcome.model.control, a_simulation_for_estimate[, X_outcome])$predictions
        e.hat.same.fold <- predict(propensity.model.on.same.fold, 
                                   newdata = a_simulation_for_estimate[,X_treatment])$predictions[,2]
        e.hat.other.fold  <- predict(propensity.model.on.other.fold, 
                                   newdata = a_simulation_for_estimate[,X_treatment])$predictions[,2]
        
        W <- a_simulation_for_estimate[, "A"]
        Y <- a_simulation_for_estimate[, "Y"]
        
        aipw.same.fold = mean(mu.hat.1 - mu.hat.0
                    + W / e.hat.same.fold * (Y -  mu.hat.1)
                    - (1 - W) / (1 - e.hat.same.fold) * (Y -  mu.hat.0))
        aipw.other.fold = mean(mu.hat.1 - mu.hat.0
                              + W / e.hat.other.fold * (Y -  mu.hat.1)
                              - (1 - W) / (1 - e.hat.other.fold) * (Y -  mu.hat.0))
        
        
        new.row <- data.frame("sample.size" = rep(sample.size, 2),
                              "estimate" = c(aipw.same.fold, aipw.other.fold),
                              "estimator" = rep(c("aipw"),2),
                              "subset" = rep(method, 2),
                              "simulation" = c("same.fold", "other.fold") ,
                              "cross-fitting" = rep(NA, 2),
                              "independence" = rep(NA, 2),
                              "nuisance" = rep("forest",2))
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}

write.csv(x=results.linear, file="./data/new.csv")

