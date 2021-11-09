# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(glmnet)
library(grf)
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
                             "nuisance" = c())

different_subset_tested <- c("outcome.wo.instruments",
                             "minimal.set")

for (sample.size in c(300, 1000, 3000, 9000, 30000, 100000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    for (simulation in c("linear", "non.linear")){
      
      # generate a simulation
      if (simulation == "linear"){
        a_simulation <- generate_simulation_linear(n_obs = sample.size, independent_covariate = FALSE)
      } else {
        a_simulation <-  generate_simulation_wager_nie(n = sample.size, setup = "C")
      }
     
      
      # choose subset
      for (method in different_subset_tested){
        if (method == "outcome.and.instruments"){
          X_treatment <- paste0("X.", 1:10)
          X_outcome <- paste0("X.", 1:10)
        } else if (method == "outcome.wo.instruments"){
          X_treatment <- paste0("X.", 4:10)
          X_outcome <- paste0("X.", 4:10)
        } else if (method == "smart"){
          X_treatment <- paste0("X.", 4:7)
          X_outcome <- paste0("X.", 4:10)
        } else if (method == "minimal.set"){
          
        } else {
          stop("error in subset.")
        }
        
        custom_aipw_linear <- aipw_linear(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 2)
        custom_aipw_forest <- aipw_forest(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 2, min.node.size.if.forest = 1)
        
        
        new.row <- data.frame("sample.size" = rep(sample.size, 6),
                              "estimate" = c(custom_aipw_linear["ipw"],
                                             custom_aipw_linear["t.learner"],
                                             custom_aipw_linear["aipw"],
                                             custom_aipw_forest["ipw"],
                                             custom_aipw_forest["t.learner"],
                                             custom_aipw_forest["aipw"]),
                              "estimator" = rep(c("ipw",
                                                  "t-learner",
                                                  "aipw"),2),
                              "subset" = rep(method, 6),
                              "simulation" = rep(simulation, 6),
                              "cross-fitting" = rep(2, 6),
                              "nuisance" = c(rep("linear", 3)), rep("forest", 3))
        
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/2021-11-01-linear-constant-ate-linear-nuisance.csv")