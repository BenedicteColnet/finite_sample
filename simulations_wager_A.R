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

for (sample.size in c(100, 200, 400, 5000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    
    # generate a simulation
    a_simulation <- generate_simulation_wager_nie(n = sample.size, setup = "A")
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", 1:6)
        X_outcome <- paste0("X.", 1:6)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 1:2)
        X_outcome <- paste0("X.", 1:6)
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", 1:2)
        X_outcome <- paste0("X.", 1:2)
      } else {
        stop("error in subset.")
      }
      
      for (number_of_folds in c(2)){

        custom_aipw <- aipw_forest(X_treatment, 
                                     X_outcome, 
                                     dataframe = a_simulation,
                                     min.node.size.if.forest = 1,
                                     n.folds = number_of_folds)
        
        new.row <- data.frame("sample.size" = rep(sample.size, 3),
                              "estimate" = c(custom_aipw["ipw"],
                                             custom_aipw["t.learner"],
                                             custom_aipw["aipw"]),
                              "estimator" = rep(c("ipw",
                                                  "t-learner",
                                                  "aipw"),1),
                              "subset" = rep(method, 3),
                              "simulation" = rep("A", 3),
                              "cross-fitting" = rep(number_of_folds, 3),
                              "independence" = rep(NA,3),
                              "nuisance" = rep("forest",3))
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}

write.csv(x=results.linear, file="./data/A.csv")
