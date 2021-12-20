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

for (sample.size in c(300, 1000, 2000, 5000, 10000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    
    # generate a simulation
    a_simulation <- generate_simulation_wager_nie(n = sample.size, setup = "D")
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", 1:5)
        X_outcome <- paste0("X.", 1:5)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 1:5)
        X_outcome <- paste0("X.", 1:2)
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", 1:2)
        X_outcome <- paste0("X.", 1:2)
      } else {
        stop("error in subset.")
      }
      
      for (number_of_folds in c(5)){
        
        aipw <- aipw_forest(X_treatment, X_outcome, dataframe = a_simulation,
                            n.folds = 2,
                            min.node.size.if.forest = 1)
        
        aipw_two_folds <- aipw_forest_two_fold(X_treatment, X_outcome, dataframe = a_simulation,
                            min.node.size.if.forest = 1)
        
        
        new.row <- data.frame("sample.size" = rep(sample.size, 6),
                              "estimate" = c(aipw["ipw"],
                                             aipw["t.learner"],
                                             aipw["aipw"],
                                             aipw_two_folds["ipw"],
                                             aipw_two_folds["t.learner"],
                                             aipw_two_folds["aipw"]),
                              "estimator" = rep(c("ipw",
                                                  "t-learner",
                                                  "aipw"),2),
                              "subset" = rep(method, 2),
                              "simulation" = c(rep("cross.fitting", 3), rep("drop.fold", 3)) ,
                              "cross-fitting" = rep(2, 6),
                              "independence" = rep(NA,6),
                              "nuisance" = rep("forest",6))
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}

write.csv(x=results.linear, file="./data/new.D.csv")