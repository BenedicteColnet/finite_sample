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

for (sample.size in c(100, 300, 1000, 3000, 10000, 30000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    
    # generate a simulation
    a_simulation <- generate_simulation_wager_nie(n = sample.size, setup = "D", all_covariates_output = TRUE)
    
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
                                   n.folds = number_of_folds,
                                   return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = sample.size,
                              "estimate" = custom_aipw["aipw"],
                              "estimator" = "aipw",
                              "subset" = method,
                              "simulation" = "D",
                              "cross-fitting" = number_of_folds,
                              "nuisance" = "forest",
                              "term.A" = custom_aipw["term.A"], 
                              "term.B" = custom_aipw["term.B"], 
                              "term.C" = custom_aipw["term.C"],
                              "term.D" = custom_aipw["term.D"], 
                              "term.E" = custom_aipw["term.E"], 
                              "term.F" = custom_aipw["term.F"])
        results.linear <- rbind(results.linear, new.row)
        
      }
      
      
      
      if (method == "extended"){
        custom_ipw <- ipw_forest(covariates_names = paste0("X.", 1:6), 
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 1,
                                 n.folds = number_of_folds,
                                 return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = sample.size,
                              "estimate" = custom_ipw,
                              "estimator" = "ipw",
                              "subset" = method,
                              "simulation" = "D",
                              "cross-fitting" = number_of_folds,
                              "nuisance" = "forest",
                              "term.A" = NA, 
                              "term.B" = NA, 
                              "term.C" = NA,
                              "term.D" = NA, 
                              "term.E" = NA, 
                              "term.F" = NA)
        results.linear <- rbind(results.linear, new.row)
        
        
        custom_tl <- t_learner_forest(covariates_names = paste0("X.", 1:6), 
                                      dataframe = a_simulation,
                                      min.node.size.if.forest = 1,
                                      n.folds = number_of_folds,
                                      return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = sample.size,
                              "estimate" = custom_tl,
                              "estimator" = "t-learner",
                              "subset" = method,
                              "simulation" = "D",
                              "cross-fitting" = number_of_folds,
                              "nuisance" = "forest",
                              "term.A" = NA, 
                              "term.B" = NA, 
                              "term.C" = NA,
                              "term.D" = NA, 
                              "term.E" = NA, 
                              "term.F" = NA)
        results.linear <- rbind(results.linear, new.row)
      } else if (method == "minimal"){
        
        custom_ipw <- ipw_forest(covariates_names = paste0("X.", 1:2), 
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 1,
                                 n.folds = number_of_folds,
                                 return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = sample.size,
                              "estimate" = custom_ipw,
                              "estimator" = "ipw",
                              "subset" = method,
                              "simulation" = "D",
                              "cross-fitting" = number_of_folds,
                              "nuisance" = "forest",
                              "term.A" = NA, 
                              "term.B" = NA, 
                              "term.C" = NA,
                              "term.D" = NA, 
                              "term.E" = NA, 
                              "term.F" = NA)
        results.linear <- rbind(results.linear, new.row)
        
        
        custom_tl <- t_learner_forest(covariates_names = paste0("X.", 1:2), 
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 1,
                                 n.folds = number_of_folds,
                                 return.decomposition = TRUE)
        
        new.row <- data.frame("sample.size" = sample.size,
                              "estimate" = custom_tl,
                              "estimator" = "t-learner",
                              "subset" = method,
                              "simulation" = "D",
                              "cross-fitting" = number_of_folds,
                              "nuisance" = "forest",
                              "term.A" = NA, 
                              "term.B" = NA, 
                              "term.C" = NA,
                              "term.D" = NA, 
                              "term.E" = NA, 
                              "term.F" = NA)
        results.linear <- rbind(results.linear, new.row)
        
      } else {
        stop("error in subset.")
      }
      
      
    }
  }
}

write.csv(x=results.linear, file="./data/Dn.csv")
