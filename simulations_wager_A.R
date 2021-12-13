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

for (sample.size in c(1000, 1500, 2000, 3000, 4500, 9000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    
    # generate a simulation
    a_simulation <- generate_simulation_wager_nie(n = sample.size, setup = "A")
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", 2:6)
        X_outcome <- paste0("X.", 2:6)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 2:6)
        X_outcome <- paste0("X.", 2:3)
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", 2:3)
        X_outcome <- paste0("X.", 2:3)
      } else {
        stop("error in subset.")
      }
      
      for (number_of_folds in c(2)){
        
        #SL.o = c("SL.mean", "SL.lm", "SL.ranger", "SL.glmnet")
        #SL.t = c("SL.glm", "SL.mean", "SL.ranger", "SL.glmnet")
        
        
        # custom_aipw <- aipw_ML(covariates_names_vector_treatment = X_treatment, 
        #                        covariates_names_vector_outcome = X_outcome, 
        #                        dataframe = a_simulation, 
        #                        n.folds = number_of_folds,
        #                        sl_libs_outcome = SL.o,
        #                        sl_libs_treatment = SL.t)
        
        # aipw.wrapper <- aipw_wrapped(covariates_names_vector_treatment = X_treatment, 
        #                              covariates_names_vector_outcome = X_outcome, 
        #                              dataframe = a_simulation, 
        #                              n.folds = number_of_folds,
        #                              sl_libs_outcome = SL.o,
        #                              sl_libs_treatment = SL.t)
        
        # tmle.wrapper <- tmle_wrapper(covariates_names_vector = X_outcome, 
        #                              dataframe = a_simulation, 
        #                              n.folds = number_of_folds,
        #                              sl_libs_outcome = SL.o,
        #                              sl_libs_treatment = SL.t)
        
        # grf.wrapper <- causal_forest_wrapper(covariates_names_vector = X_outcome, 
        #                                      dataframe = a_simulation)
        
        
        custom_aipw_2 <- aipw_forest_two_fold(X_treatment, X_outcome, dataframe = a_simulation)
        custom_aipw_3 <- aipw_forest_three_fold(X_treatment, X_outcome, dataframe = a_simulation)
        
        
        new.row <- data.frame("sample.size" = rep(sample.size, 6),
                              "estimate" = c(custom_aipw_2["ipw"],
                                             custom_aipw_2["t.learner"],
                                             custom_aipw_2["aipw"],
                                             custom_aipw_3["ipw"],
                                             custom_aipw_3["t.learner"],
                                             custom_aipw_3["aipw"]),
                              "estimator" = rep(c("ipw",
                                                  "t-learner",
                                                  "aipw"),2),
                              "subset" = rep(method, 6),
                              "simulation" = rep("wager-A", 6),
                              "cross-fitting" = c("2 folds", "2 folds", "2 folds", "3 folds", "3 folds", "3 folds"),
                              "independence" = rep(NA,6),
                              "nuisance" = rep("forest",6))
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}


write.csv(x=results.linear, file="./data/2021-12-09-wager-A-two-folds.csv")
