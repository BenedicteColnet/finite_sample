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
                             "cross-fitting" = c())

different_subset_tested <- c("all.covariates",
                             "all.covariates.wo.instruments",
                             "smart",
                             "minimal.set")

for (sample.size in c(3000, 9000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:20){
    print(paste0("Repetition:", i))
    # generate a simulation
    a_simulation <- generate_simulation_wager_nie(n = sample.size, setup = "C")
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "all.covariates"){
        X_treatment <- paste0("X.", 1:12)
        X_outcome <- paste0("X.", 1:12)
      } else if (method == "all.covariates.wo.instruments"){
        X_treatment <- paste0("X.", 2:12)
        X_outcome <- paste0("X.", 2:12)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 2:6)
        X_outcome <- paste0("X.", 2:3)
      } else if (method == "minimal.set"){
        X_treatment <- paste0("X.", 2:3)
        X_outcome <- paste0("X.", 2:3)
      } else {
        stop("error in subset.")
      }
      
      for (number_of_folds in c(20)){
        
        SL.o = c("SL.mean", "SL.lm", "SL.ranger", "SL.glmnet")
        SL.t = c("SL.glm", "SL.mean", "SL.ranger", "SL.glmnet")
        
        
        custom_aipw <- aipw_ML(covariates_names_vector_treatment = X_treatment, 
                               covariates_names_vector_outcome = X_outcome, 
                               dataframe = a_simulation, 
                               n.folds = number_of_folds,
                               sl_libs_outcome = SL.o,
                               sl_libs_treatment = SL.t)
        
        aipw.wrapper <- aipw_wrapped(covariates_names_vector_treatment = X_treatment, 
                                     covariates_names_vector_outcome = X_outcome, 
                                     dataframe = a_simulation, 
                                     n.folds = number_of_folds,
                                     sl_libs_outcome = SL.o,
                                     sl_libs_treatment = SL.t)
        
        tmle.wrapper <- tmle_wrapped(covariates_names_vector = X_outcome, 
                                     dataframe = a_simulation, 
                                     n.folds = number_of_folds,
                                     sl_libs_outcome = SL.o,
                                     sl_libs_treatment = SL.t)
        
        grf.wrapper <- causal_forest_wrapper(covariates_names_vector = X_outcome, 
                                             dataframe = a_simulation)
        
        
        
        new.row <- data.frame("sample.size" = rep(sample.size, 6),
                              "estimate" = c(custom_aipw["ipw"],
                                             custom_aipw["t.learner"],
                                             custom_aipw["aipw"],
                                             aipw.wrapper,
                                             tmle.wrapper,
                                             grf.wrapper),
                              "estimator" = c("ipw",
                                              "t-learner",
                                              "aipw",
                                              "wrapper aipw",
                                              "wrapper tmle",
                                              "wrapper grf"),
                              "subset" = rep(method, 6),
                              "simulation" = rep("wager-C", 6),
                              "cross-fitting" = c(rep(number_of_folds, 4), NA, NA))
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}


write.csv(x=results.linear, file="./data/2021-10-26-wager-C.csv")