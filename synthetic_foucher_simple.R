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
                             "nuisance" = c())

different_subset_tested <- c("extended",
                             "minimal")

for (sample.size in c(300, 1000, 3000, 10000, 30000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:100){
    
    # generate a simulation
    a_simulation <- generate_simulation_leborgne_foucher(n = sample.size, setup = "simple", all_covariates_output = FALSE)
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", c(3, 11, 15, 7, 19, 2, 10, 6, 18, 14))
        X_outcome <- paste0("X.", c(3, 11, 15, 7, 19, 2, 10, 6, 18, 14))
        
      } else if (method == "smart"){
        X_treatment <- paste0("X.", c(3, 11, 15, 7, 19))
        X_outcome <- paste0("X.", c(3, 11, 15, 7, 19, 2, 10, 6, 18, 14))
        
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", c(3, 11, 15, 7, 19))
        X_outcome <-paste0("X.", c(3, 11, 15, 7, 19))
        
        
      } else {
        stop("error in subset.")
      }
      
      
      custom_aipw <- aipw_forest(X_treatment,
                                 X_outcome,
                                 dataframe = a_simulation,
                                 min.node.size.if.forest = 5,
                                 n.folds = 2,
                                 return.decomposition = FALSE)


      custom_aipw_linear <- aipw_linear(X_treatment,
                  X_outcome,
                  dataframe = a_simulation,
                  n.folds = 2)
      
      
      wrapper_for_grf <- causal_forest_wrapper(X_outcome, dataframe = a_simulation)
      
      wrapper_for_tmle <- tmle_wrapper(covariates_names_vector = X_outcome, dataframe = a_simulation)
      
      
      estimates <- c(custom_aipw["ipw"], 
                     custom_aipw["aipw"], 
                     custom_aipw["t.learner"], 
                     custom_aipw_linear["ipw"],
                     custom_aipw_linear["t.learner"],
                     custom_aipw_linear["aipw"],
                     wrapper_for_grf,
                     wrapper_for_tmle)
      
      
      new.row <- data.frame("sample.size" = rep(sample.size, 8),
                            "estimate" = estimates,
                            "estimator" = c("ipw", "aipw", "g-formula", "ipw", "g-formula", "aipw", 'cf', 'tmle'),
                            "subset" = rep(method, 8),
                            "nuisance" = c('forest', 'forest', 'forest', "linear", "linear", "linear", "forest", "forest"))
      
                            
      results.linear <- rbind(results.linear, new.row)
    }
  }
}

write.csv(x=results.linear, file="./data/synthetic_foucher_simple.csv")
