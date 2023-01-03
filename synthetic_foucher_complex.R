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
library(DoubleML)
library(mlr3)
library(mlr3learners)

source("estimators.R")
source("generate_data_models.R")

results.linear <- data.frame("sample.size" = c(),
                             "estimate" = c(),
                             "estimator" = c(),
                             "subset" = c(),
                             "nuisance" = c())

different_subset_tested <- c("extended",
                             "minimal")

for (sample.size in c(300, 1000, 3000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    
    # generate a simulation
    a_simulation <- generate_simulation_leborgne_foucher(n = sample.size, setup = "complex", all_covariates_output = FALSE)
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", c(3, 11, 15, 7, 19, 2, 10, 8, 18, 14, 6))
        X_outcome <- paste0("X.", c(3, 11, 15, 7, 19, 2, 10, 8, 18, 14, 6))
        
      } else if (method == "smart"){
        X_treatment <- paste0("X.", c(3, 11, 15, 7, 19))
        X_outcome <- paste0("X.", c(3, 11, 15, 7, 19, 2, 10, 8, 18, 14, 6))
        
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", c(3, 11, 15, 7, 19))
        X_outcome <- paste0("X.", c(3, 11, 15, 7, 19))
        
        
      } else {
        stop("error in subset.")
      }
      
      
      custom_ipw <- ipw_forest(covariates_names = X_treatment, 
                               dataframe = a_simulation,
                               outcome_name = "Y",
                               treatment_name = "A",
                               min.node.size.if.forest = 5,
                               number.of.trees = 1000)
      
      custom_aipw <- aipw_forest(covariates_names_vector_treatment = X_treatment,
                                 covariates_names_vector_outcome = X_outcome,
                                 dataframe = a_simulation,
                                 outcome_name = "Y",
                                 treatment_name = "A",
                                 min.node.size.if.forest = 1,
                                 n.folds = 10,
                                 outcome_nature = "binary",
                                 return.decomposition = FALSE)


      custom_aipw_linear <- aipw_linear(covariates_names_vector_treatment = X_treatment,
                                        covariates_names_vector_outcome = X_outcome,
                                        dataframe = a_simulation,
                                        n.folds = 0)
      
      
      wrapper_for_grf <- causal_forest_wrapper(X_outcome, dataframe = a_simulation)
      
      wrapper_for_tmle <- tmle_wrapper(covariates_names_vector = X_outcome, dataframe = a_simulation)
      
      wrapper_for_doubleML <- DoubleML_wrapper_fast(covariates_names_vector = X_outcome, dataframe = a_simulation, n.folds = 10, outcome_nature = "binary")
      
      estimates <- c(custom_ipw, 
                     custom_aipw["aipw"], 
                     custom_aipw["t.learner"], 
                     custom_aipw_linear["ipw"],
                     custom_aipw_linear["t.learner"],
                     custom_aipw_linear["aipw"],
                     wrapper_for_grf,
                     wrapper_for_tmle,
                     wrapper_for_doubleML)
      
      
      new.row <- data.frame("sample.size" = rep(sample.size, 9),
                            "estimate" = estimates,
                            "estimator" = c("ipw", "aipw", "g-formula", "ipw", "g-formula", "aipw", 'cf', 'tmle', "doubleML"),
                            "subset" = rep(method, 9),
                            "nuisance" = c('forest', 'forest', 'forest', "linear", "linear", "linear", "forest", "superlearner", "forest"))
      
                            
      results.linear <- rbind(results.linear, new.row)
    }
  }
}

write.csv(x=results.linear, file="./data/synthetic_foucher_complex.csv")
