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

for (sample.size in c(100, 300, 1000, 3000, 9000, 30000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    # generate a simulation
    a_simulation <- generate_simulation_linear_constant_cate(n_obs = sample.size, independent_covariate = FALSE)
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "all.covariates"){
        X_treatment <- paste0("X.", 1:12)
        X_outcome <- paste0("X.", 1:12)
      } else if (method == "all.covariates.wo.instruments"){
        X_treatment <- paste0("X.", 4:12)
        X_outcome <- paste0("X.", 4:12)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 4:7)
        X_outcome <- paste0("X.", 4:10)
      } else if (method == "minimal.set"){
        X_treatment <- paste0("X.", 4:7)
        X_outcome <- paste0("X.", 4:7)
      } else {
        stop("error in subset.")
      }
      
      for (number_of_folds in c(2, 10)){
        custom_aipw <- aipw_linear(X_treatment, X_outcome, dataframe = a_simulation, n.folds = number_of_folds)
        tmle.estimate <- tmle_wrapper(covariates_names_vector = X_outcome, dataframe = a_simulation, nuisance = "linear", n.folds = number_of_folds)
        new.row <- data.frame("sample.size" = rep(sample.size, 4),
                              "estimate" = c(custom_aipw["ipw"],
                                             custom_aipw["t.learner"],
                                             custom_aipw["aipw"],
                                             tmle.estimate),
                              "estimator" = c("ipw",
                                              "t-learner",
                                              "aipw",
                                              "tmle"),
                              "subset" = rep(method, 4),
                              "simulation" = rep("linear", 4),
                              "cross-fitting" = rep(number_of_folds, 4))
        
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/2021-10-23-linear-nonindep.csv")