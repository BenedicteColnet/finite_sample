# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(glmnet)
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
                             "independence" = c())

different_subset_tested <- c("all.covariates",
                             "all.covariates.wo.instruments",
                             "smart",
                             "minimal.set")

for (sample.size in c(300, 1000, 3000, 9000, 30000, 100000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    for (independence in c(TRUE, FALSE)){
      
      # generate a simulation
      a_simulation <- generate_simulation_linear(n_obs = sample.size, independent_covariate = independence)
      
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
        print(paste0("custom aipw, rep", i))
        custom_aipw_5 <- aipw_linear(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 5)
        custom_aipw_20 <- aipw_linear(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 20)
        print(paste0("wrapper aipw", i))
        wrapper_5 <- aipw_wrapped(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 5)
        wrapper_20 <- aipw_wrapped(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 20)
        print(paste0("wrapper tmle", i))
        wrapper_tmle <- tmle_wrapper(X_treatment, dataframe = a_simulation)
        
        new.row <- data.frame("sample.size" = rep(sample.size, 9),
                              "estimate" = c(custom_aipw_5["ipw"],
                                             custom_aipw_5["t.learner"],
                                             custom_aipw_5["aipw"],
                                             custom_aipw_20["ipw"],
                                             custom_aipw_20["t.learner"],
                                             custom_aipw_20["aipw"],
                                             wrapper_5,
                                             wrapper_20,
                                             wrapper_tmle),
                                "estimator" = c("ipw",
                                                "t-learner",
                                                "aipw",
                                                "ipw",
                                                "t-learner",
                                                "aipw",
                                                "wrapper aipw",
                                                "wrapper aipw",
                                                "wrapper tmle"),
                                "subset" = rep(method, 9),
                                "simulation" = rep("linear", 9),
                                "cross-fitting" = c(5,5,5,20,20,20,5,20,NA),
                                "independence" = rep(independence, 9))
          
        results.linear <- rbind(results.linear, new.row)
      
      }
    }
  }
}

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/2021-10-26-linear-constant-ate.csv")