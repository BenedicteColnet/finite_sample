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
library(grf)

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


for (sample.size in c(1000, 3000, 9000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:100){
    for (independence in c(FALSE)){
      
      # generate a simulation
      a_simulation <- generate_simulation_linear(n_obs = sample.size, independent_covariate = independence)
      
      # choose subset
      for (method in different_subset_tested){
        if (method == "extended"){
          X_treatment <- paste0("X.", 1:30)
          X_outcome <- paste0("X.", 1:30)
        } else if (method == "minimal"){
          X_treatment <- paste0("X.", 1:6)
          X_outcome <- paste0("X.", 1:6)
        } else if (method == "smart"){
          X_treatment <- paste0("X.", 1:6)
          X_outcome <- paste0("X.", 1:30)
        } else {
          stop("error in subset.")
        }
        
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

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/test.csv")