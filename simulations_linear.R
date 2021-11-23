# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(mvtnorm) # rmvnorm
library(glmnet)

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

different_subset_tested <- c("all",
                             "outcome.and.noise",
                             "outcome",
                             "smart",
                             "minimal.set")

for (sample.size in c(300, 1000, 3000, 9000, 30000, 100000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:500){
    for (independence in c(TRUE, FALSE)){
      
      # generate a simulation
      a_simulation <- generate_simulation_linear(n_obs = sample.size, independent_covariate = independence)
      
      # choose subset
      for (method in different_subset_tested){
        if (method == "all"){
          X_treatment <- paste0("X.", 1:50)
          X_outcome <- paste0("X.", 1:50)
        } else if (method == "outcome.and.noise"){
          X_treatment <- paste0("X.", 2:50)
          X_outcome <- paste0("X.", 2:50)
        } else if (method == "outcome"){
          X_treatment <- paste0("X.", 2:30)
          X_outcome <- paste0("X.", 2:30)
        } else if (method == "minimal.set"){
          X_treatment <- paste0("X.", 2:6)
          X_outcome <- paste0("X.", 2:6)
        } else if (method == "smart"){
          X_treatment <- paste0("X.", 2:6)
          X_outcome <- paste0("X.", 2:30)
        } else {
          stop("error in subset.")
        }
        
        custom_aipw_2 <- aipw_linear(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 2)
        
        new.row <- data.frame("sample.size" = rep(sample.size, 3),
                              "estimate" = c(custom_aipw_2["ipw"],
                                             custom_aipw_2["t.learner"],
                                             custom_aipw_2["aipw"]),
                                "estimator" = rep(c("ipw",
                                                "t-learner",
                                                "aipw"),1),
                                "subset" = rep(method, 3),
                                "simulation" = rep("linear.constant.cate", 3),
                                "cross-fitting" = c(2,2,2),
                                "independence" = rep(independence, 3),
                                "nuisance" = rep("linear", 3))
          
        results.linear <- rbind(results.linear, new.row)
      
      }
    }
  }
}

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/2021-11-23-parametric-variance-simple.csv")