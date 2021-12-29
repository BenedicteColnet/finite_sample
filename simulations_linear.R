# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(ggplot2)
library(tidyr) # pivot
library(mvtnorm) # rmvnorm
library(glmnet)
library(grf)

source("estimators.R")
source("generate_data_models.R")

results.linear <- data.frame("sample.size" = c(),
                             "estimate" = c(),
                             "estimator" = c(),
                             "subset" = c(),
                             "simulation" = c(),
                             "independence" = c(),
                             "nuisance" = c())

different_subset_tested <- c("outcome",
                             "smart",
                             "minimal.set")

for (sample.size in c(500, 1000, 3000, 10000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    for (independence in c(TRUE)){
      
      # generate a simulation
      a_simulation <- generate_simulation_linear(n_obs = sample.size, independent_covariate = independence)
      
      # choose subset
      for (method in different_subset_tested){
        if (method == "outcome"){
          X_treatment <- paste0("X.", 1:30)
          X_outcome <- paste0("X.", 1:30)
        } else if (method == "minimal.set"){
          X_treatment <- paste0("X.", 1:6)
          X_outcome <- paste0("X.", 1:6)
        } else if (method == "smart"){
          X_treatment <- paste0("X.", 1:6)
          X_outcome <- paste0("X.", 1:30)
        } else {
          stop("error in subset.")
        }
        
        custom_aipw <- aipw_linear(X_treatment, 
                                   X_outcome, 
                                   dataframe = a_simulation, 
                                   n.folds = 2)
        custom_aipw_forest <- aipw_forest(X_treatment, 
                                          X_outcome, 
                                          dataframe = a_simulation, 
                                          n.folds = 2,
                                          min.node.size.if.forest = 1)
        
        new.row <- data.frame("sample.size" = rep(sample.size, 6),
                              "estimate" = c(custom_aipw["ipw"],
                                             custom_aipw["t.learner"],
                                             custom_aipw["aipw"],
                                             custom_aipw_forest["ipw"],
                                             custom_aipw_forest["t.learner"],
                                             custom_aipw_forest["aipw"]),
                                "estimator" = rep(c("ipw",
                                                "t-learner",
                                                "aipw"),2),
                                "subset" = rep(method, 6),
                                "simulation" = rep("linear.constant.cate", 6),
                                "independence" = rep(independence, 6),
                                "nuisance" = c(rep("linear", 3), rep("forest", 3)))
          
        results.linear <- rbind(results.linear, new.row)
      
      }
    }
  }
}

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/linear.csv")
