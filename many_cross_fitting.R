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

different_subset_tested <- c("outcome.wo.instruments")

for (sample.size in c(3000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    print(paste0("rep: ", i))
    for (independence in c(FALSE)){
      
      # generate a simulation
      a_simulation <- generate_simulation_linear(n_obs = sample.size, independent_covariate = independence)
      
      # choose subset
      for (method in different_subset_tested){
        if (method == "outcome.and.instruments"){
          X_treatment <- paste0("X.", 1:10)
          X_outcome <- paste0("X.", 1:10)
        } else if (method == "outcome.wo.instruments"){
          X_treatment <- paste0("X.", 4:10)
          X_outcome <- paste0("X.", 4:10)
        } else if (method == "smart"){
          X_treatment <- paste0("X.", 4:7)
          X_outcome <- paste0("X.", 4:10)
        } else if (method == "minimal.set"){
          X_treatment <- paste0("X.", 4:7)
          X_outcome <- paste0("X.", 4:7)
        } else {
          stop("error in subset.")
        }
        

        custom_aipw_100 <- aipw_forest(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 100, min.node.size.if.forest = 1)
        custom_aipw_linear <- aipw_linear(X_treatment, X_outcome, dataframe = a_simulation, n.folds = 2)
        
        
        new.row <- data.frame("sample.size" = rep(sample.size, 6),
                              "estimate" = c(custom_aipw_100["ipw"],
                                             custom_aipw_100["t.learner"],
                                             custom_aipw_100["aipw"],
                                             custom_aipw_linear["ipw"],
                                             custom_aipw_linear["t.learner"],
                                             custom_aipw_linear["aipw"]),
                              "estimator" = rep(c("ipw",
                                                  "t-learner",
                                                  "aipw"),2),
                              "subset" = rep(method, 6),
                              "simulation" = rep("linear.constant.cate", 6),
                              "cross-fitting" = c(100,100,100,2,2,2),
                              "independence" = rep(independence, 6),
                              "nuisance" = c("forest","forest","forest","linear", "linear", "linear"))
        
        results.linear <- rbind(results.linear, new.row)
        
      }
    }
  }
}

results.linear$sample.size <- as.factor(results.linear$sample.size)

write.csv(x=results.linear, file="./data/2021-11-03-linear-constant-ate-forest-vs-linear-many-crossfitting.csv")