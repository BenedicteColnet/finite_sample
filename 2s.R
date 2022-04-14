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
                             "nuisance" = c(),
                             "term.A" = c(), 
                             "term.B" = c(), 
                             "term.C" = c(),
                             "term.D" = c(), 
                             "term.E" = c(), 
                             "term.F" = c(),
                             "setup" = c())

different_subset_tested <- c("extended",
                             "smart",
                             "minimal")

X_treatment <- paste0("X.", 1:2)
X_outcome <- paste0("X.", 1:6)

for (sample.size in c(100, 300, 1000, 3000, 10000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    
    for (wager in c("A", "C", "D")){

      
      # generate a simulation
      a_simulation <- generate_simulation_wager_nie(n = sample.size, setup = wager, all_covariates_output = TRUE)
      
      
      estimate.two.steps <- aipw_forest(X_treatment, 
                                        X_treatment, 
                  dataframe = a_simulation,
                  min.node.size.if.forest = 1,
                  n.folds = 0,
                  return.decomposition = TRUE,
                  with.weights = TRUE)
      
      new.row <- data.frame("sample.size" = sample.size,
                            "estimate" = estimate.two.steps["aipw"][[1]],
                            "estimator" = "aipw - full",
                            "subset" = "minimal",
                            "nuisance" = "forest",
                            "term.A" = NA, 
                            "term.B" = NA, 
                            "term.C" = NA,
                            "term.D" = NA, 
                            "term.E" = NA, 
                            "term.F" = NA,
                            "setup" = wager)
      
       new.row.gformula <- data.frame("sample.size" = sample.size,
                            "estimate" = estimate.two.steps["t.learner"][[1]],
                            "estimator" = "gformula - full",
                            "subset" = "minimal",
                            "nuisance" = "forest",
                            "term.A" = NA, 
                            "term.B" = NA, 
                            "term.C" = NA,
                            "term.D" = NA, 
                            "term.E" = NA, 
                            "term.F" = NA,
                            "setup" = wager)
       
       new.row.ipw <- data.frame("sample.size" = sample.size,
                                      "estimate" = estimate.two.steps["ipw"][[1]],
                                      "estimator" = "ipw",
                                      "subset" = "minimal",
                                      "nuisance" = "forest",
                                      "term.A" = NA, 
                                      "term.B" = NA, 
                                      "term.C" = NA,
                                      "term.D" = NA, 
                                      "term.E" = NA, 
                                      "term.F" = NA,
                                      "setup" = wager)
      
      results.linear <- rbind(results.linear, new.row, new.row.gformula, new.row.ipw)
      
    }
  }
}

write.csv(x=results.linear, file="./data/2s.csv")
