# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
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
library(gtools) # quantcut

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
      
      
      ipw.10 <- binned_ipw(X_treatment,
                           dataframe = a_simulation,
                           nb.bin = 10)
      ipw.20 <- binned_ipw(X_treatment,
                           dataframe = a_simulation,
                           nb.bin = 20)
      ipw.30 <- binned_ipw(X_treatment,
                           dataframe = a_simulation,
                           nb.bin = 30)
      ipw.40 <- binned_ipw(X_treatment,
                           dataframe = a_simulation,
                           nb.bin = 40)
      
      estimate.cf <- causal_forest_wrapper(X_treatment, 
                                           dataframe = a_simulation)
      
      new.row <- data.frame("sample.size" = rep(sample.size, 4),
                            "estimate" = c(ipw.10, ipw.20, ipw.30, ipw.40),
                            "estimator" = paste0("ipw.bin.", c("10", "20", "30", "40")),
                            "subset" = rep("minimal", 4),
                            "nuisance" = rep("forest", 4),
                            "term.A" = rep(NA, 4), 
                            "term.B" = rep(NA, 4),  
                            "term.C" = rep(NA, 4), 
                            "term.D" = rep(NA, 4),  
                            "term.E" = rep(NA, 4),  
                            "term.F" = rep(NA, 4), 
                            "setup" = rep(wager, 4))
      
      new.row.bis <- data.frame("sample.size" = sample.size,
                            "estimate" = estimate.cf,
                            "estimator" = "causal forest",
                            "subset" = "minimal",
                            "nuisance" = "forest",
                            "term.A" = NA, 
                            "term.B" = NA, 
                            "term.C" = NA,
                            "term.D" = NA, 
                            "term.E" = NA, 
                            "term.F" = NA,
                            "setup" = wager)
      
      
      results.linear <- rbind(results.linear, new.row, new.row.bis)
      
    }
  }
}

write.csv(x=results.linear, file="./data/2s.csv")
