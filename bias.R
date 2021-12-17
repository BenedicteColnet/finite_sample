# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(tidyr) # pivot
library(mvtnorm) # rmvnorm
library(ranger)
library(grf)

source("estimators.R")
source("generate_data_models.R")

results <- data.frame("sample.size" = c(),
                     "bias.mu.1" = c(),
                     "bias.mu.0" = c(),
                     "bias.e" = c())


for (sample.size in c(500, 1000, 1500, 2000, 4000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:30){
    # generate a simulation
    simulation <- generate_simulation(n_obs = sample.size)
    
    # fit models
    outcome.model.treated <-  ranger(Y ~ .,  
                                     num.trees = 500, 
                                     mtry = 4,
                                     max.depth = NULL,
                                     min.node.size = 1, 
                                     data = simulation[simulation$A == 1, c("Y", paste0("X.", 1:12))])
    outcome.model.control <-  ranger(Y ~ .,  
                                     num.trees = 500, 
                                     mtry = 4,
                                     max.depth = NULL,
                                     min.node.size = 1, 
                                     data = simulation[simulation$A == 0, c("Y", paste0("X.", 1:12))])
    propensity.model <- probability_forest(simulation[, paste0("X.", 1:4)], 
                                           as.factor(simulation[, "A"]), 
                                           num.trees = 500, 
                                           min.node.size=1)
    
    
    # prediction and estimation
    simulation.to.estimate <- generate_simulation(n_obs = 10000, all_covariates_output = TRUE)
    mu.hat.1 <- predict(outcome.model.treated, simulation.to.estimate[, paste0("X.", 1:12)])$predictions
    bias.mu.1 <- mean(mu.hat.1-simulation.to.estimate$mu_1)
    mu.hat.1 <- predict(outcome.model.control, simulation.to.estimate[, paste0("X.", 1:12)])$predictions
    bias.mu.0 <- mean(mu.hat.0-simulation.to.estimate$mu_1)
    e.hat <- predict(propensity.model, 
                     newdata = simulation.to.estimate[,paste0("X.", 1:4)])$predictions[,2]
    bias.e <- mean(e.hat-simulation.to.estimate$e)
    
    new_row <- data.frame("sample.size" = sample.size,
                          "bias.mu.1" = bias.mu.1,
                          "bias.mu.0" = bias.mu.0,
                          "bias.e" = bias.e)
    
    results <- rbind(results, new_row)
    
  }
}
  


write.csv(x=results, file="./data/2021-12-16-bias.csv")