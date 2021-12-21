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
                     "bias.e" = c(),
                     "term.A" = c(),
                     "term.B" = c(),
                     "term.C" = c(),
                     "AIPW" = c(),
                     "subset" = c())

different_subset_tested <- c("extended",
                             "smart",
                             "minimal")


for (sample.size in c(150, 300, 500, 1000, 2000, 5000, 10000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    # generate a simulation
    simulation <- generate_simulation_wager_nie(n = sample.size, setup = "D.bis")
    
    # choose subset
    for (method in different_subset_tested){
      if (method == "extended"){
        X_treatment <- paste0("X.", 1:5)
        X_outcome <- paste0("X.", 1:5)
      } else if (method == "smart"){
        X_treatment <- paste0("X.", 1:2)
        X_outcome <- paste0("X.", 1:5)
      } else if (method == "minimal"){
        X_treatment <- paste0("X.", 1:2)
        X_outcome <- paste0("X.", 1:2)
      } else {
        stop("error in subset.")
      }
    
      # fit models
      outcome.model.treated <-  ranger(Y ~ .,  
                                       num.trees = 500, 
                                       max.depth = NULL,
                                       min.node.size = 1, 
                                       data = simulation[simulation$A == 1, c("Y", X_outcome)])
      outcome.model.control <-  ranger(Y ~ .,  
                                       num.trees = 500, 
                                       max.depth = NULL,
                                       min.node.size = 1, 
                                       data = simulation[simulation$A == 0, c("Y", X_outcome)])
      
      propensity.model <- probability_forest(simulation[, X_treatment], 
                                             as.factor(simulation[, "A"]), 
                                             num.trees = 500, 
                                             min.node.size=1)
      
      
      # prediction and estimation
      simulation.to.estimate <- generate_simulation_wager_nie(n = 10000, setup = "D", all_covariates_output = TRUE)
      mu.hat.1 <- predict(outcome.model.treated, simulation.to.estimate[, X_outcome])$predictions
      bias.mu.1 <- mean(mu.hat.1-simulation.to.estimate$mu_1)
      mu.hat.0 <- predict(outcome.model.control, simulation.to.estimate[, X_outcome])$predictions
      bias.mu.0 <- mean(mu.hat.0-simulation.to.estimate$mu_0)
      e.hat <- predict(propensity.model, 
                       newdata = simulation.to.estimate[,X_treatment])$predictions[,2]
      bias.e <- mean(e.hat-simulation.to.estimate$e)
      
      term.A <- mean( (simulation.to.estimate$mu_1 - mu.hat.1) * (1 - (simulation.to.estimate$A /simulation.to.estimate$e))  ) 
      term.B <- mean( (simulation.to.estimate$Y_1 - simulation.to.estimate$mu_1) * ((1/e.hat) - (1/simulation.to.estimate$e))  )
      term.C <- mean( (e.hat-simulation.to.estimate$e) * (mu.hat.1-simulation.to.estimate$mu_1) )
      
      W <- simulation.to.estimate$A
      Y <- simulation.to.estimate$Y
      aipw.on.second.fold <- mean(mu.hat.1 - mu.hat.0
                                  + W / e.hat * (Y -  mu.hat.1)
                                  - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
      
      new_row <- data.frame("sample.size" = sample.size,
                            "bias.mu.1" = bias.mu.1,
                            "bias.mu.0" = bias.mu.0,
                            "bias.e" = bias.e,
                            "term.A" = term.A,
                            "term.B" = term.B,
                            "term.C" = term.C,
                            "AIPW" = aipw.on.second.fold,
                            "subset" = method)
      
      results <- rbind(results, new_row)
    }
  }
}

write.csv(x=results, file="./data/2021-12-16-biasbis.csv")

