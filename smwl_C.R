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
                             "term.F" = c())

different_subset_tested <- c("extended",
                             "smart",
                             "minimal")# Reproducibility
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
                             "term.F" = c())

different_subset_tested <- c("extended",
                             "smart",
                             "minimal")

fmla.outcome <-  formula("Y ~ X.1 + X.2 + X.3")
fmla.outcome.ext <-  formula("Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6")
fmla.treat <-  formula("A ~ X.1 + X.2 + X.3")
fmla.treat.ext  <-  formula("A ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6")

for (sample.size in c(100, 300, 1000, 3000, 10000, 30000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    
    # generate a simulation
    simulation <- generate_simulation_wager_nie(n = sample.size, setup = "A", all_covariates_output = TRUE)
    
    # estimate surface responses
    mu.1.model <- lm(fmla.outcome, 
                     data = simulation[simulation$A == 1, ])
    mu.0.model <- lm(fmla.outcome, 
                     data = simulation[simulation$A == 0, ])
    
    # estimate surface responses with extended set
    mu.1.model.ext  <- lm(fmla.outcome.ext, 
                          data = simulation[simulation$A == 1, ])
    mu.0.model.ext  <- lm(fmla.outcome.ext, 
                          data = simulation[simulation$A == 0, ])
    
    # estimate propensity scores without and with extended set
    propensity.model <- glm(fmla.treat, data = simulation, family="binomial")
    propensity.model.ext <- glm(fmla.treat.ext, data = simulation, family="binomial")
    
    # Predict
    mu.hat.1 <- predict(mu.1.model, newdata = simulation, type="response")
    mu.hat.0 <- predict(mu.0.model, newdata = simulation, type="response")
    
    mu.hat.1.ext <- predict(mu.1.model.ext, newdata = simulation, type="response")
    mu.hat.0.ext <- predict(mu.0.model.ext, newdata = simulation, type="response")
    
    e.hat <- predict(propensity.model, newdata = simulation, type="response")
    e.hat.ext <- predict(propensity.model.ext, newdata = simulation, type="response")
    
    
    # compute estimates
    Y = simulation$Y
    A = simulation$A
    
    ipw <- Y * (A/e.hat - (1-A)/(1-e.hat))
    
    ipw.ext <- Y * (A/e.hat.ext - (1-A)/(1-e.hat.ext))
    
    aipw <- (mu.hat.1 - mu.hat.0
             + A / e.hat * (Y -  mu.hat.1)
             - (1 - A) / (1 - e.hat) * (Y -  mu.hat.0))
    aipw.ext <- (mu.hat.1.ext - mu.hat.0.ext
                 + A / e.hat.ext * (Y -  mu.hat.1.ext)
                 - (1 - A) / (1 - e.hat.ext) * (Y -  mu.hat.0.ext))
    
    aipw.smart <- (mu.hat.1.ext - mu.hat.0.ext
                   + A / e.hat * (Y -  mu.hat.1.ext)
                   - (1 - A) / (1 - e.hat) * (Y -  mu.hat.0.ext))
    
    gformula <- mu.hat.1 - mu.hat.0
    gformula.ext <- mu.hat.1.ext - mu.hat.0.ext
    
    new.row <- data.frame("sample.size" = rep(sample.size, 2),
                          "estimate" = c(mean(ipw),
                                         mean(ipw.ext),
                                         mean(gformula),
                                         mean(gformula.ext),
                                         mean(aipw),
                                         mean(aipw.ext),
                                         mean(aipw.smart)),
                          "estimator" = c(rep("ipw", 2), rep("t-learner", 2), rep("aipw", 3)),
                          "subset" =  c(rep(c("minimal", "extended"),3), "smart"),
                          "nuisance" = rep("linear", 7),
                          "term.A" = rep(NA, 7), 
                          "term.B" = rep(NA, 7), 
                          "term.C" = rep(NA, 7),
                          "term.D" = rep(NA, 7),
                          "term.E" = rep(NA, 7), 
                          "term.F" = rep(NA, 7))
    results.linear <- rbind(results.linear, new.row)
    
  }
}

write.csv(x=results.linear, file="./data/A_linear.csv")


fmla.outcome <-  formula("Y ~ X.1 + X.2 + X.3")
fmla.outcome.ext <-  formula("Y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6")
fmla.treat <-  formula("A ~ X.1 + X.2 + X.3")
fmla.treat.ext  <-  formula("A ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6")

for (sample.size in c(100, 300, 1000, 3000, 10000, 30000)){
  print(paste0("Starting sample size ", sample.size))
  for (i in 1:50){
    
    # generate a simulation
    simulation <- generate_simulation_wager_nie(n = sample.size, setup = "C", all_covariates_output = TRUE)
    
    # estimate surface responses
    mu.1.model <- lm(fmla.outcome, 
                     data = simulation[simulation$A == 1, ])
    mu.0.model <- lm(fmla.outcome, 
                     data = simulation[simulation$A == 0, ])
    
    # estimate surface responses with extended set
    mu.1.model.ext  <- lm(fmla.outcome.ext, 
                          data = simulation[simulation$A == 1, ])
    mu.0.model.ext  <- lm(fmla.outcome.ext, 
                          data = simulation[simulation$A == 0, ])
    
    # estimate propensity scores without and with extended set
    propensity.model <- glm(fmla.treat, data = simulation, family="binomial")
    propensity.model.ext <- glm(fmla.treat.ext, data = simulation, family="binomial")
    
    # Predict
    mu.hat.1 <- predict(mu.1.model, newdata = simulation, type="response")
    mu.hat.0 <- predict(mu.0.model, newdata = simulation, type="response")
    
    mu.hat.1.ext <- predict(mu.1.model.ext, newdata = simulation, type="response")
    mu.hat.0.ext <- predict(mu.0.model.ext, newdata = simulation, type="response")
    
    e.hat <- predict(propensity.model, newdata = simulation, type="response")
    e.hat.ext <- predict(propensity.model.ext, newdata = simulation, type="response")
    
    
    # compute estimates
    Y = simulation$Y
    A = simulation$A
    
    ipw <- Y * (A/e.hat - (1-A)/(1-e.hat))
    
    ipw.ext <- Y * (A/e.hat.ext - (1-A)/(1-e.hat.ext))
    
    aipw <- (mu.hat.1 - mu.hat.0
             + A / e.hat * (Y -  mu.hat.1)
             - (1 - A) / (1 - e.hat) * (Y -  mu.hat.0))
    aipw.ext <- (mu.hat.1.ext - mu.hat.0.ext
                 + A / e.hat.ext * (Y -  mu.hat.1.ext)
                 - (1 - A) / (1 - e.hat.ext) * (Y -  mu.hat.0.ext))
    
    aipw.smart <- (mu.hat.1.ext - mu.hat.0.ext
                 + A / e.hat * (Y -  mu.hat.1.ext)
                 - (1 - A) / (1 - e.hat) * (Y -  mu.hat.0.ext))
    
    gformula <- mu.hat.1 - mu.hat.0
    gformula.ext <- mu.hat.1.ext - mu.hat.0.ext
    
    new.row <- data.frame("sample.size" = rep(sample.size, 7),
                            "estimate" = c(mean(ipw),
                                           mean(ipw.ext),
                                           mean(gformula),
                                           mean(gformula.ext),
                                           mean(aipw),
                                           mean(aipw.ext),
                                           mean(aipw.smart)),
                            "estimator" = c(rep("ipw", 2), rep("t-learner", 2), rep("aipw", 3)),
                            "subset" =  c(rep(c("minimal", "extended"),3), "smart"),
                            "nuisance" = rep("linear", 7),
                            "term.A" = rep(NA, 7), 
                            "term.B" = rep(NA, 7), 
                            "term.C" = rep(NA, 7),
                            "term.D" = rep(NA, 7),
                            "term.E" = rep(NA, 7), 
                            "term.F" = rep(NA, 7))
    results.linear <- rbind(results.linear, new.row)

  }
}

write.csv(x=results.linear, file="./data/C_linear.csv")
