# IPW with forest with deep trees
ipw_forest <- function(covariates_names,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        min.node.size.if.forest = 1,
                        return.decomposition = FALSE) {
  
  n <- nrow(dataframe)
  Y = dataframe[, outcome_name]
  W = dataframe[, treatment_name]

  propensity.model <- probability_forest(dataframe[, covariates_names], 
                                         as.factor(W), 
                                          num.trees = 1000, 
                                          min.node.size = min.node.size.if.forest)
  e.hat <- predict(propensity.model, data = X_t)$predictions[,2]
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))

  return(ipw)
}


# g-formula with forest with deep trees
t_learner_forest <- function(covariates_names,
                       dataframe,
                       outcome_name = "Y",
                       treatment_name = "A",
                       min.node.size.if.forest = 1,
                       return.decomposition = FALSE) {
  
  n <- nrow(dataframe)
  Y = dataframe[, outcome_name]
  W = dataframe[, treatment_name]
  
  # Estimation
  outcome.model.treated <- regression_forest(X = dataframe[ dataframe[,treatment_name] == 1, covariates_names], 
                                             Y = dataframe[dataframe[,treatment_name] == 1, outcome_name], 
                                             num.trees = 1000, 
                                             min.node.size = min.node.size.if.forest)
  
  outcome.model.control <- regression_forest(X = dataframe[dataframe[,treatment_name] == 0, covariates_names], 
                                             Y = dataframe[dataframe[,treatment_name] == 0, outcome_name], 
                                             num.trees = 1000, 
                                             min.node.size = min.node.size.if.forest)

  # Prediction
  mu.hat.1 <- predict(outcome.model.treated, data = xt1)$predictions
  mu.hat.0 <- predict(outcome.model.control, data = xt0)$predictions
  
  t_learner = mean(mu.hat.1) - mean(mu.hat.0)
  
  return(t_learner)
}

# custom AIPW with forest
aipw_forest <- function(covariates_names_vector_treatment, 
                        covariates_names_vector_outcome,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        n.folds = 2,
                        min.node.size.if.forest = 1,
                        return.decomposition = FALSE,
                        with.weights = FALSE) {
  
  n <- nrow(dataframe)
  
  t0 = rep(0, n)
  t1 = rep(1, n)
  
  X_t <- dataframe[, covariates_names_vector_treatment]
  X_o <- dataframe[, covariates_names_vector_outcome]
  W <- dataframe[, treatment_name]
  Y <- dataframe[, outcome_name]
  
  xt <- cbind(X_o, W)
  xt0 <- cbind(X_o)
  xt1 <- cbind(X_o)
  
  mu.hat.1 <- rep(NA, n)
  mu.hat.0 <- rep(NA, n)
  e.hat <- rep(NA, n)
  
  if (n.folds > 1){
    indices <- split(seq(n), sort(seq(n) %% n.folds))
    
    # cross-fitting of nuisance parameters
    for (idx in indices) {
      
      
      # Estimation
      propensity.model <- probability_forest(dataframe[-idx, covariates_names_vector_treatment], 
                                             as.factor(W[-idx]), 
                                             num.trees = 1000, 
                                             min.node.size=min.node.size.if.forest)
      
      if (!with.weights){
        outcome.model.treated <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                   Y = dataframe[-idx & dataframe[,treatment_name] == 1, outcome_name], 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest)
        
        outcome.model.control <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                   Y = dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name], 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest)
      } else {
        
        weights <- dataframe[-idx, c(covariates_names_vector_treatment, treatment_name)]
        weights$w <- predict(propensity.model, newdata = dataframe[,covariates_names_vector_treatment])$predictions[,2]
        
        outcome.model.treated <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                   Y = dataframe[-idx & dataframe[,treatment_name] == 1, outcome_name], 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest,
                                                   sample.weights = 1/weights[weights$A == 1, "w"])
        
        outcome.model.control <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                   Y = dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name], 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest,
                                                   sample.weights = 1/(1-weights[weights$A == 0, "w"]))
        
      }
      
      
      
      
      # Prediction
      mu.hat.1[idx] <- predict(outcome.model.treated, newdata = xt1[idx,])$predictions
      mu.hat.0[idx] <- predict(outcome.model.control, newdata = xt0[idx,])$predictions
      e.hat[idx] <- predict(propensity.model, newdata = X_t[idx,])$predictions[,2]
      
    }
    
  } else if (n.folds == 0 | n.folds == 1){
    
    # Estimation
    outcome.model.treated <- regression_forest(X = dataframe[ dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                               Y = dataframe[dataframe[,treatment_name] == 1, outcome_name], 
                                               num.trees = 500, 
                                               min.node.size = min.node.size.if.forest)
    
    outcome.model.control <- regression_forest(X = dataframe[dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                               Y = dataframe[dataframe[,treatment_name] == 0, outcome_name], 
                                               num.trees = 500, 
                                               min.node.size = min.node.size.if.forest)
    
    propensity.model <- probability_forest(dataframe[, covariates_names_vector_treatment], 
                                           as.factor(W), 
                                           num.trees=500, 
                                           min.node.size=min.node.size.if.forest)
    
    # Prediction
    mu.hat.1 <- predict(outcome.model.treated, data = xt1)$predictions
    mu.hat.0 <- predict(outcome.model.control, data = xt0)$predictions
    e.hat <- predict(propensity.model, data = X_t)$predictions[,2]
    
  } else {
    stop("n.fold must be a positive integer")
  }

  
  # compute estimates
  aipw = mean(mu.hat.1 - mu.hat.0
              + W / e.hat * (Y -  mu.hat.1)
              - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  
  g_formula = mean(mu.hat.1) - mean(mu.hat.0)
  
  
  if(!return.decomposition){
    res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw)
  } else {
    
    semi.oracle.aipw <- mean(dataframe$mu_1 - dataframe$mu_0
                             + W / e.hat * (Y -  dataframe$mu_1)
                             - (1 - W) / (1 - e.hat) * (Y -  dataframe$mu_0))
    
    # warning, this loop requires the dataframe to contain extra-info such as mu_1 and true e
    term.A <- mean( (dataframe$mu_1 - mu.hat.1) * (1 - (dataframe$A /dataframe$e))  ) 
    term.B <- mean( dataframe$A * (dataframe$Y - dataframe$mu_1) * ((1/e.hat) - (1/dataframe$e))  )
    term.C <- mean( dataframe$A * ( (1/e.hat) - (1/dataframe$e) ) * (mu.hat.1-dataframe$mu_1) )
    
    term.D <- mean( (dataframe$mu_0 - mu.hat.0) * (1 - ( (1 - dataframe$A) / (1 - dataframe$e)))  ) 
    term.E <- mean( (1 - dataframe$A) * (dataframe$Y - dataframe$mu_0) * ((1/ (1 - e.hat)) - (1/ (1 - dataframe$e) ) )  )
    term.F <- mean(  (1 - dataframe$A) * ( (1/ (1 - e.hat)) - (1/ (1 - dataframe$e)) ) * (mu.hat.0-dataframe$mu_0) )
    
    res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw,
            "term.A" = term.A, "term.B" = term.B, "term.C" = term.C,
            "term.D" = term.D, "term.E" = term.E, "term.F" = term.F, 
            "semi.oracle.aipw" = semi.oracle.aipw)
  }
  
  
  
  return(res)
}


# Custom AIPW with linear and logit model
aipw_linear <- function(covariates_names_vector_treatment,
                        covariates_names_vector_outcome,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        n.folds = 2){
  
  n_obs <- nrow(dataframe)
  
  # Prepare formulas
  fmla.treatment <- formula(paste0(treatment_name,"~."))
  fmla.outcome <- formula(paste0(outcome_name,"~."))
  
  # Cross-fitted estimates of E[Y|X,W=1], E[Y|X,W=0] and e(X) = P[W=1|X]
  mu.hat.1 <- rep(NA, n_obs)
  mu.hat.0 <- rep(NA, n_obs)
  e.hat <- rep(NA, n_obs)
  
  if (n.folds > 1){
    
    indices <- split(seq(n_obs), sort(seq(n_obs) %% n.folds))
    
    for (idx in indices) {
      
      # Fit model on the set -idx
      mu.1.model <- lm(fmla.outcome, 
                       data = dataframe[-idx & dataframe[, treatment_name] == 1, c(outcome_name, covariates_names_vector_outcome)])
      mu.0.model <- lm(fmla.outcome, 
                       data = dataframe[-idx & dataframe[, treatment_name] == 0, c(outcome_name, covariates_names_vector_outcome)])
      propensity.model <- glm(fmla.treatment, data = dataframe[-idx, c(treatment_name, covariates_names_vector_treatment)], family="binomial")
      
      # Predict with cross-fitting
      mu.hat.1[idx] <- predict(mu.1.model, newdata = dataframe[idx, c(outcome_name, covariates_names_vector_outcome)], type="response")
      mu.hat.0[idx] <- predict(mu.0.model, newdata = dataframe[idx, c(outcome_name, covariates_names_vector_outcome)], type="response")
      e.hat[idx] <- predict(propensity.model, newdata = dataframe[idx,  c(treatment_name, covariates_names_vector_treatment)], type="response")
      
    }
  } else if (n.folds == 0 | n.folds == 1){
    
    # Fit model on all observations
    mu.1.model <- lm(fmla.outcome, 
                     data = dataframe[dataframe[, treatment_name] == 1, c(outcome_name, covariates_names_vector_outcome)])
    mu.0.model <- lm(fmla.outcome, 
                     data = dataframe[dataframe[, treatment_name] == 0, c(outcome_name, covariates_names_vector_outcome)])
    propensity.model <- glm(fmla.treatment, data = dataframe[, c(treatment_name, covariates_names_vector_treatment)], family="binomial")
    
    # Predict with same observations
    mu.hat.1 <- predict(mu.1.model)
    mu.hat.0 <- predict(mu.0.model)
    e.hat <- predict(propensity.model)
  } else {
    stop("n.fold must be a positive integer")
  }

    
  # Compute the summand in AIPW estimator
  W <- dataframe[, treatment_name]
  Y <- dataframe[, outcome_name]
  
  ## T- learner
  t.learner <- mean(mu.hat.1) - mean(mu.hat.0)
  
  ## AIPW
  aipw <- (mu.hat.1 - mu.hat.0
           + W / e.hat * (Y -  mu.hat.1)
           - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  
  aipw <- mean(aipw)
  
  ## IPW
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  
  res = c("ipw" = ipw, "t.learner" = t.learner, "aipw" = aipw)
  
  return(res)
  
}


# Custom AIPW with logit models
aipw_logit <- function(covariates_names_vector_treatment,
                        covariates_names_vector_outcome,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        n.folds = 2){
  
  n_obs <- nrow(dataframe)
  
  # Prepare formulas
  fmla.treatment <- formula(paste0(treatment_name,"~."))
  fmla.outcome <- formula(paste0(outcome_name,"~."))
  
  # Cross-fitted estimates of E[Y|X,W=1], E[Y|X,W=0] and e(X) = P[W=1|X]
  mu.hat.1 <- rep(NA, n_obs)
  mu.hat.0 <- rep(NA, n_obs)
  e.hat <- rep(NA, n_obs)
  
  if (n.folds > 1){
    
    indices <- split(seq(n_obs), sort(seq(n_obs) %% n.folds))
    
    for (idx in indices) {
      
      # Fit model on the set -idx
      mu.1.model <- glm(fmla.outcome, 
                       data = dataframe[-idx & dataframe[, treatment_name] == 1, c(outcome_name, covariates_names_vector_outcome)], family="binomial")
      mu.0.model <- glm(fmla.outcome, 
                       data = dataframe[-idx & dataframe[, treatment_name] == 0, c(outcome_name, covariates_names_vector_outcome)], family="binomial")
      propensity.model <- glm(fmla.treatment, data = dataframe[-idx, c(treatment_name, covariates_names_vector_treatment)], family="binomial")
      
      # Predict with cross-fitting
      mu.hat.1[idx] <- predict(mu.1.model, newdata = dataframe[idx, c(outcome_name, covariates_names_vector_outcome)], type="response")
      mu.hat.0[idx] <- predict(mu.0.model, newdata = dataframe[idx, c(outcome_name, covariates_names_vector_outcome)], type="response")
      e.hat[idx] <- predict(propensity.model, newdata = dataframe[idx,  c(treatment_name, covariates_names_vector_treatment)], type="response")
      
    }
  } else if (n.folds == 0 | n.folds == 1){
    
    # Fit model on all observations
    mu.1.model <- glm(fmla.outcome, 
                     data = dataframe[dataframe[, treatment_name] == 1, c(outcome_name, covariates_names_vector_outcome)], family="binomial")
    mu.0.model <- glm(fmla.outcome, 
                     data = dataframe[dataframe[, treatment_name] == 0, c(outcome_name, covariates_names_vector_outcome)], family="binomial")
    propensity.model <- glm(fmla.treatment, data = dataframe[, c(treatment_name, covariates_names_vector_treatment)], family="binomial")
    
    # Predict with same observations
    mu.hat.1 <- predict(mu.1.model)
    mu.hat.0 <- predict(mu.0.model)
    e.hat <- predict(propensity.model)
  } else {
    stop("n.fold must be a positive integer")
  }
  
  
  # Compute the summand in AIPW estimator
  W <- dataframe[, treatment_name]
  Y <- dataframe[, outcome_name]
  
  ## T- learner
  t.learner <- mean(mu.hat.1)/mean(mu.hat.0)
  
  ## AIPW
  m1 <- (W / e.hat) * (Y -  mu.hat.1) + mu.hat.1 
  m0 <- ((1-W) / (1-e.hat)) * (Y -  mu.hat.0) + mu.hat.0
  
  aipw <- mean(m1)/mean(m0)
  
  ## IPW
  ipw = (mean(W*Y / e.hat)) / (mean((1-W)*Y/(1-e.hat)))
  
  res = c("ipw" = ipw, "t.learner" = t.learner, "aipw" = aipw)
  
  return(res)
  
}


causal_forest_wrapper <- function(covariates_names_vector, 
                                  dataframe,
                                  outcome_name = "Y",
                                  treatment_name = "A"){
  
  forest <- causal_forest(
    X=dataframe[, covariates_names_vector],  
    W=dataframe[, treatment_name],
    Y=dataframe[, outcome_name],
    num.trees = 100)
  forest.ate <- average_treatment_effect(forest)
  return(forest.ate[[1]])
}

tmle_wrapper <- function(covariates_names_vector, 
                         dataframe,
                         outcome_name = "Y",
                         treatment_name = "A",
                         n.folds = 2,
                         automate = FALSE,
                         sl_libs_outcome = c("SL.mean", "SL.lm"),
                         sl_libs_treatment = c("SL.glm", "SL.mean")){
  
  
  TMLE <- tmle(Y = dataframe[,outcome_name],
               A = dataframe[,treatment_name],
               W = dataframe[,covariates_names_vector],
               family = "gaussian",
               Q.SL.library = sl_libs_outcome,
               g.SL.library = sl_libs_treatment,
               V = n.folds)
  
  return(TMLE$estimates$ATE$psi)
}