
# custom AIPW with forest
aipw_forest <- function(covariates_names_vector_treatment, 
                        covariates_names_vector_outcome,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        n.folds = 2,
                        min.node.size.if.forest = 5) {
  
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
      outcome.model.treated <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                 Y = dataframe[-idx & dataframe[,treatment_name] == 1, outcome_name], 
                                                 num.trees = 500, 
                                                 min.node.size = min.node.size.if.forest)
      
      outcome.model.control <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                 Y = dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name], 
                                                 num.trees = 500, 
                                                 min.node.size = min.node.size.if.forest)
      
      propensity.model <- probability_forest(dataframe[-idx, covariates_names_vector_treatment], 
                                             as.factor(W[-idx]), 
                                             num.trees = 500, 
                                             min.node.size=min.node.size.if.forest)
      
      # Prediction
      mu.hat.1[idx] <- predict(outcome.model.treated, newdata = xt1[idx,])$predictions
      mu.hat.0[idx] <- predict(outcome.model.control, newdata = xt0[idx,])$predictions
      e.hat[idx] <- predict(propensity.model, newdata = X_t[idx,])$predictions[,2]
      
    }
    
  } else if (n.folds == 0 | n.folds == 1){
    
    # Estimation
    outcome.model.treated <- regression_forest(X = dataframe[ dataframe[,treatment_name] == 1, covariates_names_vector_treatment], 
                                               Y = dataframe[dataframe[,treatment_name] == 1, outcome_name], 
                                               num.trees = 500, 
                                               min.node.size = min.node.size.if.forest)
    
    outcome.model.control <- regression_forest(X = dataframe[dataframe[,treatment_name] == 0, covariates_names_vector_treatment], 
                                               Y = dataframe[dataframe[,treatment_name] == 0, outcome_name], 
                                               num.trees = 500, 
                                               min.node.size = min.node.size.if.forest)
    
    propensity.model <- probability_forest(dataframe[, covariates_names_vector_outcome], 
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
  
  
  # replace extreme values if necessary
  e.hat <- replace(e.hat, e.hat < 0.01, 0.01)
  e.hat <- replace(e.hat, e.hat > 0.99, 0.99)
  
  # compute estimates
  aipw = mean(mu.hat.1 - mu.hat.0
              + W / e.hat * (Y -  mu.hat.1)
              - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  
  g_formula = mean(mu.hat.1) - mean(mu.hat.0)
  
  res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw)
  
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

  # replace extreme values if necessary
  e.hat <- replace(e.hat, e.hat < 0.01, 0.01)
  e.hat <- replace(e.hat, e.hat > 0.99, 0.99)  
    
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
                         nuisance = "glmnet",
                         n.folds = 2,
                         automate = FALSE){
  
  if (nuisance == "glmnet"){
    SL.library<- c("SL.glm", "SL.glmnet", "SL.glm.interaction")
  } else if (nuisance == "forest"){
    SL.library <- c("SL.glm", "SL.ranger")
  } else if (nuisance == "linear"){
    SL.library.outcome <- c("SL.lm")
    SL.library.treatment <- c("SL.glm")
  } else {
    stop("error in nuisance - TMLE")
  }
  
  # ranger takes way more time
  
  
  TMLE <- tmle(Y = dataframe[,outcome_name],
               A = dataframe[,treatment_name],
               W = dataframe[,covariates_names_vector],
               family = "gaussian",
               Q.SL.library = SL.library.outcome,
               g.SL.library = SL.library.treatment,
               V = n.folds)
  
  return(TMLE$estimates$ATE$psi)
}

# smart AIPW with splines
aipw_splines <- function(covariates_names_vector_treatment,
                         covariates_names_vector_outcome,
                         dataframe,
                         outcome_name = "Y",
                         treatment_name = "A",
                         n.folds = 2){
  
  n <- nrow(dataframe)
  
  # Preparing data
  W <- dataframe[,treatment_name]
  Y <- dataframe[,outcome_name]
  
  fmla.xw <- formula(paste("~ 0 +", paste0("bs(", covariates_names_vector_treatment, ", df=3)", "*", treatment_name, collapse=" + ")))
  XW <- model.matrix(fmla.xw, dataframe)
  data.1 <- dataframe
  data.1[,treatment_name] <- 1
  XW1 <- model.matrix(fmla.xw, data.1)  # setting W=1
  data.0 <- dataframe
  data.0[,treatment_name] <- 0
  XW0 <- model.matrix(fmla.xw, data.0)  # setting W=0
  fmla.x <- formula(paste(" ~ 0 + ", paste0("bs(", covariates_names_vector_outcome, ", df=3)", collapse=" + ")))
  XX <- model.matrix(fmla.x, dataframe)
  
  penalty.factor <- rep(1, ncol(XW))
  penalty.factor[colnames(XW) == treatment_name] <- 0
  
  mu.hat.1 <- rep(NA, n)
  mu.hat.0 <- rep(NA, n)
  e.hat <- rep(NA, n)
  
  if (n.folds > 1){
    indices <- split(seq(n), sort(seq(n) %% n.folds))
    
    for (idx in indices) {
      # Estimate outcome model and propensity models
      # Note how cross-validation is done (via cv.glmnet) within cross-fitting!
      outcome.model<- cv.glmnet(x=XW[-idx,], y=Y[-idx], family="gaussian", penalty.factor=penalty.factor)
      propensity.model <- cv.glmnet(x=XX[-idx,], y=W[-idx], family="binomial")
      
      # Predict with cross-fitting
      mu.hat.1[idx] <- predict(outcome.model, newx=XW1[idx,], type="response")
      mu.hat.0[idx] <- predict(outcome.model, newx=XW0[idx,], type="response")
      e.hat[idx] <- predict(propensity.model, newx=XX[idx,], type="response")
    }
    
  } else if (n.folds == 0 | n.folds == 1){
    
    # Estimation
    # Note how cross-validation is done (via cv.glmnet) within cross-fitting!
    outcome.model<- cv.glmnet(x=XW, y=Y, family="gaussian", penalty.factor=penalty.factor)
    propensity.model <- cv.glmnet(x=XX, y=W, family="binomial")
    
    # Prediction
    mu.hat.1 <- predict(outcome.model, newx=XW1, type="response")
    mu.hat.0 <- predict(outcome.model, newx=XW0, type="response")
    e.hat <- predict(propensity.model, newx=XX, type="response")
    
  } else {
    stop("n.fold must be a positive integer")
  }
  
  
  
  # Compute the summand in AIPW estimator
  aipw <- (mu.hat.1 - mu.hat.0
           + W / e.hat * (Y -  mu.hat.1)
           - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  aipw <- mean(aipw)
  
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  
  t.learner = mean(mu.hat.1) - mean(mu.hat.0)
  
  res = c("ipw" = ipw, "t.learner" = t.learner, "aipw" = aipw)
  
  return(res)
  
}