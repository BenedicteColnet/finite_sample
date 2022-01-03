# IPW with deep forest
ipw_forest <- function(covariates_names,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        n.folds = 2,
                        min.node.size.if.forest = 1,
                        return.decomposition = FALSE) {
  
  n <- nrow(dataframe)
  Y = dataframe[, outcome_name]
  W = dataframe[, treatment_name]

  propensity.model <- probability_forest(dataframe[, covariates_names], 
                                          W, 
                                          num.trees=1000, 
                                          min.node.size=min.node.size.if.forest)
  e.hat <- predict(propensity.model, data = X_t)$predictions[,2]
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))

  return(ipw)
}


# g-formula  with deep forest
t_learner_forest <- function(covariates_names,
                       dataframe,
                       outcome_name = "Y",
                       treatment_name = "A",
                       n.folds = 2,
                       min.node.size.if.forest = 1,
                       return.decomposition = FALSE) {
  
  n <- nrow(dataframe)
  Y = dataframe[, outcome_name]
  W = as.factor(dataframe[, treatment_name])
  
  # Estimation
  outcome.model.treated <- regression_forest(X = dataframe[ dataframe[,treatment_name] == 1, covariates_names_vector_treatment], 
                                             Y = dataframe[dataframe[,treatment_name] == 1, outcome_name], 
                                             num.trees = 1000, 
                                             min.node.size = min.node.size.if.forest)
  
  outcome.model.control <- regression_forest(X = dataframe[dataframe[,treatment_name] == 0, covariates_names_vector_treatment], 
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
                        return.decomposition = FALSE) {
  
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
                                                 num.trees = 1000, 
                                                 min.node.size = min.node.size.if.forest)
      
      outcome.model.control <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                 Y = dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name], 
                                                 num.trees = 1000, 
                                                 min.node.size = min.node.size.if.forest)
      
      propensity.model <- probability_forest(dataframe[-idx, covariates_names_vector_treatment], 
                                             as.factor(W[-idx]), 
                                             num.trees = 1000, 
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
  
  
  # # replace extreme values if necessary
  # e.hat <- replace(e.hat, e.hat < 0.01, 0.01)
  # e.hat <- replace(e.hat, e.hat > 0.99, 0.99)
  
  # compute estimates
  aipw = mean(mu.hat.1 - mu.hat.0
              + W / e.hat * (Y -  mu.hat.1)
              - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  
  g_formula = mean(mu.hat.1) - mean(mu.hat.0)
  
  
  if(!return.decomposition){
    res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw)
  } else {
    
    # warning, this loop requires the dataframe to contain extra-info such as mu_1 and true e
    term.A <- mean( (dataframe$mu_1 - mu.hat.1) * (1 - (dataframe$A /dataframe$e))  ) 
    term.B <- mean( dataframe$A * (dataframe$Y - dataframe$mu_1) * ((1/e.hat) - (1/dataframe$e))  )
    term.C <- mean( dataframe$A * ( (1/e.hat) - (1/dataframe$e) ) * (mu.hat.1-dataframe$mu_1) )
    
    term.D <- mean( (dataframe$mu_0 - mu.hat.0) * (1 - ( (1 - dataframe$A) / (1 - dataframe$e)))  ) 
    term.E <- mean( (1 - dataframe$A) * (dataframe$Y - dataframe$mu_0) * ((1/ (1 - e.hat)) - (1/ (1 - dataframe$e) ) )  )
    term.F <- mean(  (1 - dataframe$A) * ( (1/ (1 - e.hat)) - (1/ (1 - dataframe$e)) ) * (mu.hat.0-dataframe$mu_0) )
    
    
    
    res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw,
            "term.A" = term.A, "term.B" = term.B, "term.C" = term.C,
            "term.D" = term.D, "term.E" = term.E, "term.F" = term.F)
  }
  
  
  
  return(res)
}

# custom aipw with forest but double ml
aipw_forest_double_ml <- function(covariates_names_vector_treatment, 
                                  covariates_names_vector_outcome,
                                  dataframe,
                                  outcome_name = "Y",
                                  treatment_name = "A",
                                  n.folds = 3,
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
  
  indices <- split(seq(n), sort(seq(n) %% n.folds))
  # cross-fitting of nuisance parameters
  
  idx1 <- indices[1]$`0`
  idx2 <- indices[2]$`1`
  idx3 <- indices[3]$`2`
  
  
  # Estimation
  outcome.model.treated_k1 <- regression_forest(X = dataframe[idx1 & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                Y = dataframe[idx1 & dataframe[,treatment_name] == 1, outcome_name], 
                                                num.trees = 500, 
                                                min.node.size = min.node.size.if.forest)
  
  
  outcome.model.treated_k2 <- regression_forest(X = dataframe[idx2 & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                Y = dataframe[idx2 & dataframe[,treatment_name] == 1, outcome_name], 
                                                num.trees = 500, 
                                                min.node.size = min.node.size.if.forest)
  
  outcome.model.treated_k3 <- regression_forest(X = dataframe[idx3 & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                Y = dataframe[idx3 & dataframe[,treatment_name] == 1, outcome_name], 
                                                num.trees = 500, 
                                                min.node.size = min.node.size.if.forest)
  
  outcome.model.control_k1 <- regression_forest(X = dataframe[idx1 & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                Y = dataframe[idx1 & dataframe[,treatment_name] == 0, outcome_name], 
                                                num.trees = 500, 
                                                min.node.size = min.node.size.if.forest)
  
  outcome.model.control_k2 <- regression_forest(X = dataframe[idx2 & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                Y = dataframe[idx2 & dataframe[,treatment_name] == 0, outcome_name], 
                                                num.trees = 500, 
                                                min.node.size = min.node.size.if.forest)
  
  outcome.model.control_k3 <- regression_forest(X = dataframe[idx3 & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                Y = dataframe[idx3 & dataframe[,treatment_name] == 0, outcome_name], 
                                                num.trees = 500, 
                                                min.node.size = min.node.size.if.forest)
  
  
  propensity.model_k1 <- probability_forest(dataframe[idx1, covariates_names_vector_treatment], 
                                            as.factor(W[idx1]), 
                                            num.trees = 500, 
                                            min.node.size=min.node.size.if.forest)
  
  propensity.model_k2 <- probability_forest(dataframe[idx2, covariates_names_vector_treatment], 
                                            as.factor(W[idx2]), 
                                            num.trees = 500, 
                                            min.node.size=min.node.size.if.forest)
  
  propensity.model_k3 <- probability_forest(dataframe[idx3, covariates_names_vector_treatment], 
                                            as.factor(W[idx3]), 
                                            num.trees = 500, 
                                            min.node.size=min.node.size.if.forest)
  
  
  # Prediction
  mu.hat.1[idx1] <- predict(outcome.model.treated_k2, newdata = xt1[idx1,])$predictions
  mu.hat.0[idx1] <- predict(outcome.model.control_k2, newdata = xt0[idx1,])$predictions
  e.hat[idx1] <- predict(propensity.model_k3, newdata = X_t[idx1,])$predictions[,2]
  
  mu.hat.1[idx2] <- predict(outcome.model.treated_k3, newdata = xt1[idx2,])$predictions
  mu.hat.0[idx2] <- predict(outcome.model.control_k3, newdata = xt0[idx2,])$predictions
  e.hat[idx2] <- predict(propensity.model_k1, newdata = X_t[idx2,])$predictions[,2]
  
  mu.hat.1[idx3] <- predict(outcome.model.treated_k1, newdata = xt1[idx3,])$predictions
  mu.hat.0[idx3] <- predict(outcome.model.control_k1, newdata = xt0[idx3,])$predictions
  e.hat[idx3] <- predict(propensity.model_k2, newdata = X_t[idx3,])$predictions[,2]
  
  
  
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
  
  # Matrix of (transformed) covariates used to estimate and predict e(X) = P[W=1|X]
  fmla.x <- formula(paste(" ~ 0 + ", paste0("bs(", covariates_names_vector_outcome, ", df=3)", collapse=" + ")))
  XX <- model.matrix(fmla.x, dataframe)
  
  # (Optional) Not penalizing the main effect (the coefficient on W)
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


# Custom AIPW with super learning for nuisance parameters
aipw_ML <- function(covariates_names_vector_treatment,
                    covariates_names_vector_outcome,
                    dataframe,
                    outcome_name = "Y",
                    treatment_name = "A",
                    n.folds = 2,
                    sl_libs_outcome = c("SL.mean", "SL.lm"),
                    sl_libs_treatment =  c("SL.glm", "SL.mean")){
  
  n_obs <- nrow(dataframe)
  
  # Cross-fitted estimates of E[Y|X,W=1], E[Y|X,W=0] and e(X) = P[W=1|X]
  mu.hat.1 <- rep(NA, n_obs)
  mu.hat.0 <- rep(NA, n_obs)
  e.hat <- rep(NA, n_obs)
  
  # Useful quantities for afterward
  Y <- dataframe[,outcome_name]
  W <- dataframe[,treatment_name]
  X <- dataframe[,covariates_names_vector_treatment]
  
  
  if (n.folds > 1){
    
    indices <- split(seq(n_obs), sort(seq(n_obs) %% n.folds))
    
    for (idx in indices) {
      
      X_treated <- dataframe[-idx & dataframe[,treatment_name] == 1, covariates_names_vector_outcome]
      Y_treated <- dataframe[-idx & dataframe[,treatment_name] == 1, outcome_name]
      X_control <- dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_outcome]
      Y_control <- dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name]
      
      # Fit super learner on the set -idx
      #print("START mu.1 - super learning")
      mu.1.model <- SuperLearner(Y = Y_treated, 
                                 X = X_treated, 
                                 family = gaussian(), 
                                 SL.library = sl_libs_outcome) 
      #print("START mu.0 - super learning")
      mu.0.model <- SuperLearner(Y = Y_control, 
                                 X = X_control, 
                                 family = gaussian(), 
                                 SL.library = sl_libs_outcome) 
      #print("START e(x) - super learning")
      propensity.model <- SuperLearner(Y = W[-idx], 
                                       X = X[-idx,], 
                                       family = binomial(), 
                                       SL.library = sl_libs_treatment) 
      
      # Predict with cross-fitting
      mu.hat.1[idx] <- predict(mu.1.model, 
                               newdata = dataframe[idx,  covariates_names_vector_outcome])$pred
      mu.hat.0[idx] <- predict(mu.0.model, 
                               newdata = dataframe[idx,  covariates_names_vector_outcome])$pred
      e.hat[idx] <- predict(propensity.model, 
                            newdata = dataframe[idx, covariates_names_vector_treatment])$pred
      
    }
  } else if (n.folds == 0 | n.folds == 1){
    
    X_treated <- dataframe[dataframe[,treatment_name] == 1, covariates_names_vector_outcome]
    Y_treated <- dataframe[dataframe[,treatment_name] == 1, outcome_name]
    X_control <- dataframe[dataframe[,treatment_name] == 0, covariates_names_vector_outcome]
    Y_control <- dataframe[dataframe[,treatment_name] == 0, outcome_name]
    
    
    # Fit super learner
    mu.1.model <- SuperLearner(Y = Y_treated, 
                               X = X_treated, 
                               family = gaussian(), 
                               SL.library = sl_libs_outcome) 
    
    mu.0.model <- SuperLearner(Y = Y_control, 
                               X = X_control, 
                               family = gaussian(), 
                               SL.library = sl_libs_outcome) 
    
    propensity.model <- SuperLearner(Y = Y, 
                                     X = X, 
                                     family = binomial(), 
                                     SL.library = sl_libs_treatment) 
    
    # Predict 
    mu.hat.1 <- predict(mu.1.model, newdata = dataframe[,covariates_names_vector_outcome])$pred
    mu.hat.0 <- predict(mu.0.model, newdata = dataframe[, covariates_names_vector_outcome])$pred
    e.hat <- predict(propensity.model, newdata = dataframe[, covariates_names_vector_treatment])$pred
    
    
  } else {
    stop("n.fold must be a positive integer")
  }
  
  # replace extreme values if necessary
  e.hat <- replace(e.hat, e.hat < 0.01, 0.01)
  e.hat <- replace(e.hat, e.hat > 0.99, 0.99)  
  
  ## T-learner
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

aipw_wrapped <- function(covariates_names_vector_treatment,
                         covariates_names_vector_outcome,
                         dataframe,
                         outcome_name = "Y",
                         treatment_name = "A",
                         n.folds = 2,
                         sl_libs_outcome = c("SL.mean", "SL.lm"),
                         sl_libs_treatment =  c("SL.glm", "SL.mean")){
  
  AIPW_SL <-    aipw_wrapper(Y = dataframe[, outcome_name],
                             A = dataframe[, treatment_name],
                             W = NULL,
                             W.Q = dataframe[,covariates_names_vector_outcome], 
                             W.g = dataframe[,covariates_names_vector_treatment], 
                             Q.SL.library = sl_libs_outcome,
                             g.SL.library = sl_libs_treatment,
                             k_split = n.folds,
                             verbose=FALSE)
  return(AIPW_SL$result[3,1])
}


aipw_forest_two_fold <- function(covariates_names_vector_treatment, 
                                 covariates_names_vector_outcome,
                                 dataframe,
                                 outcome_name = "Y",
                                 treatment_name = "A",
                                 min.node.size.if.forest = 1) {
  
  # cut in two the data set
  n <- nrow(dataframe)
  idx_to_fit <- 1:floor(n/3)
  #idx_to_fit_e <- (floor(n/3)+1):(2*floor(n/3))
  idx_to_estimate <- (2*floor(n/3) + 1):n
  
  # Estimation
  outcome.model.treated <- regression_forest(X = dataframe[idx_to_fit & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                             Y = dataframe[idx_to_fit & dataframe[,treatment_name] == 1, outcome_name], 
                                             num.trees = 500, 
                                             min.node.size = min.node.size.if.forest)
  
  outcome.model.control <- regression_forest(X = dataframe[idx_to_fit & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                             Y = dataframe[idx_to_fit & dataframe[,treatment_name] == 0, outcome_name], 
                                             num.trees = 500, 
                                             min.node.size = min.node.size.if.forest)
  
  
  propensity.model <- probability_forest(dataframe[idx_to_fit, covariates_names_vector_treatment], 
                                         as.factor(dataframe[idx_to_fit, treatment_name]), 
                                         num.trees = 500, 
                                         min.node.size=min.node.size.if.forest)
  
  # Prediction
  mu.hat.1 <- predict(outcome.model.treated, 
                      newdata = dataframe[idx_to_estimate, covariates_names_vector_outcome])$predictions
  mu.hat.0 <- predict(outcome.model.control, 
                      newdata = dataframe[idx_to_estimate,covariates_names_vector_outcome])$predictions
  e.hat <- predict(propensity.model, 
                   newdata = dataframe[idx_to_estimate,covariates_names_vector_treatment])$predictions[,2]
  
  W <- dataframe[idx_to_estimate, treatment_name]
  Y <- dataframe[idx_to_estimate, outcome_name]
  
  # # replace extreme values if necessary
  # e.hat <- replace(e.hat, e.hat < 0.01, 0.01)
  # e.hat <- replace(e.hat, e.hat > 0.99, 0.99) 
  
  # compute estimates
  aipw = mean(mu.hat.1 - mu.hat.0
              + W / e.hat * (Y -  mu.hat.1)
              - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  g_formula = mean(mu.hat.1) - mean(mu.hat.0)
  res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw)
  
  return(res)
}

aipw_forest_three_fold <- function(covariates_names_vector_treatment, 
                                   covariates_names_vector_outcome,
                                   dataframe,
                                   outcome_name = "Y",
                                   treatment_name = "A",
                                   min.node.size.if.forest = 1) {
  
  # cut in three the data set
  n <- nrow(dataframe)
  idx_to_fit_mu <- 1:floor(n/3)
  idx_to_fit_e <- (floor(n/3)+1):(2*floor(n/3))
  idx_to_estimate <- (2*floor(n/3) + 1):n
  
  # Estimation
  outcome.model.treated <- regression_forest(X = dataframe[idx_to_fit_mu & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                             Y = dataframe[idx_to_fit_mu & dataframe[,treatment_name] == 1, outcome_name], 
                                             num.trees = 500, 
                                             min.node.size = min.node.size.if.forest)
  
  outcome.model.control <- regression_forest(X = dataframe[idx_to_fit_mu & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                             Y = dataframe[idx_to_fit_mu & dataframe[,treatment_name] == 0, outcome_name], 
                                             num.trees = 500, 
                                             min.node.size = min.node.size.if.forest)
  
  
  propensity.model <- probability_forest(dataframe[idx_to_fit_e, covariates_names_vector_treatment], 
                                         as.factor(dataframe[idx_to_fit_e, treatment_name]), 
                                         num.trees = 500, 
                                         min.node.size=min.node.size.if.forest)
  
  # Prediction
  mu.hat.1 <- predict(outcome.model.treated, 
                      newdata = dataframe[idx_to_estimate, covariates_names_vector_outcome])$predictions
  mu.hat.0 <- predict(outcome.model.control, 
                      newdata = dataframe[idx_to_estimate,covariates_names_vector_outcome])$predictions
  e.hat <- predict(propensity.model, 
                   newdata = dataframe[idx_to_estimate,covariates_names_vector_treatment])$predictions[,2]
  
  W <- dataframe[idx_to_estimate, treatment_name]
  Y <- dataframe[idx_to_estimate, outcome_name]
  
  # # replace extreme values if necessary
  # e.hat <- replace(e.hat, e.hat < 0.01, 0.01)
  # e.hat <- replace(e.hat, e.hat > 0.99, 0.99) 
  
  # compute estimates
  aipw = mean(mu.hat.1 - mu.hat.0
              + W / e.hat * (Y -  mu.hat.1)
              - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))
  g_formula = mean(mu.hat.1) - mean(mu.hat.0)
  res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw)
  
  return(res)
}



