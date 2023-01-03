#### IPW with forest with deep trees

ipw_forest <- function(covariates_names,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        min.node.size.if.forest = 1,
                        number.of.trees = 1000) {
  
  n <- nrow(dataframe)
  Y = dataframe[, outcome_name]
  W = dataframe[, treatment_name]

  propensity.model <- probability_forest(dataframe[, covariates_names], 
                                         as.factor(W), 
                                         num.trees = number.of.trees, 
                                         min.node.size = min.node.size.if.forest)
  
  e.hat <- predict(propensity.model, data = X_t)$predictions[,2]
  ipw = mean(Y * (W/e.hat - (1-W)/(1-e.hat)))

  return(ipw)
}



#-------------



#### g-formula with forest with deep trees

t_learner_forest <- function(covariates_names,
                       dataframe,
                       outcome_name = "Y",
                       treatment_name = "A",
                       min.node.size.if.forest = 5, #like in grf package and in particular causal forests
                       number.of.trees = 200, #like in grf package and in particular causal forests
                       return.decomposition = FALSE) {
  
  n <- nrow(dataframe)
  Y = dataframe[, outcome_name]
  W = dataframe[, treatment_name]
  
  # Estimation
  outcome.model.treated <- regression_forest(X = dataframe[ dataframe[,treatment_name] == 1, covariates_names], 
                                             Y = dataframe[dataframe[,treatment_name] == 1, outcome_name], 
                                             num.trees = number.of.trees, 
                                             min.node.size = min.node.size.if.forest)
  
  outcome.model.control <- regression_forest(X = dataframe[dataframe[,treatment_name] == 0, covariates_names], 
                                             Y = dataframe[dataframe[,treatment_name] == 0, outcome_name], 
                                             num.trees = number.of.trees, 
                                             min.node.size = min.node.size.if.forest)

  # Prediction
  mu.hat.1 <- predict(outcome.model.treated, data = xt1)$predictions
  mu.hat.0 <- predict(outcome.model.control, data = xt0)$predictions
  
  t_learner = mean(mu.hat.1) - mean(mu.hat.0)
  
  return(t_learner)
}



#-------------


### Custom AIPW with forest 

# custom AIPW with forest
aipw_forest <- function(covariates_names_vector_treatment, 
                        covariates_names_vector_outcome,
                        dataframe,
                        outcome_name = "Y",
                        treatment_name = "A",
                        n.folds = 2,
                        outcome_nature = "continuous",
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
      propensity.model <- probability_forest(dataframe[-idx, covariates_names_vector_treatment], 
                                             as.factor(W[-idx]), 
                                             num.trees = 1000, 
                                             min.node.size=min.node.size.if.forest)
      
      
      if(outcome_nature == "continuous"){
        
        
        outcome.model.treated <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                   Y = dataframe[-idx & dataframe[,treatment_name] == 1, outcome_name], 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest)
        
        outcome.model.control <- regression_forest(X = dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                   Y = dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name], 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest)
        
        
      } else if (outcome_nature == "binary"){
        
        
        outcome.model.treated <- probability_forest(X = dataframe[-idx & dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                   Y =  as.factor(dataframe[-idx & dataframe[,treatment_name] == 1, outcome_name]), 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest)
        
        outcome.model.control <- probability_forest(X = dataframe[-idx & dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                   Y =  as.factor(dataframe[-idx & dataframe[,treatment_name] == 0, outcome_name]), 
                                                   num.trees = 1000, 
                                                   min.node.size = min.node.size.if.forest)
        
        
        
      } else {
        
        print("Error, outcome_nature has to be continuous or binary")
        break
        
      }
      
      # Prediction
      
      if(outcome_nature == "continuous"){
        
        mu.hat.1[idx] <- predict(outcome.model.treated, newdata = xt1[idx,])$predictions
        mu.hat.0[idx] <- predict(outcome.model.control, newdata = xt0[idx,])$predictions
        
      } else { 
        
        mu.hat.1[idx] <- predict(outcome.model.treated, newdata = xt1[idx, covariates_names_vector_outcome])$predictions[,2]
        mu.hat.0[idx] <- predict(outcome.model.control, newdata = xt0[idx, covariates_names_vector_outcome])$predictions[,2]
        
      }
      
      
      e.hat[idx] <- predict(propensity.model, newdata = X_t[idx,])$predictions[,2]
      
    }
    
  } else if (n.folds == 0 | n.folds == 1){
    
    # Estimation
    
    if(outcome_nature == "continuous"){
      
      outcome.model.treated <- regression_forest(X = dataframe[ dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                 Y = dataframe[dataframe[,treatment_name] == 1, outcome_name], 
                                                 num.trees = 1000, 
                                                 min.node.size = min.node.size.if.forest)
      
      outcome.model.control <- regression_forest(X = dataframe[dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                 Y = dataframe[dataframe[,treatment_name] == 0, outcome_name], 
                                                 num.trees = 1000, 
                                                 min.node.size = min.node.size.if.forest)
      
    } else if (outcome_nature == "binary"){
      
      outcome.model.treated <- probability_forest(X = dataframe[dataframe[,treatment_name] == 1, covariates_names_vector_outcome], 
                                                  Y =  as.factor(dataframe[dataframe[,treatment_name] == 1, outcome_name]), 
                                                  num.trees = 1000, 
                                                  min.node.size = min.node.size.if.forest)
      
      outcome.model.control <- probability_forest(X = dataframe[dataframe[,treatment_name] == 0, covariates_names_vector_outcome], 
                                                  Y =  as.factor(dataframe[dataframe[,treatment_name] == 0, outcome_name]), 
                                                  num.trees = 1000, 
                                                  min.node.size = min.node.size.if.forest)
      
    } else {
      
      print("Error, outcome_nature has to be continuous or binary")
      break
      
    }
    
    
    propensity.model <- probability_forest(dataframe[, covariates_names_vector_treatment], 
                                           as.factor(W), 
                                           num.trees=1000, 
                                           min.node.size=min.node.size.if.forest)
    
    # Prediction
    
    if(outcome_nature == "continuous"){
      
      mu.hat.1 <- predict(outcome.model.treated, data = xt1)$predictions
      mu.hat.0 <- predict(outcome.model.control, data = xt0)$predictions
      
    } else {
      
      mu.hat.1 <- predict(outcome.model.treated, data = xt1)$predictions[,2]
      mu.hat.0 <- predict(outcome.model.control, data = xt0)$predictions[,2]
      
    }
    
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
    term.C <- - mean( dataframe$A * ( (1/e.hat) - (1/dataframe$e) ) * (mu.hat.1-dataframe$mu_1) )
    
    term.D <- - mean( (dataframe$mu_0 - mu.hat.0) * (1 - ( (1 - dataframe$A) / (1 - dataframe$e)))  ) 
    term.E <- - mean( (1 - dataframe$A) * (dataframe$Y - dataframe$mu_0) * ((1/ (1 - e.hat)) - (1/ (1 - dataframe$e) ) )  )
    term.F <- mean(  (1 - dataframe$A) * ( (1/ (1 - e.hat)) - (1/ (1 - dataframe$e)) ) * (mu.hat.0-dataframe$mu_0) )
    
    res = c("ipw" = ipw, "t.learner" = g_formula, "aipw" = aipw,
            "term.A" = term.A, "term.B" = term.B, "term.C" = term.C,
            "term.D" = term.D, "term.E" = term.E, "term.F" = term.F, 
            "semi.oracle.aipw" = semi.oracle.aipw)
  }
  
  
  
  return(res)
}


#-------------





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

binned_ipw <- function(covariates_names_vector, 
                       dataframe,
                       outcome_name = "Y",
                       treatment_name = "A",
                       nb.bin = 10){
  
  # better have a data driven number of bins?
  for (covariate.name in covariates_names_vector){
    covariate <- dataframe[, covariate.name]
    
    # if continuous
    if(any(as.integer(covariate) != covariate) || length(unique(covariate)) > 2){
      deciles <- quantcut(covariate, seq(0, 1, by = 1/nb.bin))
      dataframe[, covariate.name] <- as.factor(deciles)
    }
  }
  
  
  e.hat <- dataframe %>% 
    group_by(across(covariates_names_vector)) %>%
    summarise(e.hat = mean(A))
  
  
  final <- dataframe
  final <- merge(final, e.hat, by = covariates_names_vector)
  gamma <- ((final$Y * final$A)/final$e.hat) - ((final$Y * (1-final$A))/(1-final$e.hat))
  return(mean(gamma, na.rm = TRUE))
}



causal_forest_wrapper <- function(covariates_names_vector, 
                                  dataframe,
                                  outcome_name = "Y",
                                  treatment_name = "A"){
  
  forest <- causal_forest(
    X=dataframe[, covariates_names_vector],
    Y=dataframe[, outcome_name],
    W=dataframe[, treatment_name],
    num.trees = 1000)

  # Estimate the conditional average treatment effect on the full sample (CATE).
  forest.ate <- average_treatment_effect(forest, target.sample = "all")

  res = forest.ate[[1]]
  
  return(res)
}


causal_forest_wrapper_different_covariate_set <- function(covariates_names_vector_treatment,
                                                          covariates_names_vector_outcome,
                                                          dataframe,
                                                          outcome_name = "Y",
                                                          treatment_name = "A"){
  
  forest.W <- regression_forest(dataframe[ ,covariates_names_vector_treatment], dataframe[, treatment_name], tune.parameters = "all")
  W.hat <- predict(forest.W)$predictions
  
  forest.Y <- regression_forest(dataframe[ ,covariates_names_vector_outcome], dataframe[, outcome_name], tune.parameters = "all")
  Y.hat <- predict(forest.Y)$predictions
  
  forest.Y.varimp <- variable_importance(forest.Y)
  
  selected.vars <- which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)
  
  tau.forest <- causal_forest(dataframe[, selected.vars], dataframe[, outcome_name], dataframe[, treatment_name],
                              W.hat = W.hat, Y.hat = Y.hat,
                              tune.parameters = "all")
  
  # Estimate the conditional average treatment effect on the full sample (CATE).
  forest.ate <- average_treatment_effect(tau.forest, target.sample = "all")
  
  res = forest.ate[[1]]
  
  return(res)
  
}


tmle_wrapper <- function(covariates_names_vector, 
                         dataframe,
                         outcome_name = "Y",
                         treatment_name = "A",
                         n.folds = 2,
                         automate = FALSE,
                         sl_libs_outcome = c('SL.ranger', "SL.glmnet", "SL.xgboost", "SL.glm", "SL.lm"),
                         sl_libs_treatment = c('SL.ranger', "SL.glmnet", "SL.xgboost", "SL.glm", "SL.lm")){
  
  
  
  TMLE <- tmle(Y = dataframe[,outcome_name],
               A = dataframe[,treatment_name],
               W = dataframe[,covariates_names_vector],
               family = "gaussian",
               Q.SL.library = sl_libs_outcome,
               g.SL.library = sl_libs_treatment,
               V = n.folds)
  
  return(TMLE$estimates$ATE$psi)
}


DoubleML_wrapper <- function(covariates_names_vector, 
                             dataframe,
                             outcome_name = "Y",
                             treatment_name = "A",
                             n.folds = 10,
                             outcome_nature = "continuous"){
  
  
  # suppress messages during fitting
  lgr::get_logger("mlr3")$set_threshold("warn")
  
  # Initialize DoubleMLData (data-backend of DoubleML)
  data_dml_base = DoubleMLData$new(dataframe,
                                   y_col = outcome_name,
                                   d_cols =treatment_name,
                                   x_cols = covariates_names_vector)
  
  
  # Trees
  trees = lrn("regr.rpart", cp = 0.0047, minsplit = 203)
  trees_class = lrn("classif.rpart", cp = 0.0042, minsplit = 104)
  
  dml_plr_tree = DoubleMLPLR$new(data_dml_base,
                                 ml_l = trees,
                                 ml_m = trees_class,
                                 n_folds = n.folds)
  dml_plr_tree$fit()
  
  
  # Random Forest
  mtry.est = floor(sqrt(length(covariates_names_vector)))
  
  if (outcome_nature == "continuous"){
    randomForest_outcome = lrn("regr.ranger", max.depth = 100, mtry = mtry.est, min.node.size = 1)
  } else if (outcome_nature == "binary"){
    randomForest_outcome = lrn("classif.ranger", max.depth = 100, mtry = mtry.est, min.node.size = 1)
  }
  
  randomForest_treatment = lrn("classif.ranger", max.depth = 100, mtry = mtry.est, min.node.size = 1)
  
  
  dml_plr_forest = DoubleMLPLR$new(data_dml_base,
                                   ml_l = randomForest_outcome,
                                   ml_m = randomForest_treatment,
                                   n_folds = 10)
  dml_plr_forest$fit()
  
  
  # Lasso
  if (outcome_nature == "continuous"){
    lasso_outcome = lrn("regr.cv_glmnet", nfolds = 10, s = "lambda.min")
  } else if (outcome_nature == "binary"){
    lasso_outcome = lrn("classif.cv_glmnet", nfolds = 10, s = "lambda.min")
  }
  
  lasso_treatment = lrn("classif.cv_glmnet", nfolds = 10, s = "lambda.min")
  
  # Initialize DoubleMLPLR model
  dml_plr_lasso = DoubleMLPLR$new(data_dml_base,
                                  ml_l = lasso_outcome,
                                  ml_m = lasso_treatment,
                                  n_folds = 10)
  dml_plr_lasso$fit()
  
  
  # Boosted trees
  if (outcome_nature == "continuous"){
    boost_outcome = lrn("regr.xgboost", objective = "reg:squarederror", eta = 0.1, nrounds = 35)
  } else if (outcome_nature == "binary"){
    boost_outcome = lrn("classif.xgboost", objective = "binary:logistic", eval_metric = "logloss",eta = 0.1, nrounds = 34)
  }

  boost_treatment = lrn("classif.xgboost",
                    objective = "binary:logistic", eval_metric = "logloss",
                    eta = 0.1, nrounds = 34)
  

  dml_plr_boost = DoubleMLPLR$new(data_dml_base,
                                  ml_l = boost_outcome,
                                  ml_m = boost_treatment,
                                  n_folds = 10)
  dml_plr_boost$fit()
  
  
  # Most general model with random forests
  if (outcome_nature == "continuous"){
    randomForest_outcome = lrn("regr.ranger")
  } else if (outcome_nature == "binary"){
    randomForest_outcome = llrn("classif.ranger")
  }

  randomForest_treatment = lrn("classif.ranger")
  
  dml_irm_forest = DoubleMLIRM$new(data_dml_base,
                                   ml_g = randomForest_outcome,
                                   ml_m = randomForest_treatment,
                                   trimming_threshold = 0.01,
                                   n_folds = n.folds)
  
  # Set nuisance-part specific parameters
  dml_irm_forest$set_ml_nuisance_params("ml_g0", "A", list(max.depth = 100, mtry = mtry.est, min.node.size = 1))
  
  dml_irm_forest$fit()
  
  
  confints = rbind(dml_plr_lasso$confint(), dml_plr_forest$confint(),
                   dml_plr_tree$confint(), dml_plr_boost$confint(), dml_irm_forest$confint())
  estimates = c(dml_plr_lasso$coef, dml_plr_forest$coef,
                dml_plr_tree$coef, dml_plr_boost$coef, dml_irm_forest$coef)
  result_plr = data.frame("model" = c(rep("PLR", 4),  "IRM"),
                          "ML" = c("glmnet", "ranger", "rpart", "xgboost", "ranger"),
                          "Estimate" = estimates,
                          "lower" = confints[,1],
                          "upper" = confints[,2])
  return(result_plr)
}


# only compute the IRM with random forest
DoubleML_wrapper_fast <- function(covariates_names_vector, 
                             dataframe,
                             outcome_name = "Y",
                             treatment_name = "A",
                             n.folds = 10,
                             outcome_nature = "continuous"){
  
  
  # Initialize DoubleMLData (data-backend of DoubleML)
  data_dml_base = DoubleMLData$new(dataframe,
                                   y_col = outcome_name,
                                   d_cols =treatment_name,
                                   x_cols = covariates_names_vector)
  
  # suppress messages during fitting
  lgr::get_logger("mlr3")$set_threshold("warn")
  
  # compute an estimation of mtry
  mtry.est = floor(sqrt(length(covariates_names_vector)))
  
  if (outcome_nature == "continuous"){
    randomForest_outcome = lrn("regr.ranger")
  } else if (outcome_nature == "binary"){
    randomForest_outcome = lrn("classif.ranger")
  }
  
  randomForest_treatment = lrn("classif.ranger")
  
  dml_irm_forest = DoubleMLIRM$new(data_dml_base,
                                   ml_g = randomForest_outcome,
                                   ml_m = randomForest_treatment,
                                   trimming_threshold = 0.01,
                                   n_folds = n.folds)
  
  # Set nuisance-part specific parameters with very deep trees
  dml_irm_forest$set_ml_nuisance_params("ml_g0", "A", list(max.depth = 1000, mtry = mtry.est, min.node.size = 1))
  
  dml_irm_forest$fit()
  
  return(dml_irm_forest$coef[[1]])
  
}
