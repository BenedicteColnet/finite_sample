# Reproducibility
set.seed(123)

# Libraries
library(dplyr) # case_when and others
library(tidyr) # pivot
library(grf)
source("estimators.R")
source("generate_data_models.R")

# load data
PATH <- "../data/cohort_2019_imputed_after_2_composite_covariate.RData"
load(PATH)

data_depp <- data_depp[data_depp$Taille_Classe > 8,]
data_depp$Treatment <- ifelse(data_depp$Taille_Classe < 13, 1, 0)

# Categ etab into one hot encoder
data_depp$REP <- ifelse(data_depp$Categ_Etab_CP == "REP", 1, 0)
data_depp$REPp <- ifelse(data_depp$Categ_Etab_CP == "REP+", 1, 0)
data_depp$Public <- ifelse(data_depp$Categ_Etab_CP == "Public", 1, 0)
data_depp$Private <- ifelse(data_depp$Categ_Etab_CP == "Private", 1, 0)

categ <- c("REP", "REPp", "Public", "Private")
minimal_set <- c(categ, "IPS_Etab_CP")
extended_set <- c(minimal_set, "Age_CP", "Sexe_Num", "T1_Math", "T1_Language")

results <- data_frame("estimator" = c(),
                      "estimate" = c(),
                      "sample.size" = c(),
                      "extended.set" = c())

for (sample.size in c(1000, 10000, 50000)){
  print(paste0("starting sample size ", str(sample.size)))
  for (i in 1:30){
    if(i == 10){
      print("starting 10")
    } else if (i == 20){
      print("starting 20")
    }
    workind_df <- data_depp[sample(nrow(data_depp), sample.size), ]
    estimate.with.minimal.set <- aipw_forest(covariates_names_vector_treatment = minimal_set,
                                             covariates_names_vector_outcome = minimal_set,
                                             dataframe = workind_df,
                                             outcome_name = "T3_Math",
                                             treatment_name = "Treatment",
                                             n.folds = 5,
                                             min.node.size.if.forest = 1)
    estimate.with.extended.set <- aipw_forest(covariates_names_vector_treatment = minimal_set,
                                              covariates_names_vector_outcome = extended_set,
                                              dataframe = workind_df,
                                              outcome_name = "T3_Math",
                                              treatment_name = "Treatment",
                                              n.folds = 5,
                                              min.node.size.if.forest = 1)
  }
  
  new_row <- data.frame("estimator" = rep(c("ipw", "t-learner", "aipw"),2),
                        "estimate" = c(estimate.with.minimal.set["ipw"],
                                       estimate.with.minimal.set["t.learner"],
                                       estimate.with.minimal.set["aipw"],
                                       estimate.with.extended.set["ipw"],
                                       estimate.with.extended.set["t.learner"],
                                       estimate.with.extended.set["aipw"]),
                        "sample.size" = rep(sample.size,6),
                        "extended.set" = c("no", "no", "no", "yes", "yes", "yes"))
  
  results <- rbind(results, new_row)
}


write.csv(x=results, file="./data/semi-synthetic-education.csv")