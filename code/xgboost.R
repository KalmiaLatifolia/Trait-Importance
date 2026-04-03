
# Section 2 of Full Protocol - Trait importance manuscript
# xgboost model
# written 12 Nov 2025
# updated 3 April 2026

# https://medium.com/@dhanyahari07/feature-selection-using-xgboost-f0622fb70c4d

library(xgboost)
library(dplyr)
library(caret)
library(SHAPforxgboost)
library(ggnewscale)
library(future)
library(future.apply)
library(ggpubr)
library(tidyr)
library(purrr)
library(report)
library(readxl)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full data ---------------------------------------------------------------

siteDetections_foliarTraits_BioCube <- readRDS("data/siteDetections_foliarTraits_BioCube_20260320.rds")

# load tidy names --------------------------------------------------------------

tidy <- read_excel("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TableS1_Biocube_var_description.xlsx",
                   range = cell_cols(1:4))

# set variable groups ----------------------------------------------------------

names(siteDetections_foliarTraits_BioCube)
species <- colnames(siteDetections_foliarTraits_BioCube)[4:97]
varSets <- list(
  spatVars  = tidy$Variable,
  notTraits = tidy$Variable[tidy$Category != "Traits"],
  notDist   = tidy$Variable[tidy$Category != "Disturbance"],
  notClim   = tidy$Variable[tidy$Category != "Climate"],
  notStr    = tidy$Variable[tidy$Category != "Structure"],
  notPheno  = tidy$Variable[tidy$Category != "Phenology"],
  notTerr   = tidy$Variable[tidy$Category != "Terrain"]
)

# make sure all column names in varSets match real variables
varSets <- lapply(varSets, function(x) x[x %in% colnames(siteDetections_foliarTraits_BioCube)])


# xgboost with SHAP and iterations in parallel ---------------------------------
# ------------------------------------------------------------------------------

# Use multiple cores
plan(multisession, workers = parallel::detectCores() - 1)  

# list parallelizable tasks (6580)
tasks <- expand.grid(
  varSet = names(varSets),
  species_i = seq_along(species),
  iter = 1:100
)

# make a fxn -------------------------------------------------------------------
run_task <- function(vs, i, iter) {
  set.seed(123 + iter)
  
  # choose variable set
  spatVars <- siteDetections_foliarTraits_BioCube[, varSets[[vs]], drop = FALSE]
  
  # choose species
  y <- NFPD(siteDetections_foliarTraits_BioCube[[species[i]]])
  
  # split into train and test sets
  n <- nrow(spatVars)
  train_idx <- sample(n, floor(0.8*n))
  test_idx  <- setdiff(seq_len(n), train_idx)
  
  X_train <- spatVars[train_idx, ]
  y_train <- y[train_idx]
  X_test  <- spatVars[test_idx, ]
  y_test  <- y[test_idx]
  
  dtrain <- xgb.DMatrix(as.matrix(X_train), label = y_train)
  dtest  <- xgb.DMatrix(as.matrix(X_test),  label = y_test)
  
  # set xgboost parameters
  params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)
  
  # run 5-fold cross validation (find best number of iterations)
  cv <- xgb.cv(params = params, data = dtrain,
               nrounds = 100, nfold = 5, metrics = "rmse",
               verbose = 0, early_stopping_rounds = 10)
  
  # run xgboost
  xgb_model <- xgb.train(params = params, data = dtrain, nrounds = cv$early_stop$best_iteration)
  
  # predict on test data
  pred <- predict(xgb_model, dtest)
  
  # calculate R2 of prediction vs test data
  R2_test  <- 1 - sum((y_test  - pred)^2) / sum((y_test  - mean(y_test))^2)
  RMSE_test <- sqrt(mean((pred - y_test)^2))
  R2_train <- 1 - sum((y_train - predict(xgb_model, dtrain))^2) / sum((y_train - mean(y_train))^2)
  RMSE_train <- sqrt(mean((predict(xgb_model, dtrain) - y_train)^2))
  
  # calculate variable importance
  imp <- xgb.importance(colnames(spatVars), model = xgb_model)
  imp$species   <- species[i]
  imp$varSet    <- vs
  imp$iteration <- iter
  
  # calculate SHAP importance
  shap <- shap.values(xgb_model, X_train = X_test)
  shap_importance <- data.frame(
    Feature = names(shap$mean_shap_score),
    SHAP_importance = shap$mean_shap_score
  )
  
  # keep importance data
  imp <- merge(imp, shap_importance)
  
  # keep model data
  model_row <- data.frame(
    species = species[i],
    varSet = names(varSets)[vs],
    iteration = iter,
    zeros = sum(y == 0, na.rm = TRUE),
    best_nrounds_cv = cv$early_stop$best_iteration,
    R2_test = R2_test,
    RMSE_test = RMSE_test,
    R2_train = R2_train,
    RMSE_train = RMSE_train
  )
  
  list(model = model_row, importance = imp)
}


# run the fxn -------------------------------------- (~30 min for 10 iterations)
handlers(global = TRUE) 

with_progress({
  p <- progressor(along = seq_len(nrow(tasks)))
  
  results <- future_lapply(
    seq_len(nrow(tasks)),
    function(k) {
      res <- run_task(
        vs   = tasks$varSet[k],
        i    = tasks$species_i[k],
        iter = tasks$iter[k]
      )
      p(sprintf("Finished task %d", k))
      res
    }
  )
})

# turn off parallel processing 
plan(sequential)

# combine outputs --------------------------------------------------------------
xgb_modelParameters <- do.call(rbind, lapply(results, `[[`, "model"))
xgb_variableImportance <- do.call(rbind, lapply(results, `[[`, "importance"))

# save it ----------------------------------------------------------------------
write_rds(xgb_modelParameters, "data/xgb_modelParameters_20260325.rds")
write_csv(xgb_modelParameters, "data/xgb_modelParameters_20260325.csv")

write_rds(xgb_variableImportance, "data/xgb_variableImportance_20260325.rds")
write_csv(xgb_variableImportance, "data/xgb_variableImportance_20260325.csv")
#xgb_variableImportance <- readRDS("data/xgb_variableImportance_20260325.rds")













