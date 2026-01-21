
# xgboost
# trying out different model types
# 12 Nov 2025

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

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full bird data ----------------------------------------------------------

siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# remove bad vars
siteDetections <- siteDetections %>% select(-gsp, -GHM_bu, -GHM_en, -NTL_2020)

# variable groups --------------------------------------------------------------

spatVars  = siteDetections[97:178]
notTraits = siteDetections[109:178]
notAnthr  = siteDetections[c(97:108, 121:178)]
notClim   = siteDetections[c(97:120, 148:178)]
notFxn    = siteDetections[c(97:147, 151:174)]
notPheno  = siteDetections[c(97:150, 161:178)]
notTopo   = siteDetections[c(97:160, 175:178)]

traits  = siteDetections[97:108]
anthr = siteDetections[109:120]
climate = siteDetections[121:147]
fxn = siteDetections[c(148:150, 175:178)]
pheno = siteDetections[151:160]
topo = siteDetections[161:174]

# run model --------------------------------------------------------------------

set.seed(123)

# Fit XGBoost regression model with K-fold CV (this gives an ntree_limit warning. b/c caret) (a little slow)
xgb_model <- train(x = spatVars, y = siteDetections$Western.Kingbird,
                   method = "xgbTree",
                   trControl = trainControl(method = "cv", number = 5),
                   tuneLength = 5)

# The resampling results from CV
print(xgb_model)

varImp(xgb_model)

# Predictions on the same data (or a separate test set if desired)
pred <- predict(xgb_model)

plot(siteDetections$Western.Kingbird, pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "xgboost: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(siteDetections$Western.Kingbird, pred), 2)), bty = "n")


# Performance metrics on full data
RMSE <- sqrt(mean((preds - siteDetections$Acorn.Woodpecker)^2))
R2 <- cor(preds, siteDetections$Acorn.Woodpecker)^2

cat("RMSE:", RMSE, "\n")
cat("R2:", R2, "\n")

# run model without caret ------------------------------------------------------

dtrain <- xgb.DMatrix(data = as.matrix(spatVars), label = NFPD(siteDetections$Western.Kingbird))

params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)

cv <- xgb.cv(params = params, # no error message. Quicker.
             data = dtrain,
             nrounds = 100,
             nfold = 5,
             metrics = "rmse",
             verbose = 1,
             early_stopping_rounds = 10)

cv$best_iteration

xgb_model <- xgb.train(params = params,
                         data = dtrain,
                         nrounds = cv$best_iteration)

imp <- xgb.importance(feature_names = colnames(spatVars), model = xgb_model)

ggplot(imp, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col() +
  coord_flip() +
  labs(x = "Feature", y = "Gain", title = "XGBoost Feature Importance (Gain)") +
  theme_minimal()

dtest <- xgb.DMatrix(as.matrix(spatVars))

pred <- predict(xgb_model, dtest)

plot(NFPD(siteDetections$Western.Kingbird), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "xgboost: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(NFPD(siteDetections$Western.Kingbird), pred), 2)), bty = "n")


# SHAP -------------------------------------------------------------------------

shap <- shap.values(xgb_model, X_train = as.matrix(spatVars))
shap_values <- shap$shap_score          # SHAP values for each sample × feature
shap_importance <- shap$mean_shap_score
shap_importance <- data.frame(Feature = names(shap_importance),
                              Importance = shap_importance)
rownames(shap_importance) <- NULL

ggplot(shap_importance, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(x = "Feature", y = "Mean |SHAP|", 
       title = "Global SHAP Feature Importance") +
  theme_minimal()

sum(shap_importance$Importance)

# loop it ----------------------------------------------------------------------

xgb_modelParameters <- data.frame()

for (i in seq_along(species)) {
  
  set.seed(123)  # for reproducibility
  
  # Observed response
  y <- NFPD(siteDetections[[names(species[i])]])
  
  # Row indices
  n <- nrow(spatVars)
  train_idx <- sample(n, size = floor(0.8*n))
  test_idx  <- setdiff(seq_len(n), train_idx)
  
  # Split predictors and response
  X_train <- spatVars[train_idx, ]
  y_train <- y[train_idx]
  X_test  <- spatVars[test_idx, ]
  y_test  <- y[test_idx]
  
  # Create DMatrix objects
  dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
  dtest  <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
  
  # XGBoost parameters
  params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)
  
  # Cross-validation on training set
  cv <- xgb.cv(params = params,
               data = dtrain,
               nrounds = 100,
               nfold = 5,
               metrics = "rmse",
               verbose = 1,
               early_stopping_rounds = 10)
  
  # Train final model on training data
  xgb_model <- xgb.train(params = params,
                         data = dtrain,
                         nrounds = cv$best_iteration)
  
  # Predict on test set
  pred <- predict(xgb_model, dtest)
  
  # Compute R^2 on test set
  R2_test <- 1 - sum((y_test - pred)^2) / sum((y_test - mean(y_test))^2)
  R2_train <- 1 - sum((y_train - predict(xgb_model, dtrain))^2) / sum((y_train - mean(y_train))^2)
  
  
  xgb_modelParameters.tmp <- data.frame(species = names(species[i]),
                                        varSet = "spatVars",
                                        zeros = sum(species[i] == 0),
                                        best_iteration_cv = cv$best_iteration,
                                        R2_test = R2_test,
                                        R2_train = R2_train)
  
  xgb_modelParameters <- rbind(xgb_modelParameters, xgb_modelParameters.tmp)
}

ggplot(xgb_modelParameters, aes(y=reorder(species, R2_test), x = R2_test)) +
  geom_point() +
  theme_minimal() +
  ggtitle("XGBoost model R2 (test data)")

ggplot(xgb_modelParameters, aes(y=reorder(species, R2_train), x = R2_train)) +
  geom_point() +
  theme_minimal() +
  ggtitle("XGBoost model R2 (train data)")

ggplot(xgb_modelParameters, aes(y=reorder(species, R2_test))) +
  geom_point(aes(x = R2_train), color="red") +
  geom_point(aes(x = R2_test)) +
  theme_minimal() +
  ggtitle("XGBoost model R2")



# double loop ------------------------------------------------------------------

varSets <- list(
  spatVars  = siteDetections[97:178],
  notTraits = siteDetections[109:178],
  notAnthr  = siteDetections[c(97:108, 121:178)],
  notClim   = siteDetections[c(97:120, 148:178)],
  notFxn    = siteDetections[c(97:147, 151:174)],
  notPheno  = siteDetections[c(97:150, 161:178)],
  notTopo   = siteDetections[c(97:160, 175:178)],
  
  traits  = siteDetections[97:108],
  anthr = siteDetections[109:120],
  climate = siteDetections[121:147],
  fxn = siteDetections[c(148:150, 175:178)],
  pheno = siteDetections[151:160],
  topo = siteDetections[161:174]
)

# start loop

xgb_modelParameters <- data.frame()
xgb_variableImportance <- data.frame()

for(vs in names(varSets)) {
  
  spatVars <- varSets[[vs]]
  
  for (i in seq_along(species)) {
    
    set.seed(123)  # for reproducibility
    
    # Observed response
    y <- NFPD(siteDetections[[names(species[i])]])
    
    # Row indices
    n <- nrow(spatVars)
    train_idx <- sample(n, size = floor(0.8*n))
    test_idx  <- setdiff(seq_len(n), train_idx)
    
    # Split predictors and response
    X_train <- spatVars[train_idx, ]
    y_train <- y[train_idx]
    X_test  <- spatVars[test_idx, ]
    y_test  <- y[test_idx]
    
    # Create DMatrix objects
    dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
    dtest  <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
    
    # XGBoost parameters
    params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)
    
    # Cross-validation on training set
    cv <- xgb.cv(params = params,
                 data = dtrain,
                 nrounds = 100,
                 nfold = 5,
                 metrics = "rmse",
                 verbose = 1,
                 early_stopping_rounds = 10)
    
    # Train final model on training data
    xgb_model <- xgb.train(params = params,
                           data = dtrain,
                           nrounds = cv$best_iteration)
    
    # Predict on test set
    pred <- predict(xgb_model, dtest)
    
    # Compute R^2 on test set
    R2_test <- 1 - sum((y_test - pred)^2) / sum((y_test - mean(y_test))^2)
    R2_train <- 1 - sum((y_train - predict(xgb_model, dtrain))^2) / sum((y_train - mean(y_train))^2)
    
    # variable importance
    imp <- xgb.importance(feature_names = colnames(spatVars), model = xgb_model)
    imp$species <- names(species[i])
    imp$varSet <- vs
    
    #save results
    xgb_modelParameters.tmp <- data.frame(species = names(species[i]),
                                          varSet = vs,
                                          zeros = sum(species[i] == 0),
                                          best_iteration_cv = cv$best_iteration,
                                          R2_test = R2_test,
                                          R2_train = R2_train)
    
    xgb_modelParameters <- rbind(xgb_modelParameters, xgb_modelParameters.tmp)
    
    xgb_variableImportance <- rbind(xgb_variableImportance, imp)
  }
}
  
  
ggplot(xgb_modelParameters, aes(y=reorder(species, R2_test), x = R2_test, color=varSet)) +
  geom_point() +
  theme_minimal() +
  ggtitle("XGBoost model R2 (test data)")


# which variables are most important for which species? ------------------------

species_list <- species

i=12

temp <- subset(xgb_variableImportance, species== names(species_list)[i] & xgb_variableImportance$varSet=="spatVars")

temp$label <- ifelse(temp$Gain > 0.05, temp$Feature, "")

temp <- merge(temp, tidy, by.x= "Feature", by.y= "VariableName")

temp <- temp[order(temp$Category, -temp$Gain), ]

temp$Feature <- factor(temp$Feature, levels = temp$Feature)
temp$label <- factor(temp$label, levels = temp$Feature)

r2 <- round(xgb_modelParameters$R2_test[xgb_modelParameters$varSet == "spatVars" & xgb_modelParameters$species == names(species_list)[i]], 2)

ggplot(temp, aes(x = "", y = Gain, fill = Feature)) +
  geom_col(width = 1) +
  scale_fill_discrete(type = "gradient") +
  new_scale_fill() +
  geom_col(aes(x=1.6, fill = Category), width = 0.1) +
  scale_fill_manual(values = c("red", "green", "lightblue", "blue", "purple", "pink")) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), color="white") +
  coord_polar(theta = "y") +
  theme_void() +
  ggtitle(paste(temp$species[1], ", R2 =", r2))
  

# which categories of variables tend to be most important?

temp <- subset(xgb_variableImportance, xgb_variableImportance$varSet=="spatVars")
temp_p <- subset(xgb_modelParameters, xgb_modelParameters$varSet=="spatVars")[c("species", "R2_test")]

temp <- merge(temp, tidy, by.x= "Feature", by.y= "VariableName")
temp <- merge(temp, temp_p, by= "species")

temp <- subset(temp, temp$R2_test > 0.45)

temp <- temp %>% group_by(species, Category) %>% 
  summarise(cumGain = sum(Gain))

ggplot(temp, aes(x=cumGain, y=Category, fill=Category)) +
  geom_violin() +
  geom_point() +
  geom_vline(xintercept = 0.166) +
  theme_minimal() +
  xlab("Model Weight (%)")


ggplot(siteDetections, aes(x=LMA, y=Yellow.rumped.Warbler)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()


# which species models significantly improve with traits? ----------------------

temp <- subset(xgb_modelParameters, xgb_modelParameters$varSet=="spatVars" | xgb_modelParameters$varSet=="notClim")

ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet)) +
  geom_point() +
  theme_minimal()

# R2 increase in accuracy

temp_1 <- subset(xgb_modelParameters, xgb_modelParameters$varSet=="spatVars")
temp_1$svR2 <- temp_1$R2_test
temp_1 <- temp_1[c("species", "svR2")]

temp_2 <- subset(xgb_modelParameters, 
                 xgb_modelParameters$varSet=="notTraits" | 
                   xgb_modelParameters$varSet=="notClim" |
                   xgb_modelParameters$varSet=="notFxn" |
                   xgb_modelParameters$varSet=="notPheno" |
                   xgb_modelParameters$varSet=="notTopo" |
                   xgb_modelParameters$varSet=="notAnthr")
temp_2$subsetR2 <- temp_2$R2_test
temp_2 <- temp_2[c("species", "varSet", "subsetR2")]

temp <- merge(temp_1, temp_2)
temp$R2improvement <- temp$svR2 - temp$subsetR2

ggplot(temp, aes(x=R2improvement, y = varSet, fill=varSet)) +
  geom_violin() +
  theme_minimal()



# double loop with SHAP and iterations ---------------------------------(this has been running for hours. must parallelize)(this is the output I'm using now)

varSets <- list(
  spatVars  = siteDetections[97:178],
  notTraits = siteDetections[109:178],
  notAnthr  = siteDetections[c(97:108, 121:178)],
  notClim   = siteDetections[c(97:120, 148:178)],
  notFxn    = siteDetections[c(97:147, 151:174)],
  notPheno  = siteDetections[c(97:150, 161:178)],
  notTopo   = siteDetections[c(97:160, 175:178)],
  
  traits  = siteDetections[97:108],
  anthr = siteDetections[109:120],
  climate = siteDetections[121:147],
  fxn = siteDetections[c(148:150, 175:178)],
  pheno = siteDetections[151:160],
  topo = siteDetections[161:174]
)

# start loop

xgb_modelParameters <- data.frame()
xgb_variableImportance <- data.frame()

for(vs in names(varSets)) {
  
  spatVars <- varSets[[vs]]
  
  for (i in seq_along(species)) {
    
    for (iter in 1:10) {
      
      set.seed(123 + iter)
      
      # Observed response
      y <- NFPD(siteDetections[[names(species[i])]])
      
      # Row indices
      n <- nrow(spatVars)
      train_idx <- sample(n, size = floor(0.8*n))
      test_idx  <- setdiff(seq_len(n), train_idx)
      
      # Split predictors and response
      X_train <- spatVars[train_idx, ]
      y_train <- y[train_idx]
      X_test  <- spatVars[test_idx, ]
      y_test  <- y[test_idx]
      
      # Create DMatrix objects
      dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
      dtest  <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
      
      # XGBoost parameters
      params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)
      
      # Cross-validation on training set
      cv <- xgb.cv(params = params,
                   data = dtrain,
                   nrounds = 100,
                   nfold = 5,
                   metrics = "rmse",
                   verbose = 1,
                   early_stopping_rounds = 10)
      
      # Train final model on training data
      xgb_model <- xgb.train(params = params,
                             data = dtrain,
                             nrounds = cv$best_iteration)
      
      # Predict on test set
      pred <- predict(xgb_model, dtest)
      
      # Compute R^2 on test set
      R2_test <- 1 - sum((y_test - pred)^2) / sum((y_test - mean(y_test))^2)
      R2_train <- 1 - sum((y_train - predict(xgb_model, dtrain))^2) / sum((y_train - mean(y_train))^2)
      
      # variable importance
      imp <- xgb.importance(feature_names = colnames(spatVars), model = xgb_model)
      imp$species <- names(species[i])
      imp$varSet <- vs
      
      # SHAP
      shap <- shap.values(xgb_model, X_train = dtest)
      shap_values <- shap$shap_score
      shap_importance <- shap$mean_shap_score
      shap_importance <- data.frame(Feature = names(shap_importance),
                                    SHAP_importance = shap_importance)
      rownames(shap_importance) <- NULL
      
      imp <- merge(imp, shap_importance)
      
      imp$iteration <- iter
      
      #save results
      xgb_modelParameters.tmp <- data.frame(species = names(species[i]),
                                            varSet = vs,
                                            zeros = sum(species[i] == 0),
                                            iteration = iter,
                                            best_nrounds_cv = cv$best_iteration,
                                            R2_test = R2_test,
                                            R2_train = R2_train)
      
      xgb_modelParameters <- rbind(xgb_modelParameters, xgb_modelParameters.tmp)
      
      xgb_variableImportance <- rbind(xgb_variableImportance, imp)
      
    }
  }
}


# save it
write.csv(xgb_modelParameters, "xgb_modelParameters_20251124.csv")
write.csv(xgb_variableImportance, "xgb_variableImportance_20251124.csv")

xgb_modelParameters <- read.csv("data/xgb_modelParameters_20251124.csv")
xgb_variableImportance <- read.csv("data/xgb_variableImportance_20251124.csv")


# plot it
temp <- subset(xgb_modelParameters, xgb_modelParameters$species=="Red.breasted.Nuthatch")

temp$varSet <- factor(temp$varSet, levels = c("spatVars","notTraits","notAnthr","notClim",
                                              "notFxn","notPheno","notTopo",
                                              "traits","anthr","climate","fxn","pheno","topo"))

ggplot(temp, aes(y=R2_test, x=varSet)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "spatVars") +
  coord_flip() +
  theme_minimal()

temp <- subset(xgb_variableImportance, xgb_variableImportance$varSet=="spatVars" & xgb_variableImportance$species=="California.Scrub.Jay")

temp <- merge(temp, tidy, by.x= "Feature", by.y= "VariableName")

ggplot(temp, aes(x=Gain, y=reorder(Feature, Gain), color=Category)) +
  geom_boxplot() +
  scale_color_manual(values = c("#70A4AF", "#7C4584", "#FDC71B", "#D97C55", "#842B3B", "#A8BBA3")) +
  theme_minimal()

ggplot(temp, aes(x=SHAP_importance, y=reorder(Feature, SHAP_importance), color=Category, fill=Category)) +
  geom_boxplot() +
  scale_color_manual(values = c("#70A4AF", "#7C4584", "#FDC71B", "#D97C55", "#842B3B", "#A8BBA3")) +
  scale_fill_manual(values = c("#A4C4CB", "#AC86B0", "#FFDB7D", "#ECAB90", "#B4777C", "#C7D3C4")) +
  theme_minimal() +
  ylab("Variable") +
  xlab("SHAP Importance")


top10 <- temp |> 
  dplyr::group_by(Feature) |> 
  dplyr::summarise(mean_shap = mean(SHAP_importance, na.rm=TRUE)) |> 
  dplyr::slice_max(mean_shap, n = 10)

temp <- temp |> dplyr::filter(Feature %in% top10$Feature)

my_cols  <- c("Climate"="#70A4AF", "Disturbance"="#7C4584", "Function"="#FDC71B", "Phenology"="#D97C55", "Physiology"="#842B3B", "Terrain"="#A8BBA3")

my_fills <- c("Climate"="#A4C4CB", "Disturbance"="#AC86B0", "Function"="#FFDB7D", "Phenology"="#ECAB90", "Physiology"="#B4777C", "Terrain"="#C7D3C4")

ggplot(temp, aes(x=SHAP_importance, y=reorder(Feature, SHAP_importance), color=Category, fill=Category)) +
  geom_boxplot() +
  scale_color_manual(values = my_cols) +
  scale_fill_manual(values = my_fills) +
  theme_minimal() +
  ylab("Variable") +
  xlab("SHAP Importance")

# which species R2 significantly improves with each category? ------------------

cat_ttest <- xgb_modelParameters %>%
  group_by(species) %>%
  summarise(notTraits = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notTraits"], alternative="greater")$p.value,
            notClim = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notClim"], alternative="greater")$p.value,
            notFxn = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notFxn"], alternative="greater")$p.value,
            notAnthr = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notAnthr"], alternative="greater")$p.value,
            notPheno = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notPheno"], alternative="greater")$p.value,
            notTopo = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notTopo"], alternative="greater")$p.value,
            R2 = max(R2_test))

# plot it

ssp <- cat_ttest$species[cat_ttest$notTraits < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notTraits")
temp$varSet[temp$varSet=="notTraits"] <- "Without Leaf Physiology"
temp$varSet[temp$varSet=="spatVars"] <- "With Leaf Physiology"

p1 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#842B3B", "black")) +
  scale_fill_manual(values = c("#aa384c", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  xlab("Best Model R2") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notClim < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notClim")
temp$varSet[temp$varSet=="notClim"] <- "Without Climate"
temp$varSet[temp$varSet=="spatVars"] <- "With Climate"

p2 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#70A4AF", "black")) +
  scale_fill_manual(values = c("#91b9c1", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  xlab("Best Model R2") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notFxn < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notFxn")
temp$varSet[temp$varSet=="notFxn"] <- "Without Ecosystem Function"
temp$varSet[temp$varSet=="spatVars"] <- "With Ecosystem Function"

p3 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#A8BBA3", "black")) +
  scale_fill_manual(values = c("#c4d1c0", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  xlab("Best Model R2") +
  ylab("")

p3 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#FDC71B", "black")) +
  scale_fill_manual(values = c("#FFE090", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  xlab("Best Model R2") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notPheno < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notPheno")
temp$varSet[temp$varSet=="notPheno"] <- "Without Phenology"
temp$varSet[temp$varSet=="spatVars"] <- "With Phenology"

p4 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#D97C55", "black")) +
  scale_fill_manual(values = c("#e29c7f", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  xlab("Best Model R2") +
  ylab("")


annotate_figure(
  ggarrange(p3, p1, p2, p4, ncol = 1, nrow = 4, heights = c(4, 3, 2, 1), align="hv"),
  left = text_grob("Species with significant improvement", rot = 90, size = 14, vjust = 1))

ggsave("SpeciesBestR2_20251202.pdf", height=10, width=8)


# SHAP scores vs variables -----------------------------------------------------

# Run one model
set.seed(123)

# Observed response
y <- NFPD(siteDetections$Fox.Sparrow)

# Row indices
spatVars  <- siteDetections[97:178]
n <- nrow(spatVars)
train_idx <- sample(n, size = floor(0.8*n))
test_idx  <- setdiff(seq_len(n), train_idx)

# Split predictors and response
X_train <- spatVars[train_idx, ]
y_train <- y[train_idx]
X_test  <- spatVars[test_idx, ]
y_test  <- y[test_idx]

# Create DMatrix objects
dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
dtest  <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
dfull  <- xgb.DMatrix(data = as.matrix(spatVars), label = y)

# XGBoost parameters
params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)

# Cross-validation on training set
cv <- xgb.cv(params = params,
             data = dtrain,
             nrounds = 100,
             nfold = 5,
             metrics = "rmse",
             verbose = 1,
             early_stopping_rounds = 10)

# Train final model on training data
xgb_model <- xgb.train(params = params,
                       data = dtrain,
                       nrounds = cv$best_iteration)

# Predict on test set
pred <- predict(xgb_model, dtest)

# Compute R^2 on test set
R2_test <- 1 - sum((y_test - pred)^2) / sum((y_test - mean(y_test))^2)
R2_train <- 1 - sum((y_train - predict(xgb_model, dtrain))^2) / sum((y_train - mean(y_train))^2)

# SHAP test set
shap <- shap.values(xgb_model, X_train = dtest)
shap_values <- shap$shap_score
temp <- data.frame(species ="Fox.Sparrow",  shap_CanopyHeight = shap_values$CanopyHeight, CanopyHeight = X_test$CanopyHeight,
                   shap_LMA = shap_values$LMA, LMA = X_test$LMA,
                   shap_Phenolics = shap_values$Phenolics, Phenolics = X_test$Phenolics,
                   shap_Potassium = shap_values$Potassium, Potassium = X_test$Potassium)

# SHAP full set
shap <- shap.values(xgb_model, X_train = dfull)
shap_values <- shap$shap_score
temp <- data.frame(species ="Fox.Sparrow",  shap_CanopyHeight = shap_values$CanopyHeight, CanopyHeight = spatVars$CanopyHeight,
                   shap_LMA = shap_values$LMA, LMA = spatVars$LMA,
                   shap_Phenolics = shap_values$Phenolics, Phenolics = spatVars$Phenolics,
                   shap_Potassium = shap_values$Potassium, Potassium = spatVars$Potassium)


# canopy height
ggplot() +
  geom_point(data = siteDetections, aes(x=CanopyHeight/100, y=NFPD(Fox.Sparrow)), color="#FDC71B") +
  geom_smooth(data = subset(temp, temp$species=="Fox.Sparrow"), aes(x=CanopyHeight/100, y=((shap_CanopyHeight *40) +0.4)), color="black") + 
  scale_y_continuous(name="Fox.Sparrow \nnormalized detection rate", sec.axis = sec_axis(~ . /40 -0.4/40, name="SHAP value")) +
  theme_minimal() +
  xlab("Canopy Height") +
  xlim(c(0,80)) +
  theme(axis.title.y.left=element_text(color="#FDC71B"))
ggsave("figures/FSpar_CH.PDF", width=6, height=3)

# LMA
ggplot() +
  geom_point(data = siteDetections, aes(x=LMA, y=NFPD(Fox.Sparrow)), color="#842B3B") +
  #geom_smooth(data = siteDetections, method = "lm", aes(x=LMA, y=NFPD(Fox.Sparrow)), color="black") +
  geom_smooth(data = subset(temp, temp$species=="Fox.Sparrow"), aes(x=LMA, y=((shap_LMA *10) +0.4)), color="black") + 
  scale_y_continuous(name="Fox.Sparrow \nnormalized detection rate", sec.axis = sec_axis(~ . /10 -0.4/10, name="SHAP value")) +
  theme_minimal() +
  xlab("LMA (g/m2)") +
  xlim(c(0,400)) +
  theme(axis.title.y.left=element_text(color="#842B3B"))
ggsave("figures/FSpar_LMA.PDF", width=6, height=3)

# phenolics
ggplot() +
  geom_point(data = siteDetections, aes(x=Phenolics, y=NFPD(Fox.Sparrow)), color="#842B3B") +
  geom_smooth(data = subset(temp, temp$species=="Fox.Sparrow"), aes(x=Phenolics, y=((shap_Phenolics *5) +0.4)), color="black") + 
  scale_y_continuous(name="Fox.Sparrow \nnormalized detection rate", sec.axis = sec_axis(~ . /5 -0.4/5, name="SHAP value")) +
  theme_minimal() +
  xlab("Phenolics") +
  xlim(c(30,120)) +
  theme(axis.title.y.left=element_text(color="#842B3B"))
ggsave("figures/FSpar_Phenolics.PDF", width=6, height=3)

# potassium
ggplot() +
  geom_point(data = siteDetections, aes(x=Potassium, y=NFPD(Fox.Sparrow)), color="#842B3B") +
  geom_smooth(data = subset(temp, temp$species=="Fox.Sparrow"), aes(x=Potassium, y=((shap_Potassium *10) +0.5)), color="black") + 
  scale_y_continuous(name="Fox.Sparrow \nnormalized detection rate", sec.axis = sec_axis(~ . /10 -0.5/10, name="SHAP value")) +
  theme_minimal() +
  xlab("Potassium") +
  xlim(c(0,15)) +
  theme(axis.title.y.left=element_text(color="#842B3B"))
ggsave("figures/FSpar_Potassium.PDF", width=6, height=3)

################################################################################

# double loop with SHAP and iterations in parallel ---------------(as yet unrun)

# Use all cores
plan(multisession)   # or multicore on Linux/macOS

tasks <- expand.grid(
  varSet = names(varSets),
  species_i = seq_along(species),
  iter = 1:10
)

# make a fxn
run_task <- function(vs, i, iter) {
  set.seed(123 + iter)
  
  spatVars <- varSets[[vs]]
  y <- NFPD(siteDetections[[names(species[i])]])
  
  n <- nrow(spatVars)
  train_idx <- sample(n, floor(0.8*n))
  test_idx  <- setdiff(seq_len(n), train_idx)
  
  X_train <- spatVars[train_idx, ]
  y_train <- y[train_idx]
  X_test  <- spatVars[test_idx, ]
  y_test  <- y[test_idx]
  
  dtrain <- xgb.DMatrix(as.matrix(X_train), label = y_train)
  dtest  <- xgb.DMatrix(as.matrix(X_test),  label = y_test)
  
  params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)
  
  cv <- xgb.cv(params = params, data = dtrain,
               nrounds = 100, nfold = 5, metrics = "rmse",
               verbose = 0, early_stopping_rounds = 10)
  
  xgb_model <- xgb.train(params = params, data = dtrain, nrounds = cv$best_iteration)
  
  pred <- predict(xgb_model, dtest)
  R2_test  <- 1 - sum((y_test  - pred)^2) / sum((y_test  - mean(y_test))^2)
  R2_train <- 1 - sum((y_train - predict(xgb_model, dtrain))^2) / sum((y_train - mean(y_train))^2)
  
  imp <- xgb.importance(colnames(spatVars), model = xgb_model)
  imp$species   <- names(species[i])
  imp$varSet    <- vs
  imp$iteration <- iter
  
  shap <- shap.values(xgb_model, X_train = X_test)
  shap_importance <- data.frame(
    Feature = names(shap$mean_shap_score),
    SHAP_importance = shap$mean_shap_score
  )
  
  imp <- merge(imp, shap_importance)
  
  model_row <- data.frame(
    species = names(species[i]),
    varSet = vs,
    iteration = iter,
    zeros = sum(species[i] == 0),
    best_nrounds_cv = cv$best_iteration,
    R2_test = R2_test,
    R2_train = R2_train
  )
  
  list(model = model_row, importance = imp)
}


# run the fxn
results <- future_lapply(
  seq_len(nrow(tasks)),
  function(k) {
    run_task(
      vs  = tasks$varSet[k],
      i   = tasks$species_i[k],
      iter = tasks$iter[k]
    )
  }
)


# combine outputs
xgb_modelParameters <- do.call(rbind, lapply(results, `[[`, "model"))
xgb_variableImportance <- do.call(rbind, lapply(results, `[[`, "importance"))





