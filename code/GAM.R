
# comparing model types
# GAM
# 13 Nov 2025

library(mgcv)
library(rsample)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full bird data ----------------------------------------------------------

siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# variable groups --------------------------------------------------------------

# exclude predictors with too few unique values
x <- as.data.frame(sapply(spatVars_clean, function(x) length(unique(x))))
# gsp, GHM_bu, GHM_en, NTL_2020 all have <100 unique values

spatVars_clean  = siteDetections[c(97:110, 112:118, 121:134, 136:182)]

# run full GAM -----------------------------------------------------------------

formula <- as.formula(paste("NFPD(siteDetections$Western.Kingbird) ~", paste(paste0("s(", names(spatVars_clean), ")"), collapse = " + ")))
model <- gam(formula, data = spatVars_clean) # slow

summary(model)

pred <- predict(model)

plot(NFPD(siteDetections$Western.Kingbird), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "GAM: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(NFPD(siteDetections$Western.Kingbird), pred), 2)), bty = "n")



# add term shrinkage and selection ---------------------------------------------

formula <- as.formula(paste(
  "NFPD(Western.Kingbird) ~",
  paste(paste0("s(", names(spatVars_clean), ", bs='ts')"), collapse = " + ")
))

GAM_model <- gam(formula, data = siteDetections, select = TRUE) # slow

summary(GAM_model)

pred <- predict(GAM_model)

plot(NFPD(siteDetections$Western.Kingbird), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "GAM: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(NFPD(siteDetections$Western.Kingbird), pred), 2)), bty = "n")

# the F statistic =~ variable importance
imp <- summary(GAM_model)$s.table[ , "F"]
imp <- as.data.frame(sort(imp, decreasing = TRUE))

gam.check(GAM_model)

# add cross validation

set.seed(123)
folds <- vfold_cv(siteDetections, v=10)

cv_err <- sapply(folds$splits, function(s) {
  train <- analysis(s)
  test  <- assessment(s)
  m <- gam(formula, data=train, select=TRUE)
  mean((predict(m, newdata=test) - test[[all.vars(formula)[1]]])^2)
})
mean(cv_err)

# train/test datasets

set.seed(123)  # for reproducibility
n <- nrow(siteDetections)
train_index <- sample(seq_len(n), size = 0.8 * n)

train_data <- siteDetections[train_index, ]
test_data  <- siteDetections[-train_index, ]

# Fit the GAM on the training set
GAM_model <- gam(formula, data = train_data, select = TRUE) # slow

# Predict on the test set
pred <- predict(GAM_model, newdata = test_data)


plot(NFPD(test_data$Western.Kingbird), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "GAM: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(NFPD(test_data$Western.Kingbird), pred), 2)), bty = "n")

# Evaluate performance
obs <- test_data[[all.vars(formula)[1]]]
RMSE <- sqrt(mean((obs - pred)^2))
R2   <- 1 - sum((obs - pred)^2)/sum((obs - mean(obs))^2)

RMSE
R2
