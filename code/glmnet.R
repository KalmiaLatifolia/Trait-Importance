
# comparing model types
# ridge, lasso, elastic net models
# 4 Nov 2025

library(glmnet)
library(mgcv)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full bird data ----------------------------------------------------------

siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# variable sets ----------------------------------------------------------------

spatVars_clean  = siteDetections[c(97:134, 136:182)]

x <- as.matrix(spatVars_clean)
y <- NFPD(siteDetections$Western.Kingbird)

# ridge model ------------------------------------------------------------------

# Ridge regression (alpha = 0)
ridge_model <- cv.glmnet(x, y, alpha = 0)

# Best lambda
best_lambda <- ridge_model$lambda.min

# Coefficients at best lambda
coef(ridge_model, s = "lambda.min")

# Predicted values
pred <- predict(ridge_model, x, s = "lambda.min")

# plot
plot(y, pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Ridge: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(y, pred), 2)), bty = "n")


# elastic net ------------------------------------------------------------------

# alpha = 0.5 is a common elastic net choice
enet_model <- cv.glmnet(x, y, alpha = 0.5)   # cv picks optimal lambda

coef(enet_model, s = "lambda.min")

pred <- predict(enet_model, newx = x, s = "lambda.min")

plot(y, pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Elastic Net: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(y, pred), 2)), bty = "n")


# lasso ------------------------------------------------------------------------

lasso_model <- cv.glmnet(x, y, alpha = 1)   # alpha = 1 → LASSO

coef(lasso_model, s = "lambda.min")

pred <- predict(lasso_model, newx = x, s = "lambda.min")

plot(y, pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Lasso: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(y, pred), 2)), bty = "n")


# single species, single variable ----------------------------------------------

m <- gam(NFPD(siteDetections$Acorn.Woodpecker) ~ s(siteDetections$Nitrogen))
m <- glm(NFPD(siteDetections$Acorn.Woodpecker) ~ siteDetections$Nitrogen)

summary(m)
pred <- predict(m)
plot(NFPD(siteDetections$Acorn.Woodpecker), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Observed vs Predicted") +
  abline(0, 1, lty = 2)
cor(NFPD(siteDetections$Acorn.Woodpecker), pred)



