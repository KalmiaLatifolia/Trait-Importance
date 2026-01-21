
# random forest
# comparing model predictions of different model types
# 4 Nov 2025

library(VSURF)
library(randomForest)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full bird data ----------------------------------------------------------

siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# run variable selection -------------------------------------------------------
spatVars_clean  = siteDetections[c(97:134, 136:182)]

vsurf_result <- VSURF(spatVars_clean, NFPD(siteDetections$Western.Kingbird))
selected_vars <- vsurf_result$varselect.pred

# VSURF random forest model ----------------------------------------------------

rf_model <- randomForest(spatVars_clean[, selected_vars], NFPD(siteDetections$Western.Kingbird), importance = TRUE)
rf_model

# predict ----------------------------------------------------------------------

pred <- predict(rf_model)

plot(NFPD(siteDetections$Western.Kingbird), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Random Forest: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(NFPD(siteDetections$Western.Kingbird), pred), 2)), bty = "n")


# what if we use NFPD? ---------------------------------------------------------

FPD <- function(species_col) {species_col / max(species_col)}
meanFPD <- function(species_col) {mean(FPD(species_col)[FPD(species_col)>0])}
NFPD <- function(species_col) {(FPD(species_col))^(log(0.5)/log(meanFPD(species_col)))}

FRichness <- function(data, species_cols, ID_col){
  data.frame(ID = data[[ID_col]],
             FRich = rowSums(sapply(data[species_cols], NFPD), na.rm=TRUE))
}

#
hist(siteDetections$Acorn.Woodpecker)
hist(NFPD(siteDetections$Acorn.Woodpecker))
mean(NFPD(siteDetections$Acorn.Woodpecker))
rf_model <- randomForest(siteDetections[, 97:182], siteDetections$Acorn.Woodpecker, importance = TRUE) # 25%
rf_model <- randomForest(siteDetections[, 97:182], NFPD(siteDetections$Acorn.Woodpecker), importance = TRUE) #64%

pred <- predict(rf_model)

plot(NFPD(siteDetections$Acorn.Woodpecker), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Random Forest: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(NFPD(siteDetections$Acorn.Woodpecker), pred), 2)), bty = "n")


# NFPD, one species, one var ---------------------------------------------------

set.seed(123)
rf_model <- randomForest(NFPD(Acorn.Woodpecker) ~ Nitrogen, data=siteDetections)
rf_model
pred <- predict(rf_model)
plot(NFPD(siteDetections$Acorn.Woodpecker), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Observed vs Predicted") +
  abline(0, 1, lty = 2)
cor(NFPD(siteDetections$Acorn.Woodpecker), pred)
