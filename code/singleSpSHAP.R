
# Single Species SHAP values
# section 17 in full protocol
# Laura Berman
# 9 April 2026

library(readxl)
library(xgboost)
library(SHAPforxgboost)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full data ---------------------------------------------------------------

siteDetections_foliarTraits_BioCube <- readRDS("data/siteDetections_foliarTraits_BioCube_20260320.rds")


# prepare data -----------------------------------------------------------------

# load tidy names 
tidy <- read_excel("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TableS1_Biocube_var_description.xlsx",
                   range = cell_cols(1:4))

# NFPD function 
FPD <- function(species_col) {species_col / max(species_col)}
meanFPD <- function(species_col) {mean(FPD(species_col)[FPD(species_col)>0])}
NFPD <- function(species_col) {(FPD(species_col))^(log(0.5)/log(meanFPD(species_col)))}

# name variable sets
varSets <- list(
  spatVars  = tidy$Variable,
  notTraits = tidy$Variable[tidy$Category != "Traits"],
  notDist   = tidy$Variable[tidy$Category != "Disturbance"],
  notClim   = tidy$Variable[tidy$Category != "Climate"],
  notStr    = tidy$Variable[tidy$Category != "Structure"],
  notPheno  = tidy$Variable[tidy$Category != "Phenology"],
  notTerr   = tidy$Variable[tidy$Category != "Terrain"]
)
varSets <- lapply(varSets, function(x) x[x %in% colnames(siteDetections_foliarTraits_BioCube)])


# get SHAP scores for 1 species ------------------------------------------------

# choose species (25, 70, 32, 29, 31, 24, 61, 48)
species <- colnames(siteDetections_foliarTraits_BioCube)[4:97]
species
i <- 70
y <- NFPD(siteDetections_foliarTraits_BioCube[[species[i]]])

# use full variable set
vs <- 1
spatVars <- siteDetections_foliarTraits_BioCube[, varSets[[vs]], drop = FALSE]

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

# calculate SHAP scores
shap <- shap.values(xgb_model, X_train = X_test)
shap_values <- shap$shap_score

# build dataframe
temp <- data.frame(species = rep(species[i], nrow(shap_values)),  
                   shap_CanopyHeight = shap_values$CA_Function_Favrichon_Sierra_CanopyHeight_30m_20250508, 
                   CanopyHeight = X_test$CA_Function_Favrichon_Sierra_CanopyHeight_30m_20250508,
                   shap_LMA = shap_values$LMA, LMA = X_test$LMA,
                   shap_Phenolics = shap_values$Phenolics, Phenolics = X_test$Phenolics,
                   shap_Potassium = shap_values$Potassium, Potassium = X_test$Potassium,
                   shap_Nitrogen = shap_values$Nitrogen, Nitrogen = X_test$Nitrogen)


# plot it ----------------------------------------------------------------------

# canopy height
ggplot() +
  geom_point(data = siteDetections_foliarTraits_BioCube, aes(x=CA_Function_Favrichon_Sierra_CanopyHeight_30m_20250508/100, y=NFPD(get(species[i]))), color="#FDC71B") +
  geom_smooth(data = temp, aes(x=CanopyHeight/100, y=((shap_CanopyHeight *4) +0.4)), color="black") + 
  scale_y_continuous(name=paste(species[i], "\nnormalized detection rate"), sec.axis = sec_axis(~ . /4 -0.4/4, name="SHAP value")) +
  theme_minimal() +
  xlab("Canopy Height") +
  xlim(c(0,80)) +
  theme(axis.title.y.left=element_text(color="#FDC71B"))
ggsave("figures/LGold_CH.PDF", width=6, height=3)

# LMA
ggplot() +
  geom_point(data = siteDetections_foliarTraits_BioCube, aes(x=LMA, y=NFPD(get(species[i]))), color="#842B3B") +
  geom_smooth(data = temp, aes(x=LMA, y=((shap_LMA *60) +0.4)), color="black") + 
  scale_y_continuous(name=paste(species[i], "\nnormalized detection rate"), sec.axis = sec_axis(~ . /60 -0.4/60, name="SHAP value")) +
  theme_minimal() +
  xlab("LMA (g/m2)") +
  xlim(c(0,400)) +
  theme(axis.title.y.left=element_text(color="#842B3B"))
ggsave("figures/GCKing_LMA.PDF", width=6, height=3)

# Nitrogen
ggplot() +
  geom_point(data = siteDetections_foliarTraits_BioCube, aes(x=Nitrogen, y=NFPD(get(species[i]))), color="#842B3B") +
  geom_smooth(data = temp, aes(x=Nitrogen, y=((shap_Nitrogen *5) +0.4)), color="black") + 
  scale_y_continuous(name=paste(species[i], "\nnormalized detection rate"), sec.axis = sec_axis(~ . /5 -0.4/5, name="SHAP value")) +
  theme_minimal() +
  xlab("Nitrogen") +
  xlim(c(10,30)) +
  theme(axis.title.y.left=element_text(color="#842B3B"))
ggsave("figures/GCKing_Nitrogen.PDF", width=6, height=3)

