
# trying PLSR models
# species relative abundance
# 4 Nov 2025

# https://www.statology.org/partial-least-squares-in-r/

library(pls)
library(plsVarSel)
library(ggplot2)
library(dplyr)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full bird data ----------------------------------------------------------

siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# exclude predictors with too few unique values
x <- as.data.frame(sapply(spatVars_clean, function(x) length(unique(x))))
# gsp, GHM_bu, GHM_en, NTL_2020 all have <100 unique values

# remove bad vars
siteDetections <- siteDetections %>% select(-gsp, -GHM_bu, -GHM_en, -NTL_2020)

# run plsr ---------------------------------------------------------------------

set.seed(123)

model <- plsr(siteDetections$Acorn.Woodpecker ~ ., data = siteDetections[, c(97:178)], scale = TRUE, validation = "CV")

summary(model)

# 3 comps has the lowest cross-validated RMSEP
# Cumulative variance explained at 3 components is 33.11%

validationplot(model, val.type = "RMSEP")   # choose optimal number of components (pick the smallest)
validationplot(model, val.type = "R2")
coef(model, ncomp=6)                        # regression coefficients = strength and direction of each var. Important to specify correct number of components (not overfit)
loadings(model)[,1:3]                       # how much each variable contributes to components. strength and direction
scores(model)[,1:3]                         # transformed component scores. where each obs. falls in each comp. 
VIP(model, opt.comp = 3)                    # variable importance (cumulative across 3 components)

# variable importance ----------------------------------------------------------

variableImportance <- stack(VIP(model, opt.comp = 3))
colnames(variableImportance) <- c("variableImportance", "variable")

tidy <- read.csv("data/tidyNames.csv")
tidy$FileName <- NULL
variableImportance <- merge(variableImportance, tidy, by.x= "variable", by.y= "VariableName")

ggplot(variableImportance, aes(x=reorder(variable, variableImportance), y=variableImportance, fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_col() + coord_flip() +
  xlab("Habitat Variable") +
  ylab("Variable Importance") +
  theme_minimal()
                                  
sum(variableImportance$variableImportance)/ nrow(variableImportance)

# predict ----------------------------------------------------------------------

pred <- predict(model, ncomp = 3)

plot(siteDetections$Acorn.Woodpecker, pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "PLSR: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(siteDetections$Acorn.Woodpecker, pred), 2)), bty = "n")

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

# try other variable groups ----------------------------------------------------

set.seed(123)

model <- plsr(siteDetections$Acorn.Woodpecker ~ ., data = notTopo, scale = TRUE, validation = "CV")

which.min(RMSEP(model, estimate = "adjCV")$val) -1
validationplot(model, val.type = "RMSEP")
summary(model)

# variable importance by species -----------------------------------------------

species = siteDetections[5:96]

variableImportance <- data.frame()

# loop it
set.seed(123)
for (i in seq_along(species)) {
  
  model <- plsr(species[[i]] ~ ., data = spatVars, scale = TRUE, validation = "CV")
  
  optimalCompNo <- which.min(RMSEP(model, estimate = "adjCV")$val) -1
  
  variableImportance.tmp <- stack(VIP(model, opt.comp = optimalCompNo))
  colnames(variableImportance.tmp) <- c("variableImportance", "variable")
  variableImportance.tmp$species <- names(species[i])
  
  variableImportance <- rbind(variableImportance, variableImportance.tmp)
}

variableImportance <- merge(variableImportance, tidy, by.x= "variable", by.y= "VariableName")

# plot it

ggplot(variableImportance, aes(x=reorder(variable, variableImportance, FUN=mean), y=variableImportance, fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal()

ggplot(variableImportance, aes(x=reorder(variable, variableImportance), y=variableImportance, fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_col() + 
  coord_flip() +
  xlab("Habitat Variable") +
  ylab("Variable Importance") +
  theme_minimal()

# try NFPD correction ----------------------------------------------------------

set.seed(123)

model <- plsr(NFPD(siteDetections$Western.Kingbird) ~ ., data = spatVars, scale = TRUE, validation = "CV")

summary(model)
optimalCompNo <- which.min(RMSEP(model, estimate = "adjCV")$val[-1])
optimalCompNo

validationplot(model, val.type = "RMSEP")

cumVar <- R2(model, estimate = "train")$val[1, , ]
plot(seq_along(cumVar), cumVar, type = "b",
     xlab = "Number of Components",
     ylab = "Cumulative % Variance Explained in Y (training)")

cumVar <- R2(model)$val[1, , ]
plot(seq_along(cumVar), cumVar, type = "b",
     xlab = "Number of Components",
     ylab = "Cross-validated R-squared")

paste( "R2 =", round(R2(model)$val[1, , optimalCompNo+1], 2))

variableImportance.tmp <- stack(VIP(model, opt.comp = optimalCompNo))
colnames(variableImportance.tmp) <- c("variableImportance", "variable")
variableImportance.tmp <- merge(variableImportance.tmp, tidy, by.x= "variable", by.y= "VariableName")

ggplot(variableImportance.tmp, aes(x=reorder(variable, variableImportance), y=variableImportance, fill=Category)) +
  scale_fill_manual(values = c("#309898", "#8A2D3B", "#A0C878", "#27548A", "#F5C45E", "#644A07", "#8E7DBE")) +
  geom_col() + 
  coord_flip() +
  xlab("Habitat Variable") +
  ylab("Variable Importance") +
  theme_minimal()

pred <- predict(model, ncomp = optimalCompNo)

plot(NFPD(siteDetections$Western.Kingbird), pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "PLSR: Observed vs Predicted") +
  abline(0, 1, lty = 2) +
  legend("topright", legend = paste0("cor = ", round(cor(NFPD(siteDetections$Western.Kingbird), pred), 2)), bty = "n")


ggplot(siteDetections, aes(x=NFPD(Acorn.Woodpecker), y=gsp)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = lm)

# loop through using NFPD ------------------------------------------------------

variableImportance <- data.frame()
modelParameters <- data.frame()
set.seed(123)
for (i in seq_along(species)) {
  
  model <- plsr(NFPD(species[[i]]) ~ ., data = spatVars, scale = TRUE, validation = "CV")
  
  optimalCompNo <- which.min(RMSEP(model, estimate = "adjCV")$val[-1])
  
  varExp <- sum(explvar(model)[1:optimalCompNo])
  R2_train_opt <- R2(model, estimate = "train")$val[1, , ][optimalCompNo+1]
  R2_CV_opt <- R2(model)$val[1, , ][optimalCompNo+1]
  R2_train_max <- max(R2(model, estimate = "train")$val[1, , ])
  R2_CV_max <- max(R2(model)$val[1, , ])
  R2_CV_max_comps <- which.max(R2(model)$val[1, , ]) -1
  
  variableImportance.tmp <- stack(VIP(model, opt.comp = optimalCompNo))
  colnames(variableImportance.tmp) <- c("variableImportance", "variable")
  variableImportance.tmp$species <- names(species[i])
  variableImportance.tmp$varSet <- "spatVars"
  
  variableImportance.tmp <- variableImportance.tmp %>%
    mutate(variableRank = rank(-variableImportance, ties.method = "first"))
  
  variableImportance <- rbind(variableImportance, variableImportance.tmp)
  
  modelParameters.tmp <- data.frame(species = names(species[i]),
                                    varSet = "spatVars",
                                    zeros = sum(species[i] == 0),
                                    optimalCompNo = optimalCompNo,
                                    R2_train_opt = R2_train_opt,
                                    R2_train_max = R2_train_max,
                                    R2_CV_opt = R2_CV_opt,
                                    R2_CV_max = R2_CV_max,
                                    R2_CV_max_comps = R2_CV_max_comps)
  
  modelParameters <- rbind(modelParameters, modelParameters.tmp)
  
}

variableImportance <- merge(variableImportance, tidy, by.x= "variable", by.y= "VariableName")


ggplot(modelParameters, aes(y=R2_CV, x=zeros)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method=lm) +
  ylab("Cross-validated R-squared")

ggplot(variableImportance, aes(x=variableRank)) +
  geom_histogram() +
  facet_wrap(~variable)

# single var, single species ---------------------------------------------------

set.seed(123)

model <- plsr(NFPD(siteDetections$Acorn.Woodpecker) ~ ., data = siteDetections[97], scale = TRUE, validation = "CV")


# make a double loop -----------------------------------------------------------

# make the list of varSets
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

# make empty data frames
variableImportance <- data.frame()
modelParameters <- data.frame()

set.seed(123)

# loop through
for(vs in names(varSets)) {
  
  spatVars <- varSets[[vs]]
  
  for (i in seq_along(species)) {
    
    model <- plsr(NFPD(species[[i]]) ~ ., data = spatVars, scale = TRUE, validation = "CV")
    
    optimalCompNo <- which.min(RMSEP(model, estimate = "adjCV")$val[-1])
    
    varExp <- sum(explvar(model)[1:optimalCompNo])
    R2_train_opt <- R2(model, estimate = "train")$val[1, , ][optimalCompNo+1]
    R2_CV_opt <- R2(model)$val[1, , ][optimalCompNo+1]
    R2_train_max <- max(R2(model, estimate = "train")$val[1, , ])
    R2_CV_max <- max(R2(model)$val[1, , ])
    R2_CV_max_comps <- which.max(R2(model)$val[1, , ]) -1
    
    variableImportance.tmp <- stack(VIP(model, opt.comp = optimalCompNo))
    colnames(variableImportance.tmp) <- c("variableImportance", "variable")
    variableImportance.tmp$species <- names(species[i])
    variableImportance.tmp$varSet <- vs
    
    variableImportance.tmp <- variableImportance.tmp %>%
      mutate(variableRank = rank(-variableImportance, ties.method = "first"))
    
    variableImportance <- rbind(variableImportance, variableImportance.tmp)
    
    modelParameters.tmp <- data.frame(species = names(species[i]),
                                      varSet = vs,
                                      zeros = sum(species[i] == 0),
                                      optimalCompNo = optimalCompNo,
                                      R2_train_opt = R2_train_opt,
                                      R2_train_max = R2_train_max,
                                      R2_CV_opt = R2_CV_opt,
                                      R2_CV_max = R2_CV_max,
                                      R2_CV_max_comps = R2_CV_max_comps)
    
    modelParameters <- rbind(modelParameters, modelParameters.tmp)
    
  }
  
}

# add category names
variableImportance <- merge(variableImportance, tidy, by.x= "variable", by.y= "VariableName")

# save it
write.csv(variableImportance, "variableImportance_PLSR_NFPD_20251117.csv")
write.csv(modelParameters, "modelParameters_PLSR_NFPD_20251117.csv")

# plot it

ggplot(modelParameters, aes(y=species, x=R2_CV_max, color=varSet)) +
  geom_point() +
   theme_minimal()
  

