
# Trait Importance
# 16 July 2025
# getting started with the trait importance project.

library(randomForest)
library(ggplot2)
library(sf)
library(basemaps)
library(mapview)
library(patchwork)
library(tidyr)
library(dplyr)
library(mgcv)
library(VSURF)
library(caret)
library(doParallel)
library(ranger)
library(Boruta)
library(stringr)
library(forcats)
library(yacca)
library(eulerr)
library(purrr)



################################################################################
# IMPORTANT BEGINING THINGS. ALWAYS RUN THIS SECTION
################################################################################

# setwd ------------------------------------------------------------------------
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# load data --------------------------------------------------------------------
siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# name variable groupings ------------------------------------------------------
names(siteDetections)
species <- colnames(siteDetections)[5:96]
spatVars <- colnames(siteDetections)[97:182]
traitVars <- colnames(siteDetections)[97:108]
anthrVars <- colnames(siteDetections)[109:123]
climVars <- colnames(siteDetections)[124:151]
fxnVars <- colnames(siteDetections)[c(152:154, 179:182)]
phenoVars <- colnames(siteDetections)[155:164]
topoVars <- colnames(siteDetections)[165:178]
notTraits <- colnames(siteDetections)[109:182]

################################################################################
# END OF PRIMER
################################################################################

# run random forest, make some test plots --------------------------------------
species
selectedSpecies <- species[12]

predictors <- paste(spatVars, collapse = " + ")
rf_formula <- as.formula(paste(selectedSpecies, "~", predictors))
model <- randomForest(rf_formula, data = siteDetections, importance = TRUE)
model
percent_var_explained <- model$rsq[length(model$rsq)] * 100


# Create a dataframe for plotting
plot_df <- data.frame(
  Observed = siteDetections[[selectedSpecies]],
  Predicted = model$predicted
)

# Plot accuracy
p1 <- ggplot(plot_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Variance explained = ", round(percent_var_explained, 2), "%")) +
  labs(title = paste("Observed vs Predicted", selectedSpecies),
       x = "Observed Detection Rate",
       y = "Predicted Detection Rate") +
  theme_minimal()

# plot histogram
p2 <- ggplot(siteDetections, aes(x = .data[[selectedSpecies]])) +
  geom_histogram() +
  theme_minimal()

# plot map
p3 <- ggplot() +
  geom_sf(data = st_transform(st_as_sf(siteDetections), crs = 3857), aes(color = .data[[selectedSpecies]]), size = 2, alpha = 0.6) +
  scale_color_viridis_c(option = "plasma", name = paste("Detection:", selectedSpecies)) +
  labs(title = paste("Detection Map of", selectedSpecies)) +
  theme_minimal() +
  theme(legend.position = "none")

# combine panels
(p1 / p2 + plot_layout(heights = c(4, 1))) | p3 + plot_layout(widths = c(2, 1))

# Get the variable importance table
imp <- importance(model)

# Sort it by %IncMSE in descending order
imp[order(imp[, "%IncMSE"], decreasing = TRUE), ]

# loop it  ---------------------------------------------------------------------



# Initialize results dataframe with species column
PercentVarExplained <- data.frame(Species = species, stringsAsFactors = FALSE)

# Define function to calculate % variance explained
get_rf_rsq <- function(response, predictors) {
  rf_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  model <- randomForest(rf_formula, data = siteDetections, importance = FALSE)
  tail(model$rsq, 1) * 100
}

# First: all spatial variables
spat_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, spatVars))
PercentVarExplained$spatVars <- round(spat_rsqs, 2)

# Then: trait variables
trait_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, traitVars))
PercentVarExplained$traitVars <- round(trait_rsqs, 2)

# Anthr variables
anthr_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, anthrVars))
PercentVarExplained$anthrVars <- round(anthr_rsqs, 2)

# climVars
clim_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, climVars))
PercentVarExplained$climVars <- round(clim_rsqs, 2)

# fxnVars
fxn_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, fxnVars))
PercentVarExplained$fxnVars <- round(fxn_rsqs, 2)

# phenoVars
pheno_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, phenoVars))
PercentVarExplained$phenoVars <- round(pheno_rsqs, 2)

# topoVars
topo_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, topoVars))
PercentVarExplained$topoVars <- round(topo_rsqs, 2)

# notTraits
notTraits_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, notTraits))
PercentVarExplained$notTraits <- round(notTraits_rsqs, 2)


# Reshape the data to long format
long_df <- PercentVarExplained %>%
  pivot_longer(cols = c(spatVars, traitVars, anthrVars, climVars, fxnVars), names_to = "VariableGroup", values_to = "VarianceExplained")

# Plot with density and legend
ggplot(data = long_df, aes(x = VarianceExplained, fill = VariableGroup)) +
  geom_density(alpha = 0.4, color = NA) +
  geom_vline(data = long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE)),
             aes(xintercept = mean_val, color = VariableGroup), linetype = "dashed", size = 1) +
  scale_fill_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", anthrVars = "#901E3E", climVars = "#8DBCC7"), name = "Variable Group") +
  scale_color_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", anthrVars = "#901E3E", climVars = "#8DBCC7"), guide = "none") +
  labs(x = "% Variance Explained", y = "Density") +
  theme_minimal()

# means
long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE))

# how much do traits improve the accuracy? -------------------------------------

# Reshape the data to long format
long_df <- PercentVarExplained %>%
  pivot_longer(cols = c(spatVars, traitVars, notTraits), names_to = "VariableGroup", values_to = "VarianceExplained")

# Plot with density and legend
ggplot(data = long_df, aes(x = VarianceExplained, fill = VariableGroup)) +
  geom_density(alpha = 0.4, color = NA) +
  geom_vline(data = long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE)),
             aes(xintercept = mean_val, color = VariableGroup), linetype = "dashed", size = 1) +
  scale_fill_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", notTraits = "#901E3E"), name = "Variable Group") +
  scale_color_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", notTraits = "#901E3E"), guide = "none") +
  labs(x = "% Variance Explained", y = "Density") +
  theme_minimal()

# means
long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE))


# Presence / Absence models ----------------------------------------------------
# accuracy seems low overall
# going to try a binary presence/absence version

# Convert species columns to binary presence/absence
binaryDetections <- siteDetections
binaryDetections[species] <- lapply(binaryDetections[species], function(x) {
  factor(ifelse(x > 0, "present", "absent"), levels = c("absent", "present"))
})

# try one species

selectedSpecies <- species[8]

predictors <- paste(spatVars, collapse = " + ")
rf_formula <- as.formula(paste(selectedSpecies, "~", predictors))
model <- randomForest(rf_formula, data = binaryDetections, importance = TRUE)
model
accuracy <- 1 - model$err.rate[nrow(model$err.rate), "OOB"]
accuracy

# do the same chunk as above but with binary data

# Initialize results dataframe with species column
PercentVarExplained_binary <- data.frame(Species = species, stringsAsFactors = FALSE)

# Define function to calculate % variance explained
get_rf_acc <- function(response, predictors) {
  response_vals <- binaryDetections[[response]]
  if (length(unique(response_vals)) < 2) return(NA) # escape for species always present or absent
  rf_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  model <- randomForest(rf_formula, data = binaryDetections, importance = FALSE)
  1 - model$err.rate[nrow(model$err.rate), "OOB"]
}


# First: all spatial variables
spat_rsqs <- sapply(species, function(sp) get_rf_acc(sp, spatVars))
PercentVarExplained_binary$spatVars <- spat_rsqs

# Then: trait variables
trait_rsqs <- sapply(species, function(sp) get_rf_acc(sp, traitVars))
PercentVarExplained_binary$traitVars <- round(trait_rsqs, 2)

# Anthr variables
anthr_rsqs <- sapply(species, function(sp) get_rf_acc(sp, anthrVars))
PercentVarExplained_binary$anthrVars <- round(anthr_rsqs, 2)

# climVars
clim_rsqs <- sapply(species, function(sp) get_rf_acc(sp, climVars))
PercentVarExplained_binary$climVars <- round(clim_rsqs, 2)

# fxnVars
fxn_rsqs <- sapply(species, function(sp) get_rf_acc(sp, fxnVars))
PercentVarExplained_binary$fxnVars <- round(fxn_rsqs, 2)

# Reshape the data to long format
long_df <- PercentVarExplained_binary %>%
  pivot_longer(cols = c(spatVars, traitVars, anthrVars, climVars, fxnVars), names_to = "VariableGroup", values_to = "VarianceExplained")

# Plot with density and legend
ggplot(data = long_df, aes(x = VarianceExplained, fill = VariableGroup)) +
  geom_density(alpha = 0.4, color = NA) +
  geom_vline(data = long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE)),
             aes(xintercept = mean_val, color = VariableGroup), linetype = "dashed", size = 1) +
  scale_fill_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", anthrVars = "#901E3E", climVars = "#8DBCC7"), name = "Variable Group") +
  scale_color_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", anthrVars = "#901E3E", climVars = "#8DBCC7"), guide = "none") +
  labs(x = "% Variance Explained", y = "Density") +
  theme_minimal()

# means
long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE))




# find out which is the best predictor for each species ------------------------

PercentVarExplained$bestPredictor <- apply(PercentVarExplained[, c("traitVars", "climVars", "anthrVars", "fxnVars", "phenoVars", "topoVars")], 1, function(x) {
  names(x)[which.max(x)]
})

# Count the frequency of each best predictor type
pie_data <- PercentVarExplained %>%
  count(bestPredictor) %>%
  mutate(perc = n / sum(n) * 100,
         label = paste0(bestPredictor, "\n", round(perc, 1), "%"),
         ymax = cumsum(perc),
         ymin = lag(ymax, default = 0),
         label_pos = (ymin + ymax) / 2)

# Create the pie chart with internal labels
ggplot(pie_data, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = bestPredictor)) +
  geom_rect() +
  geom_text(aes(x = 0.5, y = label_pos, label = label), size = 6, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Best Predictor Type per Species") +
  scale_fill_brewer(palette = "Set2", guide = "none")


# among species best predicted by traits, how much do traits improve the prediction? ---------------------------------

# Reshape the data to long format
long_df <- PercentVarExplained %>%
  pivot_longer(cols = c(spatVars, traitVars, notTraits), names_to = "VariableGroup", values_to = "VarianceExplained")

# only traits favorites
long_df <- subset(long_df, bestPredictor == "traitVars")

# Plot with density and legend
ggplot(data = long_df, aes(x = VarianceExplained, fill = VariableGroup)) +
  geom_density(alpha = 0.4, color = NA) +
  geom_vline(data = long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE)),
             aes(xintercept = mean_val, color = VariableGroup), linetype = "dashed", size = 1) +
  scale_fill_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", notTraits = "#901E3E"), name = "Variable Group") +
  scale_color_manual(values = c(spatVars = "#EFE4D2", traitVars = "#347433", notTraits = "#901E3E"), guide = "none") +
  labs(x = "% Variance Explained", y = "Density") +
  theme_minimal()

# means
long_df %>% group_by(VariableGroup) %>% summarise(mean_val = mean(VarianceExplained, na.rm = TRUE))




################################################################################
# Updating methods. Variable selection
################################################################################

# run randomForest 10x
# ID which vars consistently make the model worse
# remove those
# re-run random forest

# Choose species outside the loop
species
selectedSpecies <- species[12]

# Number of repetitions
n_reps <- 10

# choose vars
names(PercentVarExplained)
selectedVars <- traitVars

# Initialize matrix to store variable importances
importance_mat <- matrix(NA, nrow = length(selectedVars), ncol = n_reps,
                         dimnames = list(selectedVars, paste0("Run_", 1:n_reps)))

# Run multiple forests and collect %IncMSE
for (i in 1:n_reps) {
  rf_formula <- as.formula(paste(selectedSpecies, "~", paste(selectedVars, collapse = " + ")))
  model_i <- randomForest(rf_formula, data = siteDetections, importance = TRUE)
  imp <- importance(model_i, type = 1)  # %IncMSE
  importance_mat[rownames(imp), i] <- imp[, 1]
}

# Calculate mean and consistency of negative importances
mean_importance <- rowMeans(importance_mat, na.rm = TRUE)
neg_freq <- rowSums(importance_mat < 0, na.rm = TRUE)

# Flag variables with negative importance in most runs
# You can change the threshold here (e.g., > 7 of 10 times)
remove_vars <- names(neg_freq[neg_freq >= 0.7 * n_reps])

# Print summary
cat("Variables consistently unimportant (negative %IncMSE):\n")
print(remove_vars)

# Rerun random forest without those variables
final_vars <- setdiff(selectedVars, remove_vars)
rf_formula_final <- as.formula(paste(selectedSpecies, "~", paste(final_vars, collapse = " + ")))
final_model <- randomForest(rf_formula_final, data = siteDetections, importance = TRUE)
percent_var_explained <- final_model$rsq[length(final_model$rsq)] * 100

# Output
list(
  removed_variables = remove_vars,
  final_model = final_model,
  importance_across_runs = importance_mat
)

# Create a dataframe for plotting
plot_df <- data.frame(
  Observed = siteDetections[[selectedSpecies]],
  Predicted = final_model$predicted
)

# Plot accuracy
ggplot(plot_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Variance explained = ", round(percent_var_explained, 2), "%")) +
  labs(title = paste("Observed vs Predicted", selectedSpecies),
       x = "Observed Detection Rate",
       y = "Predicted Detection Rate") +
  theme_minimal()

# Get the variable importance table
imp <- importance(final_model)

# Sort it by %IncMSE in descending order
imp[order(imp[, "%IncMSE"], decreasing = TRUE), ]


################################################################################
# Single traits, single species
################################################################################

# which individual variables are best at predicting the detection rates of individual species?

# Fractional detection rates
FPD <- function(species_col) {species_col / max(species_col)}
meanFPD <- function(species_col) {mean(FPD(species_col)[FPD(species_col)>0])}
NFPD <- function(species_col) {(FPD(species_col))^(log(0.5)/log(meanFPD(species_col)))}

# exploration

hist(FPD(siteDetections$American.Robin))
meanFPD(siteDetections$American.Robin)
hist(NFPD(siteDetections$American.Robin))

predictor <- spatVars[80]
species
selectedSpecies <- species[12]


ggplot(siteDetections, aes(x=get(selectedSpecies), y=get(predictor))) +
  geom_point() +
  stat_smooth(method=lm)

ggplot(siteDetections, aes(x=NFPD(get(selectedSpecies)), y=get(predictor))) +
  geom_point() +
  stat_smooth(method=lm)

ggplot(siteDetections, aes(x=NFPD(get(selectedSpecies)), y=get(predictor))) +
  geom_point() +
  stat_smooth()

ggplot(siteDetections, aes(x=NFPD(Mountain.Chickadee), y=DTM)) +
  geom_point() +
  stat_smooth() +
  theme_minimal()

ggplot(siteDetections, aes(x=Mountain.Chickadee, y=DTM)) +
  geom_point() +
  stat_smooth() +
  theme_minimal()

ggplot(siteDetections, aes(x=NFPD(Yellow.rumped.Warbler), y=Sulfur)) +
  geom_point() +
  stat_smooth() +
  theme_minimal()

ggplot(siteDetections, aes(x=Yellow.rumped.Warbler, y=Sulfur)) +
  geom_point() +
  stat_smooth() +
  theme_minimal()

rf_formula <- as.formula(paste("NFPD(", selectedSpecies, ")~", predictor))
model <- randomForest(rf_formula, data = siteDetections, importance = TRUE)
model

model <- gam(rf_formula, data = siteDetections)
summary(model)

model <- glm(rf_formula, data = siteDetections)
summary(model)

model <- lm(rf_formula, data = siteDetections)
summary(model)

################################################################################
################################################################################
################################################################################

# Create all species-variable pairs
pairs <- expand.grid(species = species, var = spatVars, stringsAsFactors = FALSE)

# Function to fit GAM and extract R2
get_r2 <- function(sp, v) {
  f <- as.formula(paste("NFPD(", sp, ")~", v))
  mod <- gam(f, data = siteDetections)
  summary(mod)$r.sq
}

# Apply function to each pair
pairs$r2 <- mapply(get_r2, pairs$species, pairs$var)

# Define category membership
varCats <- list(
  trait = colnames(siteDetections)[97:108],
  anthr = colnames(siteDetections)[109:123],
  clim  = colnames(siteDetections)[124:151],
  fxn   = colnames(siteDetections)[c(152:154, 179:182)],
  pheno = colnames(siteDetections)[155:164],
  topo  = colnames(siteDetections)[165:178]
)

# Convert to a named vector lookup
cat_lookup <- stack(varCats) |> 
  setNames(c("var", "category"))

# Join category into pairs
pairs <- pairs %>% left_join(cat_lookup, by = "var")


ggplot(pairs, aes(x = r2, y = category)) +
  geom_violin() +
  theme_minimal()

# compute means per category
pairs %>% 
  group_by(category) %>% 
  summarise(mean_r2 = mean(r2, na.rm = TRUE))


# add niche areas to dataframe
nicheAreas <- readRDS("data/speciesNicheAreas_20250624.rds")
pairs <- merge(pairs, nicheAreas)

# compare niche size to r-squared
ggplot(pairs, aes(x=area, y=r2, color=category)) +
  geom_point(aes(alpha=0.5)) +
  theme_minimal()

ggplot(pairs, aes(x=area, y=r2, color=category)) +
  stat_smooth(method=lm) +
  theme_minimal()

################################################################################
################################################################################
################################################################################

#---------------------------
# Step 0: Define variables
#---------------------------
response <- siteDetections$Woodhouse.s.Scrub.Jay
predictors <- siteDetections[ , spatVars]

#---------------------------
# Step 1: Full Random Forest
#---------------------------
set.seed(124)
rf.full <- randomForest(x = predictors, y = response, ntree = 1000, importance = TRUE)
print(rf.full)

# Out-of-bag R² (last value is the final)
full_r2 <- rf.full$rsq[length(rf.full$rsq)]
cat("Full RF OOB R²:", full_r2, "\n")

#---------------------------
# Step 2: Variable Selection with VSURF
#---------------------------
# Detect cores and register
set.seed(124)
vsurf.res <- VSURF(x = predictors, y = response) 

# Best variables for prediction
sel.vars <- vsurf.res$varselect.pred
cat("Variables kept:", length(sel.vars), "\n")

#---------------------------
# Step 3: Fit reduced RF
#---------------------------
predictors.sel <- predictors[ , sel.vars, drop=FALSE]
names(predictors.sel)

set.seed(124)
rf.sel <- randomForest(x = predictors.sel, y = response, ntree = 1000, importance = TRUE)
print(rf.sel)

sel_r2 <- rf.sel$rsq[length(rf.sel$rsq)]
cat("Reduced RF OOB R²:", sel_r2, "\n")

#---------------------------
# Step 4: Cross-validation with caret
#---------------------------
set.seed(124)
ctrl <- trainControl(method = "cv", number = 10)

cv.rf <- train(response ~ ., data = data.frame(response=response, predictors.sel),
               method = "rf",
               trControl = ctrl,
               tuneGrid = data.frame(mtry = floor(sqrt(ncol(predictors.sel)))),
               ntree = 1000)

cat("Cross-validated R²:", max(cv.rf$results$Rsquared), "\n")


################################################################################


#---------------------------
# Predictions from RF models
#---------------------------
pred.full <- predict(rf.full, predictors)         # full model predictions
pred.sel  <- predict(rf.sel, predictors.sel)      # reduced model predictions

#---------------------------
# Combine into one dataframe
#---------------------------
plot.df <- data.frame(
  Observed = response,
  PredFull = pred.full,
  PredSel  = pred.sel
)

#---------------------------
# Plot observed vs predicted (full model)
#---------------------------
ggplot(plot.df, aes(x = Observed, y = PredFull)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Observed vs Predicted (Full RF)",
       x = "Observed Detection Rate",
       y = "Predicted Detection Rate") +
  theme_minimal()

#---------------------------
# Plot observed vs predicted (reduced model)
#---------------------------
ggplot(plot.df, aes(x = Observed, y = PredSel)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Observed vs Predicted (Reduced RF)",
       x = "Observed Detection Rate",
       y = "Predicted Detection Rate") +
  theme_minimal()

################################################################################
################################################################################
################################################################################

# Parallel setup
n.cores <- parallel::detectCores() - 1
cl <- makeCluster(n.cores)
registerDoParallel(cl)

# Function to get pruned RF CV R² for a species + predictor set
get_cvR2 <- function(resp, predictors) {
  y <- siteDetections[[resp]]
  X <- siteDetections[, predictors, drop = FALSE]
  
  # VSURF for variable selection (pruned set)
  set.seed(123)
  vsurf.res <- VSURF(x = X, y = y)
  sel.vars <- vsurf.res$varselect.pred
  if (length(sel.vars) == 0) return(NA)   # no selected variables
  
  X.sel <- X[, sel.vars, drop = FALSE]
  
  # caret CV with randomForest
  set.seed(123)
  ctrl <- trainControl(method = "cv", number = 10)
  cv.rf <- train(y ~ ., data = data.frame(y=y, X.sel),
                 method = "rf", trControl = ctrl, ntree = 1000)
  
  max(cv.rf$results$Rsquared) # predictive ceiling (CV R²)
}

# Run in parallel over species
results <- foreach(sp = species, .combine = rbind, .packages = c("VSURF","caret","randomForest")) %dopar% {
  r2_spat <- get_cvR2(sp, spatVars)
  r2_notT <- get_cvR2(sp, notTraits)
  data.frame(species = sp,
             prunedRF_spatVars_R2 = r2_spat,
             prunedRF_notTraits_R2 = r2_notT)
}

stopCluster(cl)

# Final results dataframe
results <- as.data.frame(results)

# started 11:12 am
# still running 02:09
# finished 02:11

# wow that took a while. save it now D:

write.csv(results, "randomForest_VSURFpruned_traitsVnoTraits_R2_allSpecies_20250826.csv")
randomForest_VSURFpruned_traitsVnoTraits_R2_allSpecies_20250826 <- results

# plot R2 with vs R2 without traits

ggplot(results, aes(x=prunedRF_notTraits_R2, y=prunedRF_spatVars_R2)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope=1, intercept = 0) +
  xlab("Best possible R2 without traits") +
  ylab("Best possible R2 with traits") +
  theme_minimal()

################################################################################
################################################################################
################################################################################

# ok that was super slow and didn't give the results I was looking for
# well. it ran successfully. but I'm still getting lower R2 for models with more vars. so... ???


# Parallel setup
n.cores <- parallel::detectCores() - 1
cl <- makeCluster(n.cores)
registerDoParallel(cl)

# Function to get Boruta-pruned RF CV R²
get_cvR2_boruta <- function(resp, predictors) {
  y <- siteDetections[[resp]]
  X <- siteDetections[, predictors, drop = FALSE]
  df <- data.frame(y = y, X)
  
  # Boruta variable selection
  set.seed(123)
  bor <- Boruta(y ~ ., data = df, doTrace = 0)
  sel.vars <- getSelectedAttributes(bor, withTentative = FALSE)
  if (length(sel.vars) == 0) return(NA)  # no variables kept
  
  df.sel <- df[, c("y", sel.vars), drop = FALSE]
  
  # caret CV with ranger
  set.seed(123)
  ctrl <- trainControl(method = "cv", number = 10)
  cv.rf <- train(y ~ ., data = df.sel,
                 method = "ranger",
                 trControl = ctrl,
                 importance = "permutation",
                 num.trees = 1000)
  
  max(cv.rf$results$Rsquared)
}

# Run in parallel over species
results <- foreach(sp = species, .combine = rbind,
                   .packages = c("Boruta","caret","ranger")) %dopar% {
                     r2_spat <- get_cvR2_boruta(sp, spatVars)
                     r2_notT <- get_cvR2_boruta(sp, notTraits)
                     data.frame(species = sp,
                                prunedRF_spatVars_R2 = r2_spat,
                                prunedRF_notTraits_R2 = r2_notT)
                   }

stopCluster(cl)

# Final results
results <- as.data.frame(results)

# started 2:54 pm
# still running 4:03
# failed 4:07

################################################################################
################################################################################
################################################################################

# ok, so which species are worse with traits, and is it real or accidental/random?

results$difference <- results$prunedRF_spatVars_R2 - results$prunedRF_notTraits_R2

# Woodhouse.s.Scrub.Jay had some bad results
# lets check and see if it is always that way

# Species of interest
resp <- "Woodhouse.s.Scrub.Jay"
y <- siteDetections[[resp]]
X_spat <- siteDetections[, spatVars, drop = FALSE]
X_notT <- siteDetections[, notTraits, drop = FALSE]

# Function to run VSURF + RF + caret CV once
run_once <- function(X, y) {
  # VSURF for variable selection
  vsurf.res <- VSURF(x = X, y = y)
  sel.vars <- vsurf.res$varselect.pred
  if (length(sel.vars) == 0) return(NA)  # no variables kept
  
  X.sel <- X[, sel.vars, drop = FALSE]
  
  # caret CV Random Forest
  ctrl <- trainControl(method = "cv", number = 10)
  cv.rf <- train(y ~ ., data = data.frame(y=y, X.sel),
                 method = "rf", trControl = ctrl, ntree = 1000)
  
  max(cv.rf$results$Rsquared)
}

# Run 10 iterations (no seed → random variation preserved)
WoodhouseIterations <- data.frame(iteration = 1:10,
                      prunedRF_spatVars_R2 = NA,
                      prunedRF_notTraits_R2 = NA)

for (i in 1:10) {
  r2_spat <- run_once(X_spat, y)
  r2_notT <- run_once(X_notT, y)
  
  WoodhouseIterations$prunedRF_spatVars_R2[i] <- r2_spat
  WoodhouseIterations$prunedRF_notTraits_R2[i] <- r2_notT
}

WoodhouseIterations

# started 4:44
# iteration 1: 4:50 (6 minutes)
# leaving the office. 6*10=60 minutes predicted run time.

WoodhouseIterations_long <- WoodhouseIterations %>% pivot_longer(
  cols = c(prunedRF_spatVars_R2, prunedRF_notTraits_R2),
  names_to = "modelType",
  values_to = "R2")

  
ggplot(WoodhouseIterations_long, aes(x=R2, y=modelType)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Woodhouse's Scrub Jay")


################################################################################
################################################################################
################################################################################

# 27 Aug 2025
# something was cursed about my loop code
# I'm going to remove the carat cross validation step
# re-creating the loop code, now more carefully

# variables --------------------------------------------------------------------
species <- colnames(siteDetections)[5:96]
spatVars <- colnames(siteDetections)[97:182]
notTraits <- colnames(siteDetections)[109:182]
variableSets <- list(spatVars = spatVars, notTraits = notTraits)
variableSets <- c("spatVars", "notTraits")


# try it once ------------------------------------------------------------------
variableSet <- variableSets[1]
VS <- siteDetections[c(get(variableSet))]
selectedSpecies <- siteDetections[c("Acorn.Woodpecker")]

fullRFmodel <- randomForest(x = VS, y = selectedSpecies, ntree = 1000, importance = TRUE)
fullRF_R2 <- fullRFmodel$rsq[length(fullRFmodel$rsq)]
fullRF_mse <- fullRFmodel$mse[length(fullRFmodel$mse)]

VSURFvarSelection <- VSURF(x = VS, y = selectedSpecies) 

keptVariables <- VS[ , VSURFvarSelection$varselect.pred, drop=FALSE]
keptVariableNames <- names(keptVariables)

prunedRF <- randomForest(x = keptVariables, y = selectedSpecies, ntree = 1000, importance = TRUE)
prunedRF_R2 <- prunedRF$rsq[length(prunedRF$rsq)]
prunedRF_mse <- prunedRF$mse[length(prunedRF$mse)]

test <- data.frame(species = names(selectedSpecies),
           variableSet = variableSet,
           keptVariables = paste(keptVariableNames, collapse = ", "),
           prunedRF_R2 = prunedRF_R2,
           prunedRF_mse = prunedRF_mse)


# MAKE A FUNCTION --------------------------------------------------------------

get_R2 <- function(sp, vars, varset_name) {
  # Choose variables
  variableSet <- siteDetections[, vars] # choose either spatVars or notTraits
  selectedSpecies <- siteDetections[[sp]] # choose one species
  # Run VSURF
  VSURFvarSelection <- VSURF(x = variableSet, y = selectedSpecies) 
  # find best variables
  keptVariables <- variableSet[ , VSURFvarSelection$varselect.pred, drop=FALSE]
  keptVariableNames <- names(keptVariables)
  # Run RandomForest
  prunedRF <- randomForest(x = keptVariables, y = selectedSpecies, ntree = 1000, importance = TRUE)
  # Save R2 and mse
  prunedRF_R2 <- prunedRF$rsq[length(prunedRF$rsq)]
  prunedRF_mse <- prunedRF$mse[length(prunedRF$mse)]
  data.frame(species = sp,
             variableSet = varset_name,
             keptVariables = paste(keptVariableNames, collapse = ", "),
             prunedRF_R2 = prunedRF_R2,
             prunedRF_mse = prunedRF_mse)
}

# RUN IT IN PARALLEL -----------------------------------------------------------

# smol test sets
species <- colnames(siteDetections)[5:6]
species <- colnames(siteDetections)[5:96]

combos <- expand.grid(sp = species, varset_name = names(variableSets), stringsAsFactors = FALSE)
combos <- expand.grid(sp = species, varset_name = names(variableSets), iteration = 1:10, stringsAsFactors = FALSE)

# start parallelization
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)


allsp_pruned_R2 <- foreach(i = 1:nrow(combos), .combine = rbind, .packages = c("VSURF","randomForest")) %dopar% {
  sp <- combos$sp[i]
  varset_name <- combos$varset_name[i]
  vars <- variableSets[[varset_name]]
  get_R2(sp, vars, varset_name)
}

# end parallelization
stopCluster(cl)

# long loop code running 11:16 am 29 Aug 2025
# 2 Sep 2025 this loop code worked
# only issue, I forgot to add a column for "iterations". no big deal.
saveRDS(allsp_pruned_R2, "allsp_pruned_R2_20250902.rds")


# lets see what it looks like
x <- readRDS("allsp_pruned_R2_20250902.rds")

x <- x %>% mutate(traitVariables = str_count(keptVariables, paste0("\\b(", paste(traitVars, collapse = "|"), ")\\b")))
x$Traits <- ifelse(x$traitVariables>0, "TRUE", "FALSE")

AW <- subset(x, x$species=="Acorn.Woodpecker")
AW$Traits <- ifelse(AW$traitVariables>0, "TRUE", "FALSE")
ggplot(AW, aes(x=prunedRF_R2, fill = Traits)) +
  geom_histogram() +
  theme_minimal()

OCW <- subset(x, x$species=="Orange.crowned.Warbler")
ggplot(OCW, aes(x=prunedRF_mse, fill = Traits)) +
  geom_histogram() +
  theme_minimal()

test <- x %>%
  group_by(species, Traits) %>%
  summarise(prunedRF_R2_mean = mean(prunedRF_R2),
            prunedRF_mse_mean = mean(prunedRF_mse),
            prunedRF_R2_sd = sd(prunedRF_R2)) %>%
  ungroup()

test <- test %>%
  mutate(species = fct_reorder(species, prunedRF_R2_mean, .fun = min))

ggplot(test, aes(x=prunedRF_R2_mean, y=species, colour = Traits)) +
  geom_point() +
  theme_minimal() +
  geom_errorbar(aes(xmin=(prunedRF_R2_mean - prunedRF_R2_sd), xmax=(prunedRF_R2_mean + prunedRF_R2_sd)))


test <- test %>%
  select(-prunedRF_mse_mean, -prunedRF_R2_sd) %>%
  pivot_wider(
    names_from = Traits,
    values_from = c(prunedRF_R2_mean),
    names_sep = "_"
  )

test$difference <- test$`TRUE` - test$`FALSE`
mean(test$difference, na.rm = TRUE)

# how often does each variable appear?

counts <- x %>%
  subset(variableSet=="spatVars") %>%
  # split comma-separated values into long format
  separate_rows(keptVariables, sep = ",\\s*") %>%
  # count occurrences
  count(keptVariables, sort = TRUE)

counts
sum(counts$n)


counts <- x %>%
  # split comma-separated values into long format
  separate_rows(keptVariables, sep = ",\\s*") %>%
  # count occurrences
  count(keptVariables, sort = TRUE)

counts
sum(counts$n)

nrow(subset(x, x$variableSet=="spatVars"))

length(spatVars)

# I have 920 models that could include traits
# collectively those models have 4953 variables
# on average, that is ~5 variables per model
# I have 86 spatial variables
# 4953/86 = 57.59
# If variables were selected randomly, each var should appear ~58 times


NicheAreas <- readRDS("data/speciesNicheAreas_20250624.rds")

# any relationship between niche area and trait usefullness?

speciesImprovement <- merge(test, NicheAreas)
speciesImprovement$difference[is.na(speciesImprovement$difference)] <- 0

ggplot(speciesImprovement, aes(x=difference, y=area)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

ggplot(speciesImprovement, aes(y=difference, x=area)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

# what about habitat preference?

speciesPolygons <-  readRDS("data/speciesPolygons_jittered_detPerDay_20250624.rds")
speciesCentroids <- st_centroid(speciesPolygons)
spPC <- as.data.frame(st_coordinates(speciesCentroids))
speciesCentroids <- cbind(speciesCentroids, spPC)

speciesImprovement <- merge(speciesImprovement, speciesCentroids)

ggplot(speciesImprovement, aes(y=difference, x=X)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

ggplot(speciesImprovement, aes(x=difference, y=`FALSE`)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()


################################################################################
################################################################################
################################################################################

# 5 Sep 2025
# ok now I want to compare models made with only traits to models made with only climate
# I'll reuse my previous successful code, but this time use different variable sets.


# variables --------------------------------------------------------------------
species <- colnames(siteDetections)[5:96]
spatVars <- colnames(siteDetections)[97:182]
notTraits <- colnames(siteDetections)[109:182]
traitVars <- colnames(siteDetections)[97:108]
anthrVars <- colnames(siteDetections)[109:123]
climVars <- colnames(siteDetections)[124:151]
fxnVars <- colnames(siteDetections)[c(152:154, 179:182)]
phenoVars <- colnames(siteDetections)[155:164]
topoVars <- colnames(siteDetections)[165:178]
variableSets <- list(traitVars = traitVars, 
                     climVars = climVars,
                     fxnVars = fxnVars,
                     phenoVars = phenoVars,
                     topoVars = topoVars,
                     anthrVars = anthrVars)


# MAKE A FUNCTION --------------------------------------------------------------

get_R2 <- function(sp, vars, varset_name, iteration) {
  # Choose variables
  variableSet <- siteDetections[, vars] # choose either spatVars or notTraits
  selectedSpecies <- siteDetections[[sp]] # choose one species
  # Run VSURF
  VSURFvarSelection <- VSURF(x = variableSet, y = selectedSpecies) 
  kept_idx <- VSURFvarSelection$varselect.pred
  
  # If no variables were selected, return NA row
  if (length(kept_idx) == 0) {
    return(data.frame(species = sp,
                      variableSet = varset_name,
                      iteration = iteration,
                      keptVariables = NA,
                      prunedRF_R2 = NA,
                      prunedRF_mse = NA))
  }
  
  # Otherwise proceed with RF
  keptVariables <- variableSet[, kept_idx, drop = FALSE]
  keptVariableNames <- names(keptVariables)
  # Run RandomForest
  prunedRF <- randomForest(x = keptVariables, y = selectedSpecies, ntree = 1000, importance = TRUE)
  # Save R2 and mse
  prunedRF_R2 <- prunedRF$rsq[length(prunedRF$rsq)]
  prunedRF_mse <- prunedRF$mse[length(prunedRF$mse)]
  data.frame(species = sp,
             variableSet = varset_name,
             iteration = iteration,
             keptVariables = paste(keptVariableNames, collapse = ", "),
             prunedRF_R2 = prunedRF_R2,
             prunedRF_mse = prunedRF_mse)
}

# RUN IT IN PARALLEL -----------------------------------------------------------

# smol test sets
#species <- colnames(siteDetections)[5:6]
species <- colnames(siteDetections)[5:96]

combos <- expand.grid(sp = species, varset_name = names(variableSets), iteration = 1:10, stringsAsFactors = FALSE)

# start parallelization
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)


allsp_pruned_R2_sep5 <- foreach(i = 1:nrow(combos), .combine = rbind, .packages = c("VSURF","randomForest")) %dopar% {
  sp <- combos$sp[i]
  varset_name <- combos$varset_name[i]
  vars <- variableSets[[varset_name]]
  iteration <- combos$iteration[i]
  get_R2(sp, vars, varset_name, iteration)
}

# end parallelization
stopCluster(cl)

################################################################################
# Mon 8 sep 2025
# looks like my loop worked
# lets take a look


allsp_pruned_R2_sep5_summarized <- allsp_pruned_R2_sep5 %>%
  group_by(species, variableSet) %>%
  summarise(prunedRF_R2_mean = mean(prunedRF_R2, na.rm = TRUE),
            prunedRF_R2_sd = sd(prunedRF_R2, na.rm = TRUE)) %>%
  ungroup()

allsp_pruned_R2_sep5_summarized <- allsp_pruned_R2_sep5_summarized %>%
  mutate(species = fct_reorder(species, prunedRF_R2_mean, .fun = max, .na_rm = TRUE))

ggplot(allsp_pruned_R2_sep5_summarized, aes(x=prunedRF_R2_mean, y=species, colour = variableSet)) +
  geom_point() +
  theme_minimal() +
  geom_errorbar(aes(xmin=(prunedRF_R2_mean - prunedRF_R2_sd), xmax=(prunedRF_R2_mean + prunedRF_R2_sd)))


allsp_pruned_R2_sep5_summarized <- merge(allsp_pruned_R2_sep5_summarized, NicheAreas)


ggplot(allsp_pruned_R2_sep5_summarized, aes(x=prunedRF_R2_mean, y=area, colour = variableSet)) +
  geom_point() +
  theme_minimal() +
  geom_errorbar(aes(xmin=(prunedRF_R2_mean - prunedRF_R2_sd), xmax=(prunedRF_R2_mean + prunedRF_R2_sd))) +
  ylab("Niche Area")


################################################################################

# all species, single variable random forest

species <- colnames(siteDetections)[5:96]

combos <- expand.grid(sp = species, variable = spatVars, iteration = 1:10, stringsAsFactors = FALSE)

get_single_R2 <- function(sp, variable, iteration) {
  # Choose variables
  selectedVariable <- siteDetections[variable] # choose one variable
  selectedSpecies <- siteDetections[[sp]] # choose one species
  # Run RandomForest
  RF <- randomForest(x = selectedVariable, y = selectedSpecies, ntree = 1000, importance = TRUE)
  # Save R2 and mse
  RF_R2 <- RF$rsq[length(RF$rsq)]
  RF_mse <- RF$mse[length(RF$mse)]
  data.frame(species = sp,
             variableSet = variable,
             iteration = iteration,
             RF_R2 = RF_R2,
             RF_mse = RF_mse)
}


# start parallelization
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)


allsp_singleVar_R2_sep10 <- foreach(i = 1:nrow(combos), .combine = rbind, .packages = c("randomForest")) %dopar% {
  sp <- combos$sp[i]
  variable <- combos$variable[i]
  iteration <- combos$iteration[i]
  get_single_R2(sp, variable, iteration)
}

# end parallelization
stopCluster(cl)

################################################################################

# 15 Sep 2025
# Monday

# trait vars vs clim vars
TraitVClim <- allsp_pruned_R2_sep5_summarized
TraitVClim <- subset(TraitVClim, TraitVClim$variableSet=="climVars" | TraitVClim$variableSet=="traitVars")


TraitVClim <- TraitVClim %>%
  pivot_wider(
    names_from = variableSet,
    values_from = c(prunedRF_R2_mean, prunedRF_R2_sd),
    names_sep = "_"
  )

rangeShift <- read.csv("data/rangeExpansion.csv")

TraitVClim <- merge(TraitVClim, rangeShift)

ggplot(TraitVClim, aes(x=prunedRF_R2_mean_climVars, y=prunedRF_R2_mean_traitVars, color=invasive)) +
  geom_point() +
  geom_errorbar(aes(xmin=(prunedRF_R2_mean_climVars - prunedRF_R2_sd_climVars), xmax=(prunedRF_R2_mean_climVars + prunedRF_R2_sd_climVars))) +
  geom_errorbar(aes(ymin=(prunedRF_R2_mean_traitVars - prunedRF_R2_sd_traitVars), ymax=(prunedRF_R2_mean_traitVars + prunedRF_R2_sd_traitVars))) +
  geom_abline(slope=1) +
  ylab("R2 of best trait model") +
  xlab("R2 of best climate model") +
  scale_color_manual(values = c("black", "#4FB7B3")) +
  geom_text(x=0.5, y=0, label="Climate Limited", color="black") +
  geom_text(x=0, y=0.4, label="Habitat Limited", color="black") +
  theme_minimal() 



################################################################################

# Lets try some cannonical correlation analysis


traits <- siteDetections[, traitVars]
climate  <- siteDetections[, climVars]

cc <- cancor(traits, climate)
cc$cor        # canonical correlations
cc$cor^2      # shared variance per canonical dimension

cc <- cca(traits, climate)
summary(cc)

qr(traits)$rank
qr(climate)$rank

################################################################################

# participation ratios
# 16 sep 2025

participation_ratio <- function(pca) {
  eig <- pca$sdev^2
  sum(eig)^2 / sum(eig^2)
}

pca <- prcomp(siteDetections[, c(fxnVars, phenoVars)], scale. = TRUE)
participation_ratio(pca)


EffDim <- c(
  "climate" = participation_ratio(prcomp(siteDetections[, climVars], scale. = TRUE)),
  "traits" = participation_ratio(prcomp(siteDetections[, traitVars], scale. = TRUE)),
  "climate&traits" = participation_ratio(prcomp(siteDetections[, climVars], scale. = TRUE)) +
    participation_ratio(prcomp(siteDetections[, traitVars], scale. = TRUE)) -
    participation_ratio(prcomp(siteDetections[, c(climVars, traitVars)], scale. = TRUE))
)

climate <- participation_ratio(prcomp(siteDetections[, climVars], scale. = TRUE))
traits <- participation_ratio(prcomp(siteDetections[, traitVars], scale. = TRUE))
anthro <- participation_ratio(prcomp(siteDetections[, anthrVars], scale. = TRUE))

EffDim <- c(
  "climate" = climate,
  "traits" = traits,
  "climate&traits" = climate + traits - participation_ratio(prcomp(siteDetections[, c(climVars, traitVars)], scale. = TRUE)),
)

EffDim <- c(
  "climate" = climate,
  "anthro" = anthro,
  "climate&anthro" = climate + anthro - participation_ratio(prcomp(siteDetections[, c(climVars, anthrVars)], scale. = TRUE))
)

venn_diagram <- euler(EffDim, shape = "ellipse", input="union")
venn_diagram
plot(venn_diagram, quantities = TRUE)


################################################################################

# 18 sep 2025
# lets calculate all the participation ratios and put them in a table

participation_ratio <- function(pca) {
  eig <- pca$sdev^2
  sum(eig)^2 / sum(eig^2)
}

var_groups <- list(
  traits = traitVars,
  anthro = anthrVars,
  climate = climVars,
  fxn = fxnVars,
  pheno = phenoVars,
  topo = topoVars
)

pca <- prcomp(siteDetections[, climVars], scale. = TRUE)
participation_ratio(pca)

pr_df <- map_df(names(var_groups), function(g) {
  pca <- prcomp(siteDetections[, var_groups[[g]]], scale. = TRUE)
  tibble(group = g, participation_ratio = participation_ratio(pca))
})
