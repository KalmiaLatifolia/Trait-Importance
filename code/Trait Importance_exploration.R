
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



#---------------------------
# Step 0: Define variables
#---------------------------
response <- siteDetections$Yellow.rumped.Warbler
predictors <- siteDetections[ , spatVars]

#---------------------------
# Step 1: Full Random Forest
#---------------------------
set.seed(123)
rf.full <- randomForest(x = predictors, y = response, ntree = 1000, importance = TRUE)
print(rf.full)

# Out-of-bag R² (last value is the final)
full_r2 <- rf.full$rsq[length(rf.full$rsq)]
cat("Full RF OOB R²:", full_r2, "\n")

#---------------------------
# Step 2: Variable Selection with VSURF
#---------------------------
set.seed(123)
vsurf.res <- VSURF(x = predictors, y = response)

# Best variables for prediction
sel.vars <- vsurf.res$varselect.pred
cat("Variables kept:", length(sel.vars), "\n")

#---------------------------
# Step 3: Fit reduced RF
#---------------------------
predictors.sel <- predictors[ , sel.vars, drop=FALSE]
names(predictors.sel)

set.seed(123)
rf.sel <- randomForest(x = predictors.sel, y = response, ntree = 1000, importance = TRUE)
print(rf.sel)

sel_r2 <- rf.sel$rsq[length(rf.sel$rsq)]
cat("Reduced RF OOB R²:", sel_r2, "\n")

#---------------------------
# Step 4: Cross-validation with caret
#---------------------------
set.seed(123)
ctrl <- trainControl(method = "cv", number = 10)

cv.rf <- train(response ~ ., data = data.frame(response=response, predictors.sel),
               method = "rf", trControl = ctrl, ntree = 1000)

cat("Cross-validated R²:", max(cv.rf$results$Rsquared), "\n")
