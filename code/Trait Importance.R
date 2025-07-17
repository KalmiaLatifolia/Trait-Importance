
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

################################################################################
# END OF PRIMER
################################################################################

# run random forest, make some test plots --------------------------------------

selectedSpecies <- species[8]

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


# loop it  ---------------------------------------------------------------------

# Prepare dataframe
PercentVarExplained <- data.frame(Species = character(),
                      spatVars = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each species
for (sp in species) {
  rf_formula <- as.formula(paste(sp, "~", paste(spatVars, collapse = " + ")))
  # Run random forest
  model <- randomForest(rf_formula, data = siteDetections, importance = FALSE)
  # Get final % variance explained
  percent_var_explained <- model$rsq[length(model$rsq)] * 100
  # Store result
  PercentVarExplained <- rbind(PercentVarExplained,
                   data.frame(Species = sp, spatVars = round(percent_var_explained, 2)))
}


ggplot(results, aes(x = PercentVarExplained)) +
  geom_histogram()



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
get_rf_rsq <- function(response, predictors) {
  response_vals <- binaryDetections[[response]]
  if (length(unique(response_vals)) < 2) return(NA) # escape for species always present or absent
  rf_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  model <- randomForest(rf_formula, data = binaryDetections, importance = FALSE)
  1 - model$err.rate[nrow(model$err.rate), "OOB"]
}


# First: all spatial variables
spat_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, spatVars))
PercentVarExplained_binary$spatVars <- spat_rsqs

# Then: trait variables
trait_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, traitVars))
PercentVarExplained_binary$traitVars <- round(trait_rsqs, 2)

# Anthr variables
anthr_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, anthrVars))
PercentVarExplained_binary$anthrVars <- round(anthr_rsqs, 2)

# climVars
clim_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, climVars))
PercentVarExplained_binary$climVars <- round(clim_rsqs, 2)

# fxnVars
fxn_rsqs <- sapply(species, function(sp) get_rf_rsq(sp, fxnVars))
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

















