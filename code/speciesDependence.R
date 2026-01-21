
# species dependence
# how many bird species can be modeled significantly more accurately when data is available?
# 7 October 2025

library(VSURF)
library(randomForest)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full bird data ----------------------------------------------------------

siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# name variable groupings ------------------------------------------------------
names(siteDetections)
species <- colnames(siteDetections)[5:96]
varGroups <- list(
  spatVars  = colnames(siteDetections)[97:182],
  notTraits = colnames(siteDetections)[109:182],
  notAnthr  = colnames(siteDetections)[c(97:108, 124:182)],
  notClim   = colnames(siteDetections)[c(97:123, 152:182)],
  notFxn    = colnames(siteDetections)[c(97:151, 155:178)],
  notPheno  = colnames(siteDetections)[c(97:154, 165:182)],
  notTopo   = colnames(siteDetections)[c(97:164, 179:182)]
)

# calculate best random forest models ------------------------------------------

# create dataframe
speciesDependence <- data.frame(species=character(), variableGroup=character(), R2=numeric(), iteration=integer())

# Loop (THIS TAKES ALMOST A WEEK TO RUN)
for (sp in species) {
  y <- siteDetections[[sp]]
  for (grp in names(varGroups)) {
    x <- siteDetections[, varGroups[[grp]], drop=FALSE]
    for (i in 1:10) {
      vs <- VSURF(x, y, parallel=TRUE)
      selVars <- names(x)[vs$varselect.pred]
      rf <- randomForest(x[, selVars, drop=FALSE], y)
      R2 <- 1 - mean((rf$y - rf$predicted)^2) / var(rf$y)
      speciesDependence <- rbind(speciesDependence, data.frame(species=sp, variableGroup=grp, R2=R2, iteration=i))
    }
  }
}

# SAVE IT ----------------------------------------------------------------------

write_csv(speciesDependence, "speciesDependence_20251013.csv")

# do traits/climate/etc significantly improve R2 values? -----------------------

excluded_groups <- c("notTraits", "notAnthr", "notClim", "notFxn", "notPheno", "notTopo")

comparison_table <- lapply(excluded_groups, function(g) {
  speciesDependence %>%
    filter(variableGroup %in% c("spatVars", g)) %>%
    group_by(species) %>%
    summarise(
      included_r2 = mean(R2[variableGroup=="spatVars"], na.rm=TRUE),
      included_sd = sd(R2[variableGroup=="spatVars"], na.rm=TRUE),
      excluded_r2 = mean(R2[variableGroup==g], na.rm=TRUE),
      excluded_sd = sd(R2[variableGroup==g], na.rm=TRUE),
      p_value = ifelse(length(R2[variableGroup=="spatVars"])>1 & length(R2[variableGroup==g])>1,
                       t.test(R2[variableGroup=="spatVars"], R2[variableGroup==g])$p.value,
                       NA_real_),
      sigImp = p_value < 0.05 & included_r2 > excluded_r2,
      variable_group = g,
      .groups="drop"
    )
}) %>% bind_rows() %>% select(species, variable_group, included_r2, included_sd, excluded_r2, excluded_sd, p_value, sigImp)

# make a figure ----------------------------------------------------------------

custom_colors <- c(
  notTraits = "#842A3B",
  notAnthr  = "#7C4585",
  notClim   = "#6FA4AF",
  notFxn    = "#FCC61D",
  notPheno  = "#D97D55",
  notTopo   = "#A8BBA3"
)

comparison_table %>%
  group_by(variable_group) %>%
  summarise(sig_count = sum(sigImp, na.rm=TRUE)) %>%
  ggplot(aes(x=reorder(variable_group, sig_count), y=sig_count, fill=variable_group)) +
  geom_col() +
  geom_text(aes(label=c("...disturbance", "...climate", "...ecosystem function", "...phenology","...topography","...leaf physiology")), hjust=1.2, color="white", size=4) +
  scale_fill_manual(values=custom_colors) +
  coord_flip() +
  labs(x="Variable group", y="Number of species (out of 92)") +
  theme_minimal() +
  ggtitle("The R2 of relative abundance models significantly improved with the inclusion of...") +
  theme(legend.position="none")




