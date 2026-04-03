
# Partial Protocol - Trait importance manuscript
# create correlation matrix (Figure 2)

# correlations
# created 18 December 2025
# updated 3 April 2026

library(ggplot2)
library(tidyverse)
library(paletteer)
library(ggcorrplot)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# load dataset -----------------------------------------------------------------

siteDetections_foliarTraits_BioCube <- readRDS("data/siteDetections_foliarTraits_BioCube_20260320.rds")

# set up variable groups -------------------------------------------------------
names(siteDetections_foliarTraits_BioCube)
species = siteDetections_foliarTraits_BioCube[4:97]
spatVars  = siteDetections_foliarTraits_BioCube[98:234]

# NFPD function ----------------------------------------------------------------
FPD <- function(species_col) {species_col / max(species_col)}
meanFPD <- function(species_col) {mean(FPD(species_col)[FPD(species_col)>0])}
NFPD <- function(species_col) {(FPD(species_col))^(log(0.5)/log(meanFPD(species_col)))}


# spearman correlation for non-linear monotonic relationships ------------------
spTraitCors <- as.data.frame(cor(NFPD(species), spatVars, method = "spearman"))

# pivot longer -----------------------------------------------------------------
cors_df <- as.data.frame(spTraitCors) %>%
  rownames_to_column("species") %>%           
  pivot_longer(-species, names_to = "Variable", values_to = "rho") 

# hclust rows and columns ------------------------------------------------------

row_hc <- hclust(as.dist(1 - cor(t(spTraitCors), method = "spearman")), 
                 method = "centroid")$order

col_hc <- hclust(as.dist(1 - cor(spTraitCors, method = "spearman")),
                 method = "centroid")$order

# Tidy variable labels ---------------------------------------------------------

# load tidy names
tidy <- read_excel("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TableS1_Biocube_var_description.xlsx",
                   range = cell_cols(1:4))

# add variable labels to correlation dataframe 
cors_df <- cors_df %>%
  left_join(tidy, by = "Variable") %>%
  mutate(Label = factor(Label, levels = tidy$Label))  # preserve desired axis order

# data exploration - look at specific cor values
cors_df %>%
  filter(species == "Mountain Bluebird", Label == "Canopy height") %>%
  pull(rho)

# plot it ----------------------------------------------------------------------
ggplot(cors_df,
       aes(x = factor(Variable, levels = colnames(spTraitCors)[col_hc]),
           y = factor(species, levels = rownames(spTraitCors)[row_hc]),
           fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-0.85, 0.8), midpoint = 0, name = "Spearman’s ρ") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1,
                                   colour = tidy$labelColor[match(colnames(spTraitCors)[col_hc], tidy$Variable)]),
        axis.text.y = element_text(size = 6)) +
  xlab("") +
  ylab("") +
  scale_x_discrete(labels = tidy$Label[match(colnames(spTraitCors)[col_hc], tidy$Variable)])

# save it
ggsave("figures/Correlation_Matrix1.pdf", height=10, width=15)



