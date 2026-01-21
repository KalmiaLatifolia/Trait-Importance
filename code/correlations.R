
# correlations
# 18 December 2025

library(ggplot2)
library(tidyverse)
library(paletteer)
library(ggcorrplot)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read full bird data ----------------------------------------------------------

siteDetections <- readRDS("data/siteDetections_foliarTraits_BioCube_20250522.rds")

# remove bad vars
siteDetections <- siteDetections %>% select(-gsp, -GHM_bu, -GHM_en, -NTL_2020)

# variable groups --------------------------------------------------------------

species = siteDetections[5:96]
spatVars  = siteDetections[97:178]

traits  = siteDetections[97:108]
anthr = siteDetections[109:120]
climate = siteDetections[121:147]
fxn = siteDetections[c(148:150, 175:178)]
pheno = siteDetections[151:160]
topo = siteDetections[161:174]

traitsnFxn = siteDetections[c(97:108, 148:150, 175:178)]

# do a correlation -------------------------------------------------------------

sp <- names(species[1])
trait <- names(traits[1])

ggplot(siteDetections, aes(x=get(trait), y=NFPD(get(sp)))) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_minimal()

# spearman correlation for non-linear monotonic relationships
spTraitCors <- as.data.frame(cor(NFPD(species), spatVars, method = "spearman"))

# tidy variable labels

var_labels <- read.csv("data/var_labels.csv")

# publication worthy correlation matrix ----------------------------------------

# pivot longer
cors_df <- as.data.frame(spTraitCors) %>%
  rownames_to_column("species") %>%           
  pivot_longer(-species, names_to = "Variable", values_to = "rho") 

# add variable labels
cors_df <- cors_df %>%
  left_join(var_labels, by = "Variable") %>%
  mutate(Label = factor(Label, levels = var_labels$Label))  # preserve desired axis order

# hclust rows and columns
row_ord <- hclust(as.dist(1 - cor(t(spTraitCors), method = "spearman")), 
                  method = "average")$order

col_ord <- hclust(as.dist(1 - cor(spTraitCors, method = "spearman")),
                  method = "average")$order

# plot it
ggplot(cors_df,
       aes(x = factor(Variable, levels = colnames(spTraitCors)[col_ord]),
           y = factor(species, levels = rownames(spTraitCors)[row_ord]),
           fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-0.8, 0.8), midpoint = 0, name = "Spearman’s ρ") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1,
                                   colour = var_labels$lab_color[match(colnames(spTraitCors)[col_ord], var_labels$Variable)]),
        axis.text.y = element_text(size = 6)) +
  xlab("") +
  ylab("") +
  scale_x_discrete(labels = var_labels$Label[match(colnames(spTraitCors)[col_ord], var_labels$Variable)])

