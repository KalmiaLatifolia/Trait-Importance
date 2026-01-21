
# Inherent dimensionality
# 19 Sep 2025
# Laura Berman

# Building the biocube
# What is the inherent dimensionality of the dataset?
# How much unique new information exists in each data category?

library(raster)
library(terra)
library(mapview)
library(tidyverse)

# set working directory
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")


# Load Foliar Traits -----------------------------------------------------------

# Nitrogen
NitrogenYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Nitrogen_30m_polished.tif")
NitrogenTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Nitrogen_30m_polished.tif")
Nitrogen <- merge(NitrogenYosemite, NitrogenTahoe)
rm(NitrogenTahoe, NitrogenYosemite) 

# LMA
LMAYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_LMA_30m_polished.tif")
LMATahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_LMA_30m_polished.tif")
LMA <- merge(LMAYosemite, LMATahoe)
mapview(aggregate(LMA, fact = 10))
rm(LMATahoe, LMAYosemite)

# Cellulose
CelluloseYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Cellulose_30m_polished.tif")
CelluloseTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Cellulose_30m_polished.tif")
Cellulose <- merge(CelluloseYosemite, CelluloseTahoe)
rm(CelluloseTahoe, CelluloseYosemite)

# Calcium
CalciumYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Calcium_30m_polished.tif")
CalciumTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Calcium_30m_polished.tif")
Calcium <- merge(CalciumYosemite, CalciumTahoe)
rm(CalciumTahoe, CalciumYosemite)

# Chlorophylls
ChlorophyllsYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Chlorophylls_30m_polished.tif")
ChlorophyllsTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Chlorophylls_30m_polished.tif")
Chlorophylls <- merge(ChlorophyllsYosemite, ChlorophyllsTahoe)
rm(ChlorophyllsTahoe, ChlorophyllsYosemite)

# Fiber
FiberYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Fiber_30m_polished.tif")
FiberTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Fiber_30m_polished.tif")
Fiber <- merge(FiberYosemite, FiberTahoe)
rm(FiberTahoe, FiberYosemite)

# Lignin
LigninYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Lignin_30m_polished.tif")
LigninTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Lignin_30m_polished.tif")
Lignin <- merge(LigninYosemite, LigninTahoe)
rm(LigninTahoe, LigninYosemite)

# NSC
NSCYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_NSC_30m_polished.tif")
NSCTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_NSC_30m_polished.tif")
NSC <- merge(NSCYosemite, NSCTahoe)
rm(NSCTahoe, NSCYosemite)

# Phenolics
PhenolicsYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Phenolics_30m_polished.tif")
PhenolicsTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Phenolics_30m_polished.tif")
Phenolics <- merge(PhenolicsYosemite, PhenolicsTahoe)
rm(PhenolicsTahoe, PhenolicsYosemite)

# Phosphorus
PhosphorusYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Phosphorus_30m_polished.tif")
PhosphorusTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Phosphorus_30m_polished.tif")
Phosphorus <- merge(PhosphorusYosemite, PhosphorusTahoe)
rm(PhosphorusTahoe, PhosphorusYosemite)

# Potassium
PotassiumYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Potassium_30m_polished.tif")
PotassiumTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Potassium_30m_polished.tif")
Potassium <- merge(PotassiumYosemite, PotassiumTahoe)
rm(PotassiumTahoe, PotassiumYosemite)

# Starch
StarchYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Starch_30m_polished.tif")
StarchTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Starch_30m_polished.tif")
Starch <- merge(StarchYosemite, StarchTahoe)
rm(StarchTahoe, StarchYosemite)

# Sugar
SugarYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Sugar_30m_polished.tif")
SugarTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Sugar_30m_polished.tif")
Sugar <- merge(SugarYosemite, SugarTahoe)
rm(SugarTahoe, SugarYosemite)

# Sulfur
SulfurYosemite <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/yosemite_20180622_Sulfur_30m_polished.tif")
SulfurTahoe <- rast("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/tahoe_20180621_Sulfur_30m_polished.tif")
Sulfur <- merge(SulfurYosemite, SulfurTahoe)
rm(SulfurTahoe, SulfurYosemite)


# Make a trait stack -----------------------------------------------------------

# tidy names
nitrogen_only <- Nitrogen$s6v3_mean
names(nitrogen_only) <- "Nitrogen"
LMA_only <- LMA$s6v3_mean
names(LMA_only) <- "LMA"
Calcium_only <- Calcium$s6v3_mean
names(Calcium_only) <- "Calcium"
Cellulose_only <- Cellulose$s6v3_mean
names(Cellulose_only) <- "Cellulose"
Chlorophylls_only <- Chlorophylls$s6v3_mean
names(Chlorophylls_only) <- "Chlorophylls"
Fiber_only <- Fiber$s6v3_mean
names(Fiber_only) <- "Fiber"
Lignin_only <- Lignin$s6v3_mean
names(Lignin_only) <- "Lignin"
NSC_only <- NSC$s6v3_mean
names(NSC_only) <- "NSC"
Phenolics_only <- Phenolics$s6v3_mean
names(Phenolics_only) <- "Phenolics"
Phosphorus_only <- Phosphorus$s6v3_mean
names(Phosphorus_only) <- "Phosphorus"
Potassium_only <- Potassium$s6v3_mean
names(Potassium_only) <- "Potassium"
Starch_only <- Starch$s6v3_mean
names(Starch_only) <- "Starch"
Sugar_only <- Sugar$s6v3_mean
names(Sugar_only) <- "Sugar"
Sulfur_only <- Sulfur$s6v3_mean
names(Sulfur_only) <- "Sulfur"


#stack
trait_stack <- c(nitrogen_only, LMA_only, Calcium_only, Cellulose_only, Chlorophylls_only, 
       Fiber_only, Lignin_only, NSC_only, Phenolics_only, Phosphorus_only,
       Potassium_only, Starch_only, Sugar_only, Sulfur_only)

# save the trait stack
saveRDS(trait_stack, "trait_stack_20251007.RDS")

# Assemble the biocube ---------------------------------------------------------

# get file paths
CAbiocube_files <- list.files("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/California biocube layers",  
                              full.names = TRUE)

# Use the first raster as template
biocube_template <- rast(CAbiocube_files[1])

# Function to align each raster to template
align_to_template <- function(f, template) {
  r <- rast(f)
  
  # Ensure CRS matches
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    # If CRS differs, reproject first
    if (crs(r) != crs(template)) {
      r <- project(r, template)
    }
    # Then resample to template grid
    r <- resample(r, template)
  }
  return(r)
}

# Combine into one stack
biocube_stack <- rast(lapply(CAbiocube_files, align_to_template, template = biocube_template))

# tidy names
names(biocube_stack[[86]]) <- "CA_Soil_SoilGrids_bdod_0_200"
names(biocube_stack[[89]]) <- "CA_Soil_SoilGrids_clay_0_200"
names(biocube_stack[[90]]) <- "CA_Soil_SoilGrids_nitrogen_0_200"
names(biocube_stack[[93]]) <- "CA_Soil_SoilGrids_sand_0_200"

# Merge biocube and traits -----------------------------------------------------

# match crs and resolution (this step needs to be fixed) (updating terra fixed it)
trait_stack_projected <- project(trait_stack, biocube_stack, method = "bilinear", mask = FALSE)

# merge
full_stack <- c(trait_stack_projected, biocube_stack)

# save
saveRDS(full_stack, "Full_CA_Biocube_20251007.RDS")

# do any layers have false zeros? Looks fine.

# Make a mask: TRUE if all layers have non-NA values, FALSE otherwise
valid_mask <- app(full_stack, fun = function(x) all(!is.na(x)))

# Apply mask: set NA anywhere the mask is FALSE
full_stack_masked <- mask(full_stack, valid_mask, maskvalues = 0)

# trim to remove outer NA only pixels
full_stack_trimmed <- trim(full_stack_masked)

# convert to dataframe
full_data <- as.data.frame(full_stack_trimmed, xy=TRUE)


# Create subsets ---------------------------------------------------------------

names(full_data)
spatVars <- colnames(full_data)[3:120]
traitVars <- colnames(full_data)[3:16]
anthrVars <- colnames(full_data)[17:31]
climVars <- colnames(full_data)[32:80]
fxnVars <- colnames(full_data)[c(81:83, 117:120)]
phenoVars <- colnames(full_data)[84:101]
topoVars <- colnames(full_data)[102:116]
notTraits <- colnames(full_data)[17:120]
notAnthr <- colnames(full_data)[c(3:16, 32:120)]
notClim <- colnames(full_data)[c(3:31, 81:120)]
notFxn <- colnames(full_data)[c(3:80, 84:116)]
notPheno <- colnames(full_data)[c(3:83, 102:120)]
notTopo <- colnames(full_data)[c(3:101, 117:120)]

var_groups <- list(
  spatVars = spatVars,
  traitVars = traitVars,
  anthrVars = anthrVars,
  climVars = climVars,
  fxnVars = fxnVars,
  phenoVars = phenoVars,
  topoVars = topoVars,
  notTraits = notTraits,
  notAnthr = notAnthr,
  notClim = notClim,
  notFxn = notFxn,
  notPheno = notPheno,
  notTopo = notTopo
)

# make a function - how many PCs explain 95% of variation? ---------------------

n_pcs_95 <- function(pca, threshold = 0.95) {
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  cum_var <- cumsum(var_exp)
  which(cum_var >= threshold)[1]
}

# calculate inherent dimensionality of subsets ---------------------------------

inherent_dimensionality <- map_df(names(var_groups), function(g) {
  pca <- prcomp(full_data[, var_groups[[g]]], scale. = TRUE)
  tibble(
    group = g,
    n_variables = length(var_groups[[g]]),
    n_pcs_95 = n_pcs_95(pca, 0.95)
  )
})


# make a converging bar chart --------------------------------------------------

manual_positioning <- tibble(
  group = c("anthrVars", "climVars", "topoVars",  "phenoVars", "fxnVars",   "traitVars", "spatVars", 
            "notAnthr",  "notClim", "notTopo", "notPheno",  "notFxn", "notTraits"),
  description = c("disturbance", "climate", "terrain", "phenology", "function", "traits", "All variables", 
                  "All variables except disturbance", "All variables except climate", "All variables except terrain", 
                  "All variables except phenology", "All variables except function", "All variables except traits"),
  side = c("right","right","right","right","right","right",
           "full", "left", "left","left","left","left","left"),
  y = c(7,6,5,4,3,2,1,
        7,6,5,4,3,2) 
)

inherent_dimensionality <- merge(inherent_dimensionality, manual_positioning)

full_width <- inherent_dimensionality$n_pcs_95[inherent_dimensionality$group == "spatVars"]
  
inherent_dimensionality <- inherent_dimensionality %>%
  mutate(
    xmin = case_when(side == "left"  ~ 0,
                     side == "right" ~ full_width - n_pcs_95,
                     side == "full"  ~ 0),
    xmax = case_when(side == "left"  ~ n_pcs_95,
                     side == "right" ~ full_width,
                     side == "full"  ~ full_width),
    width = xmax - xmin,
    # vertical offsets
    y_offset = case_when(side == "left"  ~ +0.05,
                         side == "right" ~ -0.05,
                         TRUE ~ 0),
    ymin = y - 0.28 + y_offset,
    ymax = y + 0.28 + y_offset,
    # improved label logic
    label_x = if_else(width > 6, (xmin + xmax)/2, (xmin + xmax)/2),  # center small bars too
    hjust_val = 0.5,
    label_col = if_else(width > 36, "white", "black"),
    label = paste0(description, " (", n_pcs_95, ")")
  ) %>%
  mutate(fill_group = case_when(
    side == "left" ~ "left",
    side == "full" ~ "full",
    side == "right" ~ group  # unique by group for different colors
  )) %>%
  mutate(side = factor(side, levels = c("left", "right", "full"))) %>%
  arrange(side)

# save it ----------------------------------------------------------------------

write.csv(inherent_dimensionality, "inherent_dimensionality_20251007.csv")
inherent_dimensionality <- read.csv("inherent_dimensionality_20251007.csv")

# plot it ----------------------------------------------------------------------

inherent_dimensionality$label[inherent_dimensionality$label=="traits (5)"] <- "physiology (5)"
inherent_dimensionality$label[inherent_dimensionality$label=="All variables except traits (34)"] <- "All variables except leaf physiology (34)"
inherent_dimensionality$label[inherent_dimensionality$label=="All variables except function (33)"] <- "All variables except ecosystem function (33)"

ggplot(inherent_dimensionality) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_group),
            color = NA, alpha = 0.7) +
  geom_text(aes(x = label_x, y = y, label = label, hjust = hjust_val, color = label_col),
            vjust = 0.5, size = 3.6, fontface = "bold") +
  scale_x_continuous(limits = c(0, full_width), expand = c(0,0), breaks = seq(0, full_width, 5)) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = c(left = "grey75", full = "black", traitVars = "#842A3B", climVars = "#6FA4AF", topoVars = "#A8BBA3", phenoVars = "#D97D55", fxnVars = "#FCC61D", anthrVars = "#7C4585")) +
  scale_color_identity() +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  labs(x = "Inherent Dimensionality \n (number of principal components required to explain 95% of the variation)", y = NULL)








