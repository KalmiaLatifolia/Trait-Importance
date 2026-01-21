
# merging TRY and SHIFT data
# 16 Jan 2026
# Laura Berman

library(rtry)
library(tidyverse)

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read in TRY data -------------------------------------------------------------

TRYdata <- rtry_import("data/TRY_data_20251208/45681.txt")

# remove BAAD data
TRYdata <- subset(TRYdata, TRYdata$DatasetID!=431) # BAAD data

# read in SHIFT data -----------------------------------------------------------

SHIFTdata <- read.csv("data/shift_sample_traits_cleaned4Laura.csv")

# merge ------------------------------------------------------------------------

# trim the try data
TRYdata <- subset(TRYdata, !is.na(TRYdata$TraitID))
TRYdata <- TRYdata[, c("AccSpeciesName", "TraitID", "TraitName", "StdValue", "UnitName")]

# trim the shift data
SHIFTdata <- SHIFTdata[, c("species_or_type", "LMA", "nitrogen", "phosphorus", "potassium", "sulfur", "lignin", "phenolics")]

# reformat SHIFT to match
SHIFT_long <- pivot_longer(
  SHIFTdata,
  cols = -species_or_type,
  names_to = "TraitName",
  values_to = "StdValue"
)

SHIFT_long$UnitName <- "mg/g"

# match trait names
unique(TRYdata$TraitID)
unique(TRYdata$TraitName)
TRYdata[TRYdata$TraitID==14]$TraitName <- "nitrogen"
TRYdata[TRYdata$TraitID==15]$TraitName <- "phosphorus"
TRYdata[TRYdata$TraitID==18]$TraitName
TRYdata[TRYdata$TraitID==44]$TraitName <- "potassium"
TRYdata[TRYdata$TraitID==87]$TraitName <- "lignin"
TRYdata[TRYdata$TraitID==92]$TraitName <- "cellulose"
TRYdata[TRYdata$TraitID==262]$TraitName <- "sulfur"
TRYdata[TRYdata$TraitID==773]$TraitName <- "canopy height"
TRYdata[TRYdata$TraitID==3106]$TraitName
TRYdata[TRYdata$TraitID==3107]$TraitName
TRYdata[TRYdata$TraitID==3117]$TraitName <- "SLA"

# convert SLA to LMA
TRY_SLA <- subset(TRYdata, TRYdata$TraitID==3117)
TRY_SLA$LMA_value <- 1/TRY_SLA$StdValue * 1000
TRY_LMA <- TRY_SLA
TRY_LMA$StdValue <- NULL
TRY_LMA$TraitName <- "LMA"
TRY_LMA$UnitName <- "g/m"
names(TRY_LMA)[names(TRY_LMA) == 'LMA_value'] <- 'StdValue'

# merge LMA
TRYdata <- rbind(TRYdata, TRY_LMA)

# match species names
SHIFT_long$AccSpeciesName <- sub("^((\\S+\\s+){1}\\S+).*", "\\1", SHIFT_long$species_or_type)

# remove unmatched columns
TRYdata$TraitID <- NULL
SHIFT_long$species_or_type <- NULL

# merge
TraitData <- rbind(TRYdata, SHIFT_long)

# combine subspecies
TraitData$AccSpeciesName[grepl("^Artemisia tridentata", TraitData$AccSpeciesName)] <- "Artemisia tridentata"
TraitData$AccSpeciesName[grepl("^Pinus ponderosa", TraitData$AccSpeciesName)] <- "Pinus ponderosa"
TraitData$AccSpeciesName[grepl("^Populus balsamifera", TraitData$AccSpeciesName)] <- "Populus balsamifera"
TraitData$AccSpeciesName[grepl("^Pinus contorta", TraitData$AccSpeciesName)] <- "Pinus contorta"
TraitData$AccSpeciesName[grepl("^Alnus incana", TraitData$AccSpeciesName)] <- "Alnus incana"
TraitData$AccSpeciesName[grepl("^Ceanothus cuneatus", TraitData$AccSpeciesName)] <- "Ceanothus cuneatus"
TraitData$AccSpeciesName[grepl("^Juniperus occidentalis", TraitData$AccSpeciesName)] <- "Juniperus occidentalis"


# make some nice plots ---------------------------------------------------------

# plant height
temp <- subset(TraitData, TraitData$TraitName== "Plant height vegetative")
temp <- subset(temp, temp$AccSpeciesName=="Pinus jeffreyi" |
                 temp$AccSpeciesName=="Pinus contorta" |
                 temp$AccSpeciesName=="Pinus ponderosa" |
                 temp$AccSpeciesName=="Pinus monticola" |
                 temp$AccSpeciesName=="Abies concolor" |
                 temp$AccSpeciesName=="Quercus chrysolepis" |
                 temp$AccSpeciesName=="Quercus douglasii" |
                 temp$AccSpeciesName=="Salix scouleriana" |
                 temp$AccSpeciesName=="Quercus kelloggii" |
                 temp$AccSpeciesName=="Populus fremontii" |
                 temp$AccSpeciesName=="Salix exigua" |
                 temp$AccSpeciesName=="Populus tremuloides" |
                 temp$AccSpeciesName=="Alnus incana" )
ggplot(temp, aes(x=StdValue, y=reorder(AccSpeciesName, StdValue, FUN = mean))) +
  geom_boxplot(color="#FDC71B", fill="#FFE090") + 
  xlab("Plant Height (m)") +
  xlim(0,80) +
  ylab("") + 
  theme_minimal() +
  scale_y_discrete(position = "right") +
  geom_text(data = temp %>% group_by(AccSpeciesName) %>% summarise(n = n()), 
            aes(x = 80 * 0.95, y = AccSpeciesName, label = paste0("n=", n)), 
            vjust = 0, size = 3)
ggsave("figures/PlantHeight.PDF", height=4.5, width=5)


temp <- subset(TraitData, TraitData$TraitName== "LMA")
temp <- subset(temp, temp$AccSpeciesName=="Pinus jeffreyi" |
                 temp$AccSpeciesName=="Pinus contorta" |
                 temp$AccSpeciesName=="Pinus ponderosa" |
                 temp$AccSpeciesName=="Pinus monticola" |
                 temp$AccSpeciesName=="Abies concolor" |
                 temp$AccSpeciesName=="Quercus chrysolepis" |
                 temp$AccSpeciesName=="Quercus douglasii" |
                 temp$AccSpeciesName=="Salix scouleriana" |
                 temp$AccSpeciesName=="Quercus kelloggii" |
                 temp$AccSpeciesName=="Populus fremontii" |
                 temp$AccSpeciesName=="Salix exigua" |
                 temp$AccSpeciesName=="Populus tremuloides" |
                 temp$AccSpeciesName=="Alnus incana" )
ggplot(temp, aes(x=StdValue, y=reorder(AccSpeciesName, StdValue, FUN = function(x) mean(x, na.rm = TRUE)))) +
  geom_boxplot(color="#842B3B", fill="#B4777C") + 
  xlab("LMA") +
  xlim(0,400) +
  ylab("") + 
  theme_minimal() +
  scale_y_discrete(position = "right") +
  geom_text(data = temp %>% group_by(AccSpeciesName) %>% summarise(n = n()), 
            aes(x = 400 * 0.95, y = AccSpeciesName, label = paste0("n=", n)), 
            vjust = 0, size = 3) 
ggsave("figures/LMA.PDF", height=4.5, width=5)
