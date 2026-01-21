
# TRY data
# 4 Dec 2025

library(rtry)


# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read in TRY data -------------------------------------------------------------

TRYdata <- rtry_import("data/TRY_data_20251208/45681.txt")

# Boxplots by trait ------------------------------------------------------------

unique(TRYdata$TraitID)

temp <- subset(TRYdata, TRYdata$TraitID==3106)

temp$AccSpeciesName[grepl("^Artemisia tridentata", temp$AccSpeciesName)] <- "Artemisia tridentata"
temp$AccSpeciesName[grepl("^Pinus ponderosa", temp$AccSpeciesName)] <- "Pinus ponderosa"
temp$AccSpeciesName[grepl("^Populus balsamifera", temp$AccSpeciesName)] <- "Populus balsamifera"
temp$AccSpeciesName[grepl("^Pinus contorta", temp$AccSpeciesName)] <- "Pinus contorta"
temp$AccSpeciesName[grepl("^Alnus incana", temp$AccSpeciesName)] <- "Alnus incana"
temp$AccSpeciesName[grepl("^Ceanothus cuneatus", temp$AccSpeciesName)] <- "Ceanothus cuneatus"
temp$AccSpeciesName[grepl("^Juniperus occidentalis", temp$AccSpeciesName)] <- "Juniperus occidentalis"

temp <- subset(temp, temp$AccSpeciesName=="Pinus jeffreyi" |
                 temp$AccSpeciesName=="Pinus contorta" |
                 temp$AccSpeciesName=="Pinus ponderosa" |
                 temp$AccSpeciesName=="Tsuga mertensiana" |
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
  geom_text(data = temp %>% group_by(AccSpeciesName) %>% summarise(n = n()), 
            aes(x = 80 * 0.95, y = AccSpeciesName, label = paste0("n=", n)), 
            vjust = -0.2)
ggsave("TRY_CH.PDF", height=4, width=6)


temp <- subset(TRYdata, TRYdata$TraitID==3117)
temp <- subset(temp, temp$DatasetID!=431) # BAAD data

temp$LMA <- 1/(temp$StdValue) * 1000

temp <- subset(temp, temp$AccSpeciesName!="Pinus contorta var. latifolia" & 
                 temp$AccSpeciesName!="Pinus lambertiana" &
                 temp$AccSpeciesName!="Juniperus occidentalis" &
                 temp$AccSpeciesName!="Populus balsamifera")

ggplot(temp, aes(x=LMA, y=reorder(AccSpeciesName, LMA))) +
  geom_boxplot(color="#842B3B", fill="#aa384c") + 
  xlab("LMA (g m-2)") +
  ylab("") + 
  xlim(c(0,400)) +
  theme_minimal() +
  geom_text(data = temp %>% group_by(AccSpeciesName) %>% summarise(n = n()), 
            aes(x = 400 * 0.95, y = AccSpeciesName, label = paste0("n=", n)), 
            vjust = -0.2)
ggsave("TRY_LMA.PDF", height=4, width=6)
