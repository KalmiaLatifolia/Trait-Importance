
# Section 1 of Full Protocol - Trait importance manuscript
# Laura Berman
# 13 April 2026
# Data cleaning and processing for analysis


library(tidyverse)
library(vroom)
library(geosphere)
library(sf)
library(mapview)
library(terra)
library(readxl)


# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

################################################################################
# Table of Contents
################################################################################

# Data processing --------------------------------------------------------------
# 1) Load and format bioacoustic data - 37
# 2) Remove burned areas - 150
# 3) Remove sites outside study area - 192
# 4) Get foliar traits -216
# 5) Get biocube variables - 273
# 6) Exclude duplicate/no variance variables - 323
# 7) Exclude rarely detected species - 377
# 8) Exclude NA variables + sites - 390


################################################################################
# Load and format bioacoustic data
################################################################################

# load ARU locations -----------------------------------------------------------

ARUlocations <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/focal_aru_metadata.csv") # 2100 locations
ARUlocations <- ARUlocations %>% distinct(survey_year, cell_unit, box, .keep_all = TRUE) # 2073 unique locations

# load detections --------------------------------------------------------------

# index detection files
files <- list.files("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean", 
                    pattern="filtered_detections\\.csv$", full.names=TRUE)
#load
det <- vroom(files)

# det is 108,747,136 obs (correct)
# names(det) = "filename" "cell_unit" "year" "survey_date" "month" "day" "hour" "min" "second" "relative_sec" "species_code" "logit_score"
# still missing common names and exact times


# load common names ------------------------------------------------------------

commonNames <- read.csv("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean/species_codes_thresholds.csv") # 241 species


# load effort ------------------------------------------------------------------

# index effort files
files <- list.files("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/laura_berman_diurnal_birds/data-clean", 
                    pattern="filtered_effort\\.csv$", full.names=TRUE)
#load
effort <- vroom(files)

# effort is 371,673 obs (correct)
# effort includes rows for days with 0 effort
# nrow(subset(effort, effort$aru_hours==7)) = 66,349 obs = total full recording days (7 hours of recording per day)

# calculate unit effort
unit_effort <- effort %>%
  subset(aru_hours==7) %>% # sampling schedule records 7 hours a day   
  group_by(cell_unit) %>%
  dplyr::summarize(startDate = min(date),
                   endDate = max(date),
                   active_days = n())
# unit_effort = 895 sites
# min(unit_effort$active_days) = 10
# max(unit_effort$active_days) = 108

write_rds(unit_effort, "data/unit_effort.rds")

# summarize detections by site -------------------------------------------------

siteDetections <- det %>% 
  group_by(cell_unit, species_code) %>%
  dplyr::summarise(n_detections = n())
# siteDetections = 59,626 obs = one row per site x species

# combine detections, locations, effort, and common names ----------------------

siteDetections <- merge(siteDetections, commonNames) 
siteDetections <- merge(siteDetections, unit_effort) 

# check max dist between redeployments -----------------------------------------

max_dists_by_cell <- ARUlocations %>%
  group_by(cell_unit) %>%
  summarise(max_distance_m = max(distm(cbind(long, lat))))

hist(max_dists_by_cell$max_distance_m)
mean(max_dists_by_cell$max_distance_m)
# max distance between redeployments is 985 meters
# on average redeployments are within 27 m 

# combine years with the same cell_unit ----------------------------------------

# average coordinates between redeployments
ARUlocations <- ARUlocations %>% 
  group_by(cell_unit) %>%
  dplyr::summarise(long = mean(long),
                   lat = mean(lat))
# 901 unique cell units

# add coordinates to siteDetections --------------------------------------------

siteDetections <- merge(siteDetections, ARUlocations)

# remove rows where active days < 30 -------------------------------------------

siteDetections <- subset(siteDetections, siteDetections$active_days > 30) # 59,342

# pivot wider ------------------------------------------------------------------
# [I will want active_days back later]

# calculate detection rates
siteDetections$obs_per_day <- siteDetections$n_detections / siteDetections$active_days

# remove extra columns
siteDetections <- siteDetections[c("cell_unit", "common_name", "obs_per_day", "long", "lat")]

# pivot wider
siteDetections <- siteDetections %>%
  group_by(cell_unit) %>%
  pivot_wider(names_from = common_name,
              values_from = obs_per_day)
# 889 sites, 104 obs

# No detections (NAs) become 0s
siteDetections[is.na(siteDetections)] <- 0



################################################################################
# Remove burned areas
################################################################################


# get burn years ---------------------------------------------------------------

# load fire data 2024
# fire data comes from https://gis.data.ca.gov/datasets/CALFIRE-Forestry::california-fire-perimeters-all/explore 
fire <- read_sf(dsn = "/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Fire perimeters/California_Fire_Perimeters_2024") # 22,261 fires
mapview(fire)

# format fire data
names(fire)[names(fire) == 'YEAR_'] <- 'burnYear'
fire <- fire[c("burnYear", "geometry")]
fire$burnYear <- as.numeric(fire$burnYear)

# make sf of site locations
siteDetections <- st_as_sf(siteDetections, coords = c("long", "lat"), crs=4326, remove=FALSE)

# match CRS
siteDetections <- st_transform(siteDetections, crs = 3857)

# find site locations that burned
siteDetections <- st_join(siteDetections, fire, multiple="all")
# this creates duplicate rows for sites that burned multiple times. 1133 obs

# keep only the most recent burn for each point
siteDetections <- siteDetections %>%
  group_by(cell_unit) %>% 
  slice_max(burnYear) %>% 
  ungroup() # 890
siteDetections <- siteDetections %>% distinct(cell_unit, .keep_all = TRUE) # 889 obs (correct)

# remove points that burned ----------------------------------------------------

siteDetections <- subset(siteDetections, is.na(siteDetections$burnYear) | siteDetections$burnYear < 2018) # 649 obs

# tidy up
rm(fire)


################################################################################
# Remove sites outside study area
################################################################################

# get tahoe and yosemite boxes -------------------------------------------------

tahoe <- read_sf("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/tahoe box/tahoe_shp.shp")
yosemite <- read_sf("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/yosemite box/polygon_yosemite_box.shp")
mapview(tahoe) +
  mapview(yosemite)

# exclude ARUs not in boxes ----------------------------------------------------
mapview(siteDetections) + mapview(tahoe) + mapview(yosemite)
# not necessary

# save siteDetections ----------------------------------------------------------

write_rds(siteDetections, "data/siteDetections_20260311.rds")
x <- as.data.frame(siteDetections)
x$geometry <- NULL
write_csv(x, "data/siteDetections_20260311.csv")
# siteDetections <- readRDS("data/siteDetections_20260311.rds")


################################################################################
# Get foliar traits
################################################################################

# buffer site locations
sites.buffer <- st_buffer(siteDetections, dist = 150) # 649
# a 150 m radius buffer = about 9 hectares, avg home range size of small passerines
sites.buffer <- sites.buffer[c("cell_unit", "long", "lat")]

# make a fxn to merge tiffs and extract values
get_trait <- function(trait){
  base_path <- "/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/mosaics/"
  yos <- terra::rast(paste0(base_path,"yosemite_20180622_",trait,"_30m_polished.tif"))
  tah <- terra::rast(paste0(base_path,"tahoe_20180621_",trait,"_30m_polished.tif"))
  r <- terra::merge(yos,tah)
  x <- terra::extract(r$s6v3_mean, sites.buffer, method="simple", exact=TRUE, na.rm=TRUE, fun=mean)
  x <- cbind(sites.buffer, x["s6v3_mean"])
  names(x)[names(x)=="s6v3_mean"] <- trait
  x[[trait]][x[[trait]]==0] <- NA
  as.data.frame(x)
}

# list traits
traits <- c("Nitrogen","LMA", "Cellulose", "Calcium", "Chlorophylls", "Fiber", 
            "Lignin", "NSC", "Phenolics", "Phosphorus", "Potassium", "Starch",
            "Sugar", "Sulfur")  # extend to all 10 traits

# do all trait extractions
foliarTraits <- Reduce(function(a,b) merge(a,b), lapply(traits, get_trait))

# save foliar traits -----------------------------------------------------------

write_rds(foliarTraits, "data/foliarTraits_20260313.rds")
x <- as.data.frame(foliarTraits)
x$geometry <- NULL 
write_csv(x, "data/foliarTraits_20260313.csv")
#foliarTraits <- readRDS("data/foliarTraits_20260313.rds")

# merge foliar traits with siteDetections --------------------------------------

foliarTraits <- as.data.frame(foliarTraits)
foliarTraits$geometry <- NULL 
siteDetections_foliarTraits <- merge(siteDetections, foliarTraits) #649

# subset to remove NA rows (no foliar trait data) ------------------------------

siteDetections_foliarTraits <- siteDetections_foliarTraits[!is.na(siteDetections_foliarTraits$Nitrogen), ] #578

# save foliar traits x siteDetections ------------------------------------------

write_rds(siteDetections_foliarTraits, "data/siteDetections_foliarTraits_20260313.rds")
x <- as.data.frame(siteDetections_foliarTraits)
x$geometry <- NULL
write_csv(x, "data/siteDetections_foliarTraits_20260313.csv")
# siteDetections_foliarTraits <- readRDS("data/siteDetections_foliarTraits_20260313.rds")


################################################################################
# Get biocube variables
################################################################################

# update buffer -----------------------------------------------------------------

# buffer site locations (we've removed some sites)
sites.buffer <- st_buffer(siteDetections_foliarTraits, dist = 150) # 578
# a 150 m radius buffer = about 9 hectares, avg home range size of small passerines
sites.buffer <- sites.buffer[c("cell_unit", "long", "lat")]


# get biocube paths and files --------------------------------------------------

path_CA <- "/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/BioCube_CA_Layers"
CAbiocube_files <- list.files(path_CA, full.names=TRUE)

# function to extract one layer ------------------------------------------------

get_biocube <- function(file){
  r <- terra::rast(file) # load the biocube layer
  x <- terra::extract(r, sites.buffer, method="simple", exact=TRUE, na.rm=TRUE, fun=mean) # extract from the site buffer
  names(x)[2] <- tools::file_path_sans_ext(basename(file)) # ensure column name is descriptive
  x <- cbind(sites.buffer, x) # merge them
  x[c("ID","geometry","lat", "long")] <- NULL # remove extra rows that cause issues
  as.data.frame(x)
}

# run all extractions ----------------------------------------------------------

BioCube_vars <- Reduce(function(a,b) merge(a,b, by="cell_unit"), lapply(CAbiocube_files, get_biocube)) # 578 sites, 105 vars

# save it ----------------------------------------------------------------------

write_rds(BioCube_vars, "data/BioCube_vars_20260317.rds")
# BioCube_vars <- readRDS("data/BioCube_vars_20260317.rds")


# merge with detections and foliar traits --------------------------------------

siteDetections_foliarTraits_BioCube <- merge(siteDetections_foliarTraits, BioCube_vars)


# save -------------------------------------------------------------------------

write_rds(siteDetections_foliarTraits_BioCube, "data/siteDetections_foliarTraits_BioCube_20260313.rds")
write_csv(siteDetections_foliarTraits_BioCube, "data/siteDetections_foliarTraits_BioCube_20260313.csv")



################################################################################
# Exclude duplicate/no variance variables
################################################################################

# correlation matrix of variables ----------------------------------------------
cor_matrix <- cor(st_drop_geometry(siteDetections_foliarTraits_BioCube[-1]), use = "complete.obs")

# These variables have too many NA values in the cor matrix and should be dropped:
sum(siteDetections_foliarTraits_BioCube$`Energy production cover` != 0, na.rm = TRUE)
# Wild Turkey - only 13 non-zero values
# Tree Swallow - only 2 non-zero values
# Brown-headed Cowbird - only 15 non-zero values
# California Thrasher - only 2 non-zero values
# Ice cover - only zero values
# Min nightlight - only 14 non-zero values
# Built-up cover - only 3 non-zero values
# Energy production cover - 58 non-zero values (essentially binary)

# remove NA variables ----------------------------------------------------------
siteDetections_foliarTraits_BioCube$`Wild Turkey` <- NULL
siteDetections_foliarTraits_BioCube$`Tree Swallow` <- NULL
siteDetections_foliarTraits_BioCube$`Brown-headed Cowbird` <- NULL
siteDetections_foliarTraits_BioCube$`California Thrasher` <- NULL
siteDetections_foliarTraits_BioCube$CA_Anthropocene_Hoskins_et_al_LU2005_ICE <- NULL
siteDetections_foliarTraits_BioCube$CA_Anthropocene_Li_et_al_NTL_min <- NULL
siteDetections_foliarTraits_BioCube$CA_Anthropocene_Theobald_et_al_tGHM_bu <- NULL
siteDetections_foliarTraits_BioCube$CA_Anthropocene_Theobald_et_al_tGHM_en <- NULL


# these variables are perfectly identical --------------------------------------
which(abs(cor_matrix) == 1 & row(cor_matrix) != col(cor_matrix), arr.ind = TRUE)
# structural richness and functional richness are identical variables
# CA_Diversity_GEDI_FRic.tiff and CA_Function_GEDI_Cali_StructuralRichness_1km_20240120.tiff
# removing the older one from the file folder: CA_Function_GEDI_Cali_StructuralRichness_1km_20240120

# remove duplicate variables ---------------------------------------------------
siteDetections_foliarTraits_BioCube$CA_Function_GEDI_Cali_StructuralRichness_1km_20240120 <- NULL

# check for variables with outliers/artifacts ----------------------------------
df_num <- st_drop_geometry(siteDetections_foliarTraits_BioCube[-1])
boxplot(df_num, las=2)
m <- as.matrix(df_num)
idx <- which.max(m)
list(
  row = row(m)[idx],
  column = colnames(m)[col(m)[idx]],
  value = m[idx])
hist(siteDetections_foliarTraits_BioCube$CA_Climate_CHELSA_gsp_v2)

# remove outlier variables -----------------------------------------------------
siteDetections_foliarTraits_BioCube$CA_Anthropocene_NASA_DARTE <- NULL

# siteDetections_foliarTraits_BioCube has 578 obs of 246 vars

################################################################################
# Exclude rarely detected species
################################################################################

# count number of detection sites ----------------------------------------------
species = siteDetections_foliarTraits_BioCube[4:100]
x <- as.data.frame(sapply(species, function(x) sum(x != 0, na.rm = TRUE)))

# remove species detected at fewer than 50 sites -------------------------------
siteDetections_foliarTraits_BioCube$`Song Sparrow` <- NULL
siteDetections_foliarTraits_BioCube$`Vaux's Swift` <- NULL
siteDetections_foliarTraits_BioCube$`Yellow-breasted Chat` <- NULL

################################################################################
# Exclude NA variables + sites
################################################################################

# How many NAs in each column? -------------------------------------------------
x <- as.data.frame(colSums(is.na(siteDetections_foliarTraits_BioCube)))
# burnYear = 348
# Snow cover days = 315
# Snow water equivalent = 315
# Frost change frequency = 71
# Tree sp richness = 25

# remove variables with more than 20 NA values ---------------------------------
siteDetections_foliarTraits_BioCube$burnYear <- NULL
siteDetections_foliarTraits_BioCube$CA_Climate_CHELSA_scd_v2 <- NULL
siteDetections_foliarTraits_BioCube$CA_Climate_CHELSA_swe_v2 <- NULL
siteDetections_foliarTraits_BioCube$CA_Climate_CHELSA_fcf_v2 <- NULL
siteDetections_foliarTraits_BioCube$CA_Diversity_Liang_et_al_2022_Tree_SR <- NULL

# check how many rows still have NAs -------------------------------------------
x <- as.data.frame(rowSums(is.na(siteDetections_foliarTraits_BioCube)))

# Remove rows that still have NA values ----------------------------------------
siteDetections_foliarTraits_BioCube <- siteDetections_foliarTraits_BioCube[complete.cases(st_drop_geometry(siteDetections_foliarTraits_BioCube)), ]

# Fabian suggests excluding Kling layers ---------------------------------------
siteDetections_foliarTraits_BioCube$CA_Diversity_Kling_et_al_2018_EndemicSpecies_PD <- NULL
siteDetections_foliarTraits_BioCube$CA_Diversity_Kling_et_al_2018_Species_Diversity <- NULL
siteDetections_foliarTraits_BioCube$CA_Diversity_Kling_et_al_2018_Species_Endemism <- NULL

# siteDetections_foliarTraits_BioCube has 553 obs of 234 vars

# save it ----------------------------------------------------------------------
write_rds(siteDetections_foliarTraits_BioCube, "data/siteDetections_foliarTraits_BioCube_20260320.rds")
write_csv(siteDetections_foliarTraits_BioCube, "data/siteDetections_foliarTraits_BioCube_20260320.csv")
# siteDetections_foliarTraits_BioCube <- readRDS("data/siteDetections_foliarTraits_BioCube_20260320.rds")
