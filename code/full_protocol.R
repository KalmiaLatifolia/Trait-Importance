
# Full Protocol - Trait importance manuscript
# Laura Berman
# 11 March 2026

# Full data processing and analysis pipeline.

library(tidyverse)
library(vroom)
library(geosphere)
library(sf)
library(mapview)
library(terra)
library(readxl)
library(basemaps)
library(ggspatial)

# set working directory

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

################################################################################
# Table of Contents
################################################################################
# 1) Load and format Bioacoustic data - 33
# 2) Remove burned areas - 146
# 3) Remove sites outside study area - 188
# 4) Get foliar traits -211
# 5) Get biocube variables - 266
# 6) Make a nice map - 342
# 7)
# 8)
# 9)

################################################################################
# Load and format Bioacoustic data
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
# Get biocube variables (will need to re-run from here once I get Fabian's metadata)
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


# tidy variable names ----------------------------------------------------------

# load tidy names
tidy <- read_excel("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TableS1_Biocube_var_description.xlsx",
                   range = cell_cols(1:3))

# Create a named vector for renaming
rename_map <- setNames(tidy$VariableName, tidy$FileName)

# Rename only columns in BioCube that are in tidy$file
names(BioCube_vars)[names(BioCube_vars) %in% names(rename_map)] <- rename_map[names(BioCube_vars)[names(BioCube_vars) %in% names(rename_map)]]


# merge with detections and foliar traits --------------------------------------

siteDetections_foliarTraits_BioCube <- merge(siteDetections_foliarTraits, BioCube_vars)


# save -------------------------------------------------------------------------

write_rds(siteDetections_foliarTraits_BioCube, "data/siteDetections_foliarTraits_BioCube_20260313.rds")
write_csv(siteDetections_foliarTraits_BioCube, "data/siteDetections_foliarTraits_BioCube_20260313.csv")


################################################################################
#  make a nice map
################################################################################

# update sites
sites <- siteDetections_foliarTraits_BioCube[c("cell_unit", "long", "lat")]

# Create bounding box
ext <- st_as_sfc(st_bbox(c(
  xmin = -121, xmax = -118.5,
  ymin = 36.5, ymax = 39.5
), crs = 4326))

# Transform to EPSG:3857 (Web Mercator)
ext <- st_transform(ext, 3857)

#plot it
ggplot() +
  basemap_gglayer(ext, map_service = "esri", map_type = "world_dark_gray_base") +
  scale_fill_identity() + 
  geom_sf(data = sites, color = "navy", size = 1, show.legend = FALSE, alpha=0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = st_bbox(ext_3857)[c("xmin", "xmax")],
           ylim = st_bbox(ext_3857)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")


ggplot() +
  basemap_gglayer(ext, map_service = "esri", map_type = "world_physical_map") +
  scale_fill_identity() + 
  geom_sf(data = sites, color = "navy", size = 1, show.legend = FALSE, alpha=0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = st_bbox(ext_3857)[c("xmin", "xmax")],
           ylim = st_bbox(ext_3857)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")

ggsave("point_map.pdf", height=7, width=5)


### LEFT OFF HERE - Monday MARCH 19TH 2026 

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



################################################################################
# check for colinearity
################################################################################

# 90% cutoff

cor_matrix <- cor(siteDetections_foliarTraits_BioCube[,97:182], use = "complete.obs")
corrplot(cor_matrix,
         method = "color",
         order = "hclust",         
         hclust.method = "complete", 
         addrect = 3,             
         tl.cex = 0.8)

high_corr_matrix <- cor_matrix
high_corr_matrix[abs(high_corr_matrix) < 0.95 | diag(ncol(high_corr_matrix)) == 1] <- 0  
corrplot(high_corr_matrix,
         method = "color",
         order = "hclust",
         hclust.method = "complete", 
         tl.cex = 0.8)

high_corr_pairs <- as.data.frame(as.table(cor_matrix)) %>%
  filter(Var1 != Var2, abs(Freq) > 0.9)

# in Cali, there are dry summers and wet winters, so some bioclim vars are duplicate.
# bio16 and bio19 are 0.999 identical. precip of wettest/coldest quarter. remove bio16(wettest)?
# bio09 and bio10 are 0.999 identical. temp of dryest/warmest quarter. remove bio09(dryest)?
# ggd5 and gdd10 are 0.999 identical. remove gdd5?
# bio08 and bio11 are 0.998 identical. temp of wettest/coldest quarter. remove bio11(coldest)?
# swe and ngd0 are 0.98 identical. remove ngd0
# swe and scd 0.98 identical

# All these are > 95% correlated with DTM
# bio01, bio05, bio06, bio08, bio10, gdd10, ngd10, ngd5, annualPET, PETDryQuart, PETseas, PETWarmQuart, thermInd

# remove
siteDetections_foliarTraits_BioCube$bio16 <- NULL
siteDetections_foliarTraits_BioCube$bio09 <- NULL
siteDetections_foliarTraits_BioCube$gdd5 <- NULL
siteDetections_foliarTraits_BioCube$bio11 <- NULL
siteDetections_foliarTraits_BioCube$bio01 <- NULL
siteDetections_foliarTraits_BioCube$bio05 <- NULL
siteDetections_foliarTraits_BioCube$bio06 <- NULL
siteDetections_foliarTraits_BioCube$bio08 <- NULL
siteDetections_foliarTraits_BioCube$bio10 <- NULL
siteDetections_foliarTraits_BioCube$gdd0 <- NULL
siteDetections_foliarTraits_BioCube$gdd10 <- NULL
siteDetections_foliarTraits_BioCube$ngd10 <- NULL
siteDetections_foliarTraits_BioCube$ngd5 <- NULL
siteDetections_foliarTraits_BioCube$annualPET <- NULL
siteDetections_foliarTraits_BioCube$PETDryQuart <- NULL
siteDetections_foliarTraits_BioCube$PETseas <- NULL
siteDetections_foliarTraits_BioCube$PETWarmQuart <- NULL
siteDetections_foliarTraits_BioCube$ngd0 <- NULL
siteDetections_foliarTraits_BioCube$thermInd <- NULL
siteDetections_foliarTraits_BioCube$TRI <- NULL
siteDetections_foliarTraits_BioCube$scd <- NULL
siteDetections_foliarTraits_BioCube$bio19 <- NULL
siteDetections_foliarTraits_BioCube$Starch <- NULL

# save updated version
write_rds(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.rds")
write_csv(siteDetections_foliarTraits_BioCube, "siteDetections_foliarTraits_BioCube_20250522.csv")


