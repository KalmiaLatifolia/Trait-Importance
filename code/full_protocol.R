
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
library(randomForest)
library(xgboost)
library(SHAPforxgboost)
library(future)
library(future.apply)
library(progressr)
library(ggpubr)

# set working directory

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

################################################################################
# Table of Contents
################################################################################

# Data processing --------------------------------------------------------------
# 1) Load and format bioacoustic data - 45
# 2) Remove burned areas - 158
# 3) Remove sites outside study area - 200
# 4) Get foliar traits -224
# 5) Get biocube variables - 281
# 6) Exclude duplicate/no variance variables - 331
# 7) Exclude rarely detected species - 385
# 8) Exclude NA variables + sites - 405

# Create Figures ---------------------------------------------------------------
# 9) Make a nice map (Figure 1) - 433
# 10) Drop geometry - 482
# 11) Make a correlation matrix plot (Figure 2) - 490
# 12) run XGBoost models - 565
# 13) which species R2 significantly improves with each category? (Figure 3) - 710
# 14) which species RMSE significantly improves with each category? (Figure 3 alt) - 852
# 15) Heat map - variable importance - 977
# 16) stacked bar chart (Figure 4) - 1021
# 17) Single species SHAP scores - 1089


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

# plot it - grey base
ggplot() +
  basemap_gglayer(ext, map_service = "esri", map_type = "world_dark_gray_base") +
  scale_fill_identity() + 
  geom_sf(data = sites, color = "navy", size = 1, show.legend = FALSE, alpha=0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = st_bbox(ext)[c("xmin", "xmax")],
           ylim = st_bbox(ext)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")

# plot it - terrain base
ggplot() +
  basemap_gglayer(ext, map_service = "esri", map_type = "world_physical_map") +
  scale_fill_identity() + 
  geom_sf(data = sites, color = "navy", size = 1, show.legend = FALSE, alpha=0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = st_bbox(ext)[c("xmin", "xmax")],
           ylim = st_bbox(ext)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")

ggsave("point_map.pdf", height=7, width=5)


################################################################################
# Drop geometry 
################################################################################

siteDetections_foliarTraits_BioCube <- st_drop_geometry(siteDetections_foliarTraits_BioCube)


################################################################################
# Make a correlation matrix plot (Figure 2)
################################################################################

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



################################################################################
# run XGBoost models
################################################################################

# name variable grouping -------------------------------------------------------
names(siteDetections_foliarTraits_BioCube)
species <- colnames(siteDetections_foliarTraits_BioCube)[4:97]
varSets <- list(
  spatVars  = tidy$Variable,
  notTraits = tidy$Variable[tidy$Category != "Traits"],
  notDist   = tidy$Variable[tidy$Category != "Disturbance"],
  notClim   = tidy$Variable[tidy$Category != "Climate"],
  notStr    = tidy$Variable[tidy$Category != "Structure"],
  notPheno  = tidy$Variable[tidy$Category != "Phenology"],
  notTerr   = tidy$Variable[tidy$Category != "Terrain"]
)

# make sure all column names in varSets match real variables
varSets <- lapply(varSets, function(x) x[x %in% colnames(siteDetections_foliarTraits_BioCube)])


# xgboost with SHAP and iterations in parallel ---------------------------------
# ------------------------------------------------------------------------------

# Use multiple cores
plan(multisession, workers = parallel::detectCores() - 1)  

# list parallelizable tasks (6580)
tasks <- expand.grid(
  varSet = names(varSets),
  species_i = seq_along(species),
  iter = 1:100
)

# make a fxn -------------------------------------------------------------------
run_task <- function(vs, i, iter) {
  set.seed(123 + iter)
  
  # choose variable set
  spatVars <- siteDetections_foliarTraits_BioCube[, varSets[[vs]], drop = FALSE]
  
  # choose species
  y <- NFPD(siteDetections_foliarTraits_BioCube[[species[i]]])
  
  # split into train and test sets
  n <- nrow(spatVars)
  train_idx <- sample(n, floor(0.8*n))
  test_idx  <- setdiff(seq_len(n), train_idx)
  
  X_train <- spatVars[train_idx, ]
  y_train <- y[train_idx]
  X_test  <- spatVars[test_idx, ]
  y_test  <- y[test_idx]
  
  dtrain <- xgb.DMatrix(as.matrix(X_train), label = y_train)
  dtest  <- xgb.DMatrix(as.matrix(X_test),  label = y_test)
  
  # set xgboost parameters
  params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)
  
  # run 5-fold cross validation (find best number of iterations)
  cv <- xgb.cv(params = params, data = dtrain,
               nrounds = 100, nfold = 5, metrics = "rmse",
               verbose = 0, early_stopping_rounds = 10)
  
  # run xgboost
  xgb_model <- xgb.train(params = params, data = dtrain, nrounds = cv$early_stop$best_iteration)
  
  # predict on test data
  pred <- predict(xgb_model, dtest)
  
  # calculate R2 of prediction vs test data
  R2_test  <- 1 - sum((y_test  - pred)^2) / sum((y_test  - mean(y_test))^2)
  RMSE_test <- sqrt(mean((pred - y_test)^2))
  R2_train <- 1 - sum((y_train - predict(xgb_model, dtrain))^2) / sum((y_train - mean(y_train))^2)
  RMSE_train <- sqrt(mean((predict(xgb_model, dtrain) - y_train)^2))
  
  # calculate variable importance
  imp <- xgb.importance(colnames(spatVars), model = xgb_model)
  imp$species   <- species[i]
  imp$varSet    <- vs
  imp$iteration <- iter
  
  # calculate SHAP importance
  shap <- shap.values(xgb_model, X_train = X_test)
  shap_importance <- data.frame(
    Feature = names(shap$mean_shap_score),
    SHAP_importance = shap$mean_shap_score
  )
  
  # keep importance data
  imp <- merge(imp, shap_importance)
  
  # keep model data
  model_row <- data.frame(
    species = species[i],
    varSet = names(varSets)[vs],
    iteration = iter,
    zeros = sum(y == 0, na.rm = TRUE),
    best_nrounds_cv = cv$early_stop$best_iteration,
    R2_test = R2_test,
    RMSE_test = RMSE_test,
    R2_train = R2_train,
    RMSE_train = RMSE_train
  )
  
  list(model = model_row, importance = imp)
}


# run the fxn -------------------------------------- (~30 min for 10 iterations)
handlers(global = TRUE) 

with_progress({
  p <- progressor(along = seq_len(nrow(tasks)))
  
  results <- future_lapply(
    seq_len(nrow(tasks)),
    function(k) {
      res <- run_task(
        vs   = tasks$varSet[k],
        i    = tasks$species_i[k],
        iter = tasks$iter[k]
      )
      p(sprintf("Finished task %d", k))
      res
    }
  )
})

# turn off parallel processing 
plan(sequential)

# combine outputs --------------------------------------------------------------
xgb_modelParameters <- do.call(rbind, lapply(results, `[[`, "model"))
xgb_variableImportance <- do.call(rbind, lapply(results, `[[`, "importance"))

# save it ----------------------------------------------------------------------
write_rds(xgb_modelParameters, "data/xgb_modelParameters_20260325.rds")
write_csv(xgb_modelParameters, "data/xgb_modelParameters_20260325.csv")

write_rds(xgb_variableImportance, "data/xgb_variableImportance_20260325.rds")
write_csv(xgb_variableImportance, "data/xgb_variableImportance_20260325.csv")
#xgb_variableImportance <- readRDS("data/xgb_variableImportance_20260325.rds")

################################################################################
# which species R2 significantly improves with each category? (Figure 3)
################################################################################

# here is the varsets chunk again as a reminder
varSets <- list(
  spatVars  = tidy$Variable,
  notTraits = tidy$Variable[tidy$Category != "Traits"],
  notDist   = tidy$Variable[tidy$Category != "Disturbance"],
  notClim   = tidy$Variable[tidy$Category != "Climate"],
  notStr    = tidy$Variable[tidy$Category != "Structure"],
  notPheno  = tidy$Variable[tidy$Category != "Phenology"],
  notTerr   = tidy$Variable[tidy$Category != "Terrain"]
)
varSets <- lapply(varSets, function(x) x[x %in% colnames(siteDetections_foliarTraits_BioCube)])

# run t.tests for each category
cat_ttest <- xgb_modelParameters %>%
  group_by(species) %>%
  summarise(notTraits = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notTraits"], alternative="greater")$p.value,
            notClim = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notClim"], alternative="greater")$p.value,
            notStr = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notStr"], alternative="greater")$p.value,
            notDist = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notDist"], alternative="greater")$p.value,
            notPheno = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notPheno"], alternative="greater")$p.value,
            notTerr = t.test(R2_test[varSet=="spatVars"], R2_test[varSet=="notTerr"], alternative="greater")$p.value,
            R2 = mean(R2_test))

# exclude species with R2 < 0 
cat_ttest <- cat_ttest[cat_ttest$R2 > 0, ]


# plot it

ssp <- cat_ttest$species[cat_ttest$notTraits < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notTraits")
temp$varSet[temp$varSet=="notTraits"] <- "Without Foliar Traits"
temp$varSet[temp$varSet=="spatVars"] <- "With Foliar Traits"

p1 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#842B3B", "black")) +
  scale_fill_manual(values = c("#aa384c", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model R2") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notStr < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notStr")
temp$varSet[temp$varSet=="notStr"] <- "Without Forest Structure"
temp$varSet[temp$varSet=="spatVars"] <- "With Forest Structure"

p2 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#FDC71B", "black")) +
  scale_fill_manual(values = c("#FFE090", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model R2") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notClim < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notClim")
temp$varSet[temp$varSet=="notClim"] <- "Without Climate"
temp$varSet[temp$varSet=="spatVars"] <- "With Climate"

p3 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#70A4AF", "black")) +
  scale_fill_manual(values = c("#91b9c1", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model R2") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notPheno < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notPheno")
temp$varSet[temp$varSet=="notPheno"] <- "Without Phenology"
temp$varSet[temp$varSet=="spatVars"] <- "With Phenology"

p4 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#D97C55", "black")) +
  scale_fill_manual(values = c("#e29c7f", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model R2") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notTerr < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notTerr")
temp$varSet[temp$varSet=="notTerr"] <- "Without Terrain"
temp$varSet[temp$varSet=="spatVars"] <- "With Terrain"

p5 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#A8BBA3", "black")) +
  scale_fill_manual(values = c("#c4d1c0", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model R2") +
  ylab("")

ssp <- cat_ttest$species[cat_ttest$notDist < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notDist")
temp$varSet[temp$varSet=="notDist"] <- "Without Disturbance"
temp$varSet[temp$varSet=="spatVars"] <- "With Disturbance"

p6 <- ggplot(temp, aes(x=R2_test, y=reorder(species, R2_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#7C4584", "black")) +
  scale_fill_manual(values = c("#AC86B0", "grey30")) +
  xlim(-0.25,0.9) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model R2") +
  ylab("")

annotate_figure(
  ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, heights = c(4, 1), align="hv"),
  left = text_grob("Species with significant improvement", rot = 90, size = 14, vjust = 1))

ggsave("SpeciesBestR2_20260326.pdf", height=10, width=15)




################################################################################
# which species RMSE significantly improves with each category? (Figure 3 alt)
################################################################################

# run t.tests for each category
cat_ttest <- xgb_modelParameters %>%
  group_by(species) %>%
  summarise(notTraits = t.test(RMSE_test[varSet=="spatVars"], RMSE_test[varSet=="notTraits"], alternative="less")$p.value,
            notClim = t.test(RMSE_test[varSet=="spatVars"], RMSE_test[varSet=="notClim"], alternative="less")$p.value,
            notStr = t.test(RMSE_test[varSet=="spatVars"], RMSE_test[varSet=="notStr"], alternative="less")$p.value,
            notDist = t.test(RMSE_test[varSet=="spatVars"], RMSE_test[varSet=="notDist"], alternative="less")$p.value,
            notPheno = t.test(RMSE_test[varSet=="spatVars"], RMSE_test[varSet=="notPheno"], alternative="less")$p.value,
            notTerr = t.test(RMSE_test[varSet=="spatVars"], RMSE_test[varSet=="notTerr"], alternative="less")$p.value,
            RMSE = mean(RMSE_test))

# plot it

ssp <- cat_ttest$species[cat_ttest$notTraits < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notTraits")
temp$varSet[temp$varSet=="notTraits"] <- "Without Foliar Traits"
temp$varSet[temp$varSet=="spatVars"] <- "With Foliar Traits"

p1 <- ggplot(temp, aes(x=RMSE_test, y=reorder(species, RMSE_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#842B3B", "black")) +
  scale_fill_manual(values = c("#aa384c", "grey30")) +
  xlim(0.05, 0.25) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model RMSE") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notStr < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notStr")
temp$varSet[temp$varSet=="notStr"] <- "Without Forest Structure"
temp$varSet[temp$varSet=="spatVars"] <- "With Forest Structure"

p2 <- ggplot(temp, aes(x=RMSE_test, y=reorder(species, RMSE_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#FDC71B", "black")) +
  scale_fill_manual(values = c("#FFE090", "grey30")) +
  xlim(0.05, 0.25) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model RMSE") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notClim < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notClim")
temp$varSet[temp$varSet=="notClim"] <- "Without Climate"
temp$varSet[temp$varSet=="spatVars"] <- "With Climate"

p3 <- ggplot(temp, aes(x=RMSE_test, y=reorder(species, RMSE_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#70A4AF", "black")) +
  scale_fill_manual(values = c("#91b9c1", "grey30")) +
  xlim(0.05, 0.25) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model RMSE") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notPheno < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notPheno")
temp$varSet[temp$varSet=="notPheno"] <- "Without Phenology"
temp$varSet[temp$varSet=="spatVars"] <- "With Phenology"

p4 <- ggplot(temp, aes(x=RMSE_test, y=reorder(species, RMSE_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#D97C55", "black")) +
  scale_fill_manual(values = c("#e29c7f", "grey30")) +
  xlim(0.05, 0.25) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model RMSE") +
  ylab("")


ssp <- cat_ttest$species[cat_ttest$notTerr < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notTerr")
temp$varSet[temp$varSet=="notTerr"] <- "Without Terrain"
temp$varSet[temp$varSet=="spatVars"] <- "With Terrain"

p5 <- ggplot(temp, aes(x=RMSE_test, y=reorder(species, RMSE_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#A8BBA3", "black")) +
  scale_fill_manual(values = c("#c4d1c0", "grey30")) +
  xlim(0.05, 0.25) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model RMSE") +
  ylab("")

ssp <- cat_ttest$species[cat_ttest$notDist < 0.05]
temp <- subset(xgb_modelParameters, species %in% ssp)
temp <- subset(temp, temp$varSet=="spatVars" | temp$varSet=="notDist")
temp$varSet[temp$varSet=="notDist"] <- "Without Disturbance"
temp$varSet[temp$varSet=="spatVars"] <- "With Disturbance"

p6 <- ggplot(temp, aes(x=RMSE_test, y=reorder(species, RMSE_test), color=varSet, fill=varSet)) +
  geom_boxplot() +
  scale_color_manual(values = c("#7C4584", "black")) +
  scale_fill_manual(values = c("#AC86B0", "grey30")) +
  xlim(0.05, 0.25) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position = "bottom") +
  xlab("Best Model RMSE") +
  ylab("")

annotate_figure(
  ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, heights = c(4, 1), align="hv"),
  left = text_grob("Species with significant improvement", rot = 90, size = 14, vjust = 1))

ggsave("SpeciesBestRMSE_20260326.pdf", height=10, width=15)



################################################################################
# Heat map - variable importance
################################################################################


# average across iterations ----------------------------------------------------

temp <- xgb_variableImportance %>%
  subset(varSet == "spatVars") %>%
  group_by(Feature, species) %>%
  summarise(SHAP_importance = mean(SHAP_importance, na.rm=TRUE), .groups="drop") %>%
  rename(Variable = Feature) %>%
  left_join(tidy, by="Variable")

# order y (variables) by total importance --------------------------------------
var_order <- temp %>% 
  group_by(Label) %>% 
  summarise(total = sum(SHAP_importance, na.rm=TRUE), .groups="drop") %>% 
  arrange(total) %>% 
  pull(Label)

# order x (species) by total importance ----------------------------------------
sp_order <- temp %>% 
  group_by(species) %>% 
  summarise(total = sum(SHAP_importance, na.rm=TRUE), .groups="drop") %>% 
  arrange(desc(total)) %>% 
  pull(species)

# plot it ----------------------------------------------------------------------

ggplot(temp, aes(x = factor(species, levels = sp_order), 
                 y = factor(Label, levels = var_order), 
                 fill = SHAP_importance)) +
  geom_tile() +
  scale_fill_continuous(palette = c("white", "#FC9272", "#DE2D26"), na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=6, angle=45, hjust=1),
        axis.text.y = element_text(size = 6, colour = tidy$labelColor[match(var_order, tidy$Label)])) +
  ylab("Variable") +
  xlab("Species")

ggsave("SHAPimportance_20260329.pdf", height=11, width=11)


################################################################################
# stacked bar chart - Figure 4
################################################################################

# format dataset ---------------------------------------------------------------
temp <- xgb_variableImportance %>%
  subset(varSet == "spatVars") %>% # only full models
  group_by(Feature, species) %>%
  summarise(SHAP_importance = mean(SHAP_importance, na.rm=TRUE), .groups="drop") %>%
  rename(Variable = Feature) %>%
  left_join(tidy, by="Variable") %>%
  group_by(species) %>%
  mutate(total = sum(SHAP_importance, na.rm = TRUE), # cumulative SHAP values for each model
         prop_T = sum(SHAP_importance[Category == "Traits"], na.rm = TRUE) / total, # % attributed to Traits
         prop_bar = SHAP_importance / total) %>% # % attributed to each var
  ungroup() %>%
  mutate(
    species = reorder(species, prop_T), # order by trait importance
    Category = factor(Category, levels = c("Climate", "Disturbance", "Terrain", "Phenology", "Structure", "Traits"))) %>% # preferred Category order
  arrange(species, Category)

# what is the top variable per species?
temp1 <- temp %>%
  group_by(species) %>%
  slice_max(order_by = prop_bar, n = 1, with_ties = FALSE) %>%
  ungroup()

# create panels ----------------------------------------------------------------

p1 <- temp %>%
  ggplot(aes(y = species, x = SHAP_importance, fill = Category)) +
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes( label = ifelse(prop_bar > 0.1, Label, ""), x = SHAP_importance / 2), 
            position = position_fill(vjust = 0.5), size = 2) +
  scale_fill_manual(values = c("#70A4AF", "#7C4584", "#A8BBA3", "#D97C55", "#FDC71B", "#842B3B")) +
  theme_minimal() +
  xlab("SHAP importance (%)") +
  ylab("Species")


p2 <- temp %>%
  group_by(Category, species) %>%
  summarise(cat_size = sum(prop_bar)) %>%
  group_by(Category) %>%
  summarise(mean_cat_prop = mean(cat_size), y="Average") %>%
  ggplot(aes(y=y, x = mean_cat_prop, fill = Category)) +
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = paste(round(mean_cat_prop, digits = 2)*100, "%"), x = mean_cat_prop), position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c("#70A4AF", "#7C4584", "#A8BBA3", "#D97C55", "#FDC71B", "#842B3B")) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Average SHAP importance (%)") +
  ylab("")

# plot it ----------------------------------------------------------------------

ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(9, 1), align="hv")

# save it ----------------------------------------------------------------------

ggsave("SHAPimportance_20260330.png", height=12, width=12)



################################################################################
# Single species SHAP scores
################################################################################

# get SHAP scores for 1 species ------------------------------------------------

# choose species (25, 70, 32, 29, 31, 24, 61, 48)
species <- colnames(siteDetections_foliarTraits_BioCube)[4:97]
species
i <- 26
y <- NFPD(siteDetections_foliarTraits_BioCube[[species[i]]])

# use full variable set
vs <- 1
spatVars <- siteDetections_foliarTraits_BioCube[, varSets[[vs]], drop = FALSE]

# split into train and test sets
n <- nrow(spatVars)
train_idx <- sample(n, floor(0.8*n))
test_idx  <- setdiff(seq_len(n), train_idx)

X_train <- spatVars[train_idx, ]
y_train <- y[train_idx]
X_test  <- spatVars[test_idx, ]
y_test  <- y[test_idx]

dtrain <- xgb.DMatrix(as.matrix(X_train), label = y_train)
dtest  <- xgb.DMatrix(as.matrix(X_test),  label = y_test)

# set xgboost parameters
params <- list(objective = "reg:squarederror", max_depth = 6, eta = 0.1)

# run 5-fold cross validation (find best number of iterations)
cv <- xgb.cv(params = params, data = dtrain,
             nrounds = 100, nfold = 5, metrics = "rmse",
             verbose = 0, early_stopping_rounds = 10)

# run xgboost
xgb_model <- xgb.train(params = params, data = dtrain, nrounds = cv$early_stop$best_iteration)

# calculate SHAP scores
shap <- shap.values(xgb_model, X_train = X_test)
shap_values <- shap$shap_score

# build dataframe
temp <- data.frame(species = rep(species[i], nrow(shap_values)),  
                   shap_CanopyHeight = shap_values$CA_Function_Favrichon_Sierra_CanopyHeight_30m_20250508, 
                   CanopyHeight = X_test$CA_Function_Favrichon_Sierra_CanopyHeight_30m_20250508,
                   shap_LMA = shap_values$LMA, LMA = X_test$LMA,
                   shap_Phenolics = shap_values$Phenolics, Phenolics = X_test$Phenolics,
                   shap_Potassium = shap_values$Potassium, Potassium = X_test$Potassium)


# plot it ----------------------------------------------------------------------

# canopy height
ggplot() +
  geom_point(data = siteDetections_foliarTraits_BioCube, aes(x=CA_Function_Favrichon_Sierra_CanopyHeight_30m_20250508/100, y=NFPD(get(species[i]))), color="#FDC71B") +
  geom_smooth(data = temp, aes(x=CanopyHeight/100, y=((shap_CanopyHeight *4) +0.4)), color="black") + 
  scale_y_continuous(name=paste(species[i], "\nnormalized detection rate"), sec.axis = sec_axis(~ . /4 -0.4/4, name="SHAP value")) +
  theme_minimal() +
  xlab("Canopy Height") +
  xlim(c(0,80)) +
  theme(axis.title.y.left=element_text(color="#FDC71B"))
ggsave("figures/LGold_CH.PDF", width=6, height=3)

# LMA
ggplot() +
  geom_point(data = siteDetections_foliarTraits_BioCube, aes(x=LMA, y=NFPD(get(species[i]))), color="#842B3B") +
  geom_smooth(data = temp, aes(x=LMA, y=((shap_LMA *60) +0.4)), color="black") + 
  scale_y_continuous(name=paste(species[i], "\nnormalized detection rate"), sec.axis = sec_axis(~ . /60 -0.4/60, name="SHAP value")) +
  theme_minimal() +
  xlab("LMA (g/m2)") +
  xlim(c(0,400)) +
  theme(axis.title.y.left=element_text(color="#842B3B"))
ggsave("figures/GCKing_LMA.PDF", width=6, height=3)



















