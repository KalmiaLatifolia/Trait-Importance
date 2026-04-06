# Trait Importance
This repository contains code to accompany the future manuscript tentatively titled **The predictive power of new remote sensing data for biological distributions**

In its current state this is a working repository. Code is in flux and subject to change.

The goal of this project is to quantify the extent to which new remote sensing data 
like hyperspectrally derived foliar traits and GEDI lidar improve our ability to 
model species distributions and relative abundances.
A dataset of 140 remote sensing variables was used to predict the relative abundances of 94 bird species. 
We compared the accuracy of species models with and without different predictor variables to see how those predictors improved model accuracy.


![Protocol figure](figures/Figure1_Map.png)


# Data Availability

All of the derived datasets used in this project are available within this repository inside the "Data" folder.
Most files are available in both csv and rds formats.

| File Name | Description |
|-----------|-------------|
| siteDetections_foliarTraits_BioCube | The core dataset. One row per study site. Columns for each species and geospatial variable. |
| TableS1_Biocube_var_description | One row for each geospatial variable. Columns for variable descriptions, DOI of original publication, and link to dataset. |
| tidyNames | Reference sheet for full variable names used in data files versus pretty names used in figure axes. |
| xgb_modelParameters | xgboost outputs for each species. 100 iterations per species per model type. Model R2 and RMSE |


