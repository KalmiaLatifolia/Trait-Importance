################################################################################
# Figure4.R
# Trait importance manuscript
#
# Stacked bar chart of SHAP importance by species and BioCube category.
#
# Replaces the "stacked bar chart - Figure 4" section of full_protocol.R
# (approximately lines 1053-1112). One change: species are restricted to those
# whose mean full-model test R2 (NSE) exceeds zero, matching the exclusion rule
# stated in Methods 2.6 and applied in Figure3.R. This drops 8 of 94 species,
# leaving 86.
#
# Dropping those 8 has almost no effect on the average row: climate moves from
# 31.9 to 32.3 percent, foliar traits from 16.6 to 16.7, structure from 8.4 to
# 8.7, phenology 14.0 to 13.7, terrain 13.2 to 13.1, disturbance 15.9 to 15.5.
# The purpose of the change is consistency with Figure 3, not a different result.
#
# Written 10 August 2026
################################################################################

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)

# ---------------------------------------------------------------- settings ----

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

MODEL_PARAMS <- "data/xgb_modelParameters_20260622.rds"
VAR_IMP      <- "data/xgb_variableImportance_20260622.rds"
TABLE_S1     <- "data/TableS1_Biocube_var_description.xlsx"
OUTDIR       <- "figures"
DATESTAMP    <- format(Sys.Date(), "%Y%m%d")

# Set FALSE to reproduce the original all-species version of the figure
FILTER_R2 <- TRUE

dir.create(OUTDIR, showWarnings = FALSE)

# ------------------------------------------------------------------- data ----

xgb_modelParameters    <- readRDS(MODEL_PARAMS)
xgb_variableImportance <- readRDS(VAR_IMP)          # large file, slow to load

tidy <- read_excel(TABLE_S1, range = cell_cols(1:4))

# ------------------------------------------------------- species retention ----

keep <- xgb_modelParameters %>%
  filter(varSet == "spatVars") %>%
  group_by(species) %>%
  summarise(meanNSE = mean(R2_test), .groups = "drop") %>%
  filter(meanNSE > 0) %>%
  pull(species)

if (!FILTER_R2) keep <- unique(xgb_modelParameters$species)

message(sprintf("Figure 4 will show %d of %d species",
                length(keep), length(unique(xgb_modelParameters$species))))

# ----------------------------------------------------------- format data ----

temp <- xgb_variableImportance %>%
  subset(varSet == "spatVars") %>%          # only full models
  subset(species %in% keep) %>%             # <-- retained species only
  group_by(Feature, species) %>%
  summarise(SHAP_importance = mean(SHAP_importance, na.rm = TRUE), .groups = "drop") %>%
  rename(Variable = Feature) %>%
  left_join(tidy, by = "Variable") %>%
  group_by(species) %>%
  mutate(total   = sum(SHAP_importance, na.rm = TRUE),
         prop_T  = sum(SHAP_importance[Category == "Traits"], na.rm = TRUE) / total,
         prop_bar = SHAP_importance / total) %>%
  ungroup() %>%
  mutate(
    species  = reorder(species, prop_T),
    Category = factor(Category,
                      levels = c("Climate", "Disturbance", "Terrain",
                                 "Phenology", "Structure", "Traits"))) %>%
  arrange(species, Category)

catcols <- c("#70A4AF", "#7C4584", "#A8BBA3", "#D97C55", "#FDC71B", "#842B3B")

# ----------------------------------------------------------- create panels ----

p1 <- temp %>%
  ggplot(aes(y = species, x = SHAP_importance, fill = Category)) +
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = ifelse(prop_bar > 0.1, Label, ""), x = SHAP_importance / 2),
            position = position_fill(vjust = 0.5), size = 2) +
  scale_fill_manual(values = catcols) +
  theme_minimal() +
  xlab("SHAP importance (%)") +
  ylab("Species")

p2 <- temp %>%
  group_by(Category, species) %>%
  summarise(cat_size = sum(prop_bar), .groups = "drop") %>%
  group_by(Category) %>%
  summarise(mean_cat_prop = mean(cat_size), y = "Average", .groups = "drop") %>%
  ggplot(aes(y = y, x = mean_cat_prop, fill = Category)) +
  geom_bar(stat = "identity", position = position_fill()) +
  geom_text(aes(label = paste(round(mean_cat_prop, digits = 2) * 100, "%"),
                x = mean_cat_prop),
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = catcols) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("SHAP importance (%)") +
  ylab("")

# ------------------------------------------------------------------ output ----

ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(9, 1), align = "hv")

ggsave(file.path(OUTDIR, paste0("Figure4_SHAPimportance_", DATESTAMP, ".pdf")),
       height = 12, width = 12)
ggsave(file.path(OUTDIR, paste0("Figure4_SHAPimportance_", DATESTAMP, ".png")),
       height = 12, width = 12, dpi = 300)

# print the average-row percentages for the Results text
temp %>%
  group_by(Category, species) %>%
  summarise(cat_size = sum(prop_bar), .groups = "drop") %>%
  group_by(Category) %>%
  summarise(mean_pct = round(100 * mean(cat_size), 1), .groups = "drop") %>%
  arrange(desc(mean_pct)) %>%
  as.data.frame() %>%
  print()

################################################################################
# end
################################################################################
