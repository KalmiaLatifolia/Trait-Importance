
# Supplemental figures
# species SHAP importance
# 9 April 2026
# Laura Berman

# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# load data --------------------------------------------------------------------

xgb_variableImportance <- readRDS("data/xgb_variableImportance_20260325.rds")

siteDetections_foliarTraits_BioCube <- readRDS("data/siteDetections_foliarTraits_BioCube_20260320.rds")
species <- colnames(siteDetections_foliarTraits_BioCube)[4:97]

# load tidy names
tidy <- read_excel("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TableS1_Biocube_var_description.xlsx",
                   range = cell_cols(1:4))

# pick one species
A <- subset(xgb_variableImportance, xgb_variableImportance$species=="Acorn Woodpecker")

# set order
feat_order <- A %>%
  group_by(Feature) %>%
  summarise(m = mean(SHAP_importance, na.rm = TRUE), .groups="drop") %>%
  arrange(m) %>%
  pull(Feature)

A$Feature <- factor(A$Feature, levels = feat_order)

# plot it 
ggplot(A, aes(x=SHAP_importance, y=Feature)) + 
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.y = element_text(colour = tidy$labelColor[match(feat_order, tidy$Variable)])) +
  scale_y_discrete(labels = tidy$Label[match(feat_order, tidy$Variable)]) +
  ggtitle("Acorn Woodpecker")
ggsave("AcornWoodpecker_SHAP.pdf", height=15, width=6)


lapply(species, function(sp){
  
  A <- subset(xgb_variableImportance, species == sp)
  
  feat_order <- A %>% group_by(Feature) %>%
    summarise(m = mean(SHAP_importance, na.rm=TRUE), .groups="drop") %>%
    arrange(m) %>% pull(Feature)
  
  A$Feature <- factor(A$Feature, levels=feat_order)
  
  p <- ggplot(A, aes(x=SHAP_importance, y=Feature)) +
    geom_boxplot() + theme_minimal() +
    theme(axis.text.y = element_text(colour = tidy$labelColor[match(feat_order, tidy$Variable)])) +
    scale_y_discrete(labels = tidy$Label[match(feat_order, tidy$Variable)]) +
    ggtitle(sp)
  
  ggsave(paste0(gsub(" ", "", sp), "_SHAP.pdf"), plot=p, height=15, width=6)
  
})


