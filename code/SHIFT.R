
# SHIFT
# SHIFT data shared by Ting
# 11 Dec 2025


# set working directory --------------------------------------------------------

setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Trait importance/TraitImportance_GIT")

# read in SHIFT data -----------------------------------------------------------

SHIFTdata <- read.csv("data/shift_sample_traits_cleaned4Laura.csv")

# plot -------------------------------------------------------------------------

unique(temp$species_or_type)
temp <- subset(temp, temp$species_or_type=="Artemisia californica - California sagebrush" |
                 temp$species_or_type=="Salvia leucophylla - Purple sage" |
                 temp$species_or_type=="Salvia mellifera - Black sage" |
                 temp$species_or_type=="Heteromeles arbutifolia - Toyon" |
                 temp$species_or_type=="Salix lasiolepis - Arroyo willow" |
                 temp$species_or_type=="Populus trichocarpa - Black cottonwood" |
                 temp$species_or_type=="Quercus berberidifolia - Scrub oak" |
                 temp$species_or_type=="Quercus douglasii - Blue oak" |
                 temp$species_or_type=="Quercus crysolepis - Canyon live oak" |
                 temp$species_or_type=="Juniperus californica - California juniper" |
                 temp$species_or_type=="Pseudotsuga macrocarpa - Bigcone douglas fir" |
                 temp$species_or_type=="Pinus ponderosa - Ponderosa pine" |
                 temp$species_or_type=="Pinus jeffreyi - Jeffrey pine" |
                 temp$species_or_type=="Pinus sabiniana - Gray pine" |
                 temp$species_or_type=="Pinus coulteri - Coulter pine")

#sulfur
temp <- subset(SHIFTdata, !is.na(SHIFTdata$sulfur))
ggplot(temp, aes(x=sulfur, y=reorder(species_or_type, sulfur, FUN = mean))) +
  geom_boxplot(color="#842B3B", fill="#B4777C") +
  theme_minimal() +
  ylab("") +
  xlab("Sulfur") +
  xlim(0.5,3) +
  geom_text(data = temp %>% group_by(species_or_type) %>% summarise(n = n()), 
            aes(x = 3, y = species_or_type, label = paste0("n=", n)), 
            vjust = 0)
ggsave("SHIFT_Sulfur.PDF", height=4, width=8)

# potassium
temp <- subset(SHIFTdata, !is.na(SHIFTdata$potassium))
ggplot(temp, aes(x=potassium, y=reorder(species_or_type, potassium, FUN = mean))) +
  geom_boxplot(color="#842B3B", fill="#B4777C") +
  theme_minimal() +
  ylab("") +
  xlab("Potassium") +
  xlim(0,15) +
  geom_text(data = temp %>% group_by(species_or_type) %>% summarise(n = n()), 
            aes(x = 50, y = species_or_type, label = paste0("n=", n)), 
            vjust = 0)
ggsave("SHIFT_potassium.PDF", height=4, width=8)

# phosphorus
temp <- subset(SHIFTdata, !is.na(SHIFTdata$phosphorus))
ggplot(temp, aes(x=phosphorus, y=reorder(species_or_type, phosphorus, FUN = mean))) +
  geom_boxplot(color="#842B3B", fill="#B4777C") +
  theme_minimal() +
  ylab("") +
  xlab("Phosphorus") +
  geom_text(data = temp %>% group_by(species_or_type) %>% summarise(n = n()), 
            aes(x = 10, y = species_or_type, label = paste0("n=", n)), 
            vjust = 0)

# phenolics
temp <- subset(SHIFTdata, !is.na(SHIFTdata$phenolics))
ggplot(temp, aes(x=phenolics, y=reorder(species_or_type, phenolics, FUN = mean))) +
  geom_boxplot(color="#842B3B", fill="#B4777C") +
  theme_minimal() +
  ylab("") +
  xlab("Phenolics") +
  xlim(30, 120) +
  geom_text(data = temp %>% group_by(species_or_type) %>% summarise(n = n()), 
            aes(x = 200, y = species_or_type, label = paste0("n=", n)), 
            vjust = 0)
ggsave("SHIFT_phenolics.PDF", height=4, width=8)
