# project: IDE (international drought experiment)
# code summary: get data prepared to calc position within range or niche for IDE sites
# written by: A Louthan
# modified: 241121

# libraries and functions -----
library(tidyr)
library(dplyr)

# load global info data-----
r <- raster::stack(
  geodata::worldclim_global("bio", 2.5, path=tempdir()))
r.matrix <- raster::getValues(r)
r.matrix <- r.matrix[complete.cases(r.matrix), ]
res <-prcomp(r.matrix, scale = TRUE, retx=TRUE) # run a PCA on the climate data; note that the PCA is computed across the globe, not just on one species' range

# which PCA axes correspond to mean-associated bioclim variables?
means.indices <- c(1, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19) 
r.matrix.means <- r.matrix[,means.indices]
res.means <-prcomp(r.matrix.means, scale = TRUE, retx=TRUE) # calc a global PCA with only mean-associated bioclim variables
rm(r.matrix, r.matrix.means,r) # remove these giant files

# which species are at which site?
species_site_counts <- read.csv("data/species_site_consecyears.csv")
species_site_counts_Control <- species_site_counts[,c(1, which(grepl("_Control", x= colnames(species_site_counts))))]
species_site_counts_Drought <- species_site_counts[,c(1, which(grepl("_Drought", x= colnames(species_site_counts))))]
colnames(species_site_counts_Control) <- gsub("\\_.*", "", colnames(species_site_counts_Control)) # remove _Drought or _Control suffix from colnames
colnames(species_site_counts_Drought) <- gsub("\\_.*", "", colnames(species_site_counts_Drought)) 
long_species_site_counts_Control <- species_site_counts_Control %>% 
  pivot_longer(
    cols= 'allmendb.ch':'ethadn.au',
    names_to = "site", 
    values_to = "count"
  )
long_species_site_counts_Drought <- species_site_counts_Drought %>% 
  pivot_longer(
    cols= 'allmendb.ch':'ethadn.au',
    names_to = "site", 
    values_to = "count"
  )
species_sites <- bind_rows(long_species_site_counts_Control, long_species_site_counts_Drought) %>%
  mutate_if(is.numeric, tidyr::replace_na, 0) %>% #in case of having NAs
  group_by(Taxon, site) %>%
  summarise_all(., sum, na.rm = TRUE)
species_sites <- species_sites[which(species_sites$count >0), c("Taxon", "site")]
colnames(species_sites) <- c("species","site")

# load IDE site locations into species_sites file
site_info <- read.csv("data/Site_Elev-Disturb241104.csv")[,c("site_code", "latitud", "longitud")]
colnames(site_info) <- c("site", "Lat", "Long")
species_sites <- left_join(species_sites, site_info, by= "site")
# IDE metadata says: Latitude	== latitude from -90 (S) to +90 (N) in decimal degrees	
# Longitude	== longitude from -180 (W) to +180 (E) in decimal degrees

save.image(file= "data/01_input_data.RData")
