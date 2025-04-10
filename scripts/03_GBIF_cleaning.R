#### PROJECT: International Drought Experiment x Range Position
#### PURPOSE: Cleaning script for extracting species occurrence from GBIF
#### AUTHOR: Ryan Ng (ryanng3001@gmail.com)/ Allison Louthan
#### DATE LAST MODIFIED: 250223
# start here: need to make sure this runs
# setwd() as IDE_RangePosition
rm(list = ls())

# load packages & data-----


# make vector of packages needed
packages_needed <- c("rgbif", "dplyr", "CoordinateCleaner", "sp", "tidyverse", "tidyr", 
                     "raster", "ggplot2", "maps","stringr", "foreach", "doParallel","rWCVP", "rWCVPdata")

# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}



# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

destination_folder <- file.path("data", "cleaned_occurrence")
if (!file.exists(destination_folder)) {
  dir.create(destination_folder, recursive = TRUE)
}

# load data---- 
downloadkeys <- read.csv("data/02_GBIF_download_keys.csv")[,2] 
sn_cap <- read.csv("data/02_GBIF_species_names.csv")[,2]
checklist <- rWCVP::wcvp_checklist(synonyms = TRUE)

# cleaning occurence data---- 
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
good_species <- # good species is a vector that tells you which species are ok to use (their GBIF name maps to one WCVP name)
  foreach(i= 1:length(sn_cap), .packages="tidyverse",.combine=cbind) %dopar% {
    good_species_i <- TRUE
  
  # download GBIF data from the GBIF website
  d <- rgbif::occ_download_get(key = downloadkeys[i], path= "data/raw_occurence", overwrite=FALSE) %>%
    rgbif::occ_download_import()


# remove rows with occurrence issues 
# reference: https://data-blog.gbif.org/post/gbif-filtering-guide/ (Jenn's student had used GBIF but GBIF recommends using coordinate clearer)

  d <- d %>%
  filter(coordinatePrecision < 0.01 | is.na(coordinatePrecision)) %>% 
  filter(!coordinateUncertaintyInMeters %in% c(301,3036,999,9999)) %>% # remove any records taht were assigned the default value for coordinate uncertainty-- often erroneous
  filter(!decimalLatitude == 0 | !decimalLongitude == 0) %>%
  CoordinateCleaner::cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
    CoordinateCleaner::cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
    CoordinateCleaner::cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
    CoordinateCleaner::cc_sea() %>% # remove from ocean 
  distinct(decimalLongitude,decimalLatitude,speciesKey,datasetKey, .keep_all = TRUE) # removes potentially duplicated recrods

  if (dim(d)[1]== 0 | i==134 | i == 182 |i == 186
      # 134 contains >1 species name because GBIF download could not find the species.
      # For 182, there was no homotypic synonym in the WCVP checklist.
      # 186 contains >1 species name (all in the same genus), but those species names map to different non-homotypic synonyms in the WCVP checklist
      
      # In the methods, should say "we used species who we were able to map to a singular species name on GBIF and in the WCVP checklist"
      ) {good_species_i <- FALSE} else {
# select relevant columns 
d <- d[ , c("species","scientificName", "decimalLongitude", 
                                         "decimalLatitude")]

# now, remove all non-native occurences, using this protocol: https://matildabrown.github.io/rWCVP/articles/coordinate-cleaning.html
if (length(unique(d$species)) != 1) {stop("GBIF download contains multiple species names")} 

rWCVP_name_df <- rWCVP::wcvp_match_names(as.data.frame(unique(d$species)), 
                                      name_col ='unique(d$species)')
rWCVP_name <- rWCVP_name_df[which(rWCVP_name_df$wcvp_status == "Accepted"), "wcvp_name"] 
if (length(rWCVP_name)>1) {good_species_i <- FALSE} # FALSE species cannot be used b/c their GBIF name maps onto many WCVP names, so we cannot deliniate its native range
if (length(rWCVP_name) == 0){ # if there are no "Accepted" WCVP matches in terms of species names, then match the Synonym
 if (length(which(rWCVP_name_df$wcvp_status == "Synonym" & rWCVP_name_df$wcvp_homotypic)) == 1){
   homotypic_synonyms <- which(rWCVP_name_df$wcvp_status == "Synonym" & rWCVP_name_df$wcvp_homotypic)
   rWCVP_name_acceptedID <- rWCVP_name_df[homotypic_synonyms, "wcvp_accepted_id"] 
   rWCVP_name <- unique(checklist[which(checklist$plant_name_id == rWCVP_name_acceptedID),"accepted_name"])
   }
   else {stop("too few or too many appropriate synonyms; check rWCVP_name_df")}
}

native_range <- rWCVP::wcvp_distribution(rWCVP_name, taxon_rank="species",
                                  introduced=FALSE, extinct=FALSE, 
                                  location_doubtful=FALSE) %>% 
  terra::vect() %>% terra::buffer(width= 2*1000) %>%  # 1 km buffer (as recommended by rWCVP), using terra instead of st_buffer because terra works with unprojected (CRS 4326) coordinates
  terra::union()

occs <- # now back to the GBIF data-- transofrm into sf
  d %>% 
  dplyr::select(scientificName, decimalLatitude, decimalLongitude) %>%
  sf::st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), crs = sf::st_crs(4326)) %>% 
  terra::vect()

occs <- occs[native_range] # subsets occ by the native range (so, keeps only those occ who intersect the native range)
occs <- terra::crds(occs) # extracts the coordinates
# save into data/cleaned_occurrence folder
output_file <- file.path(destination_folder, paste0(gsub(" ", "_", unique(d$species)), ".csv"))
write.csv(occs, file = output_file, row.names = FALSE)
}

good_species_i
}

stopCluster(cl)


write.csv(good_species, file= "data/03_good_species.csv")
