#### PROJECT: International Drought Experiment x Range Position
#### PURPOSE: Cleaning script for extracting species occurrence from GBIF as .csv
#### AUTHOR: Ryan Ng
#### DATE LAST MODIFIED: 20240411

# setwd() as IDE_RangePosition


# make vector of packages needed
packages_needed <- c("rgbif", "dplyr", "CoordinateCleaner")


# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}


# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}


# obtain species record from gbif 
species_name <- "BOUTELOUA GRACILIS" # TODO
species_search <- occ_data(scientificName = species_name, 
                           country = "US",
                           hasCoordinate = TRUE, 
                           year = "1970,2000",
                           hasGeospatialIssue = FALSE,
                           coordinateUncertaintyInMeters = "0,1000",
                           limit = 99999) 


# remove rows with occurrence issues 
# reference: https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html
species_issues_clean <- species_search %>%
  occ_issues(-bri, -ccm, -cdiv, -cdout, -cdpi, -cdrepf, -cdreps,
             -cdumi, -conti, -cucdmis, -cum, -gdativ, -geodi, 
             -geodu, -iccos, -iddativ, -iddatunl, -indci, -mdativ, 
             -mdatunl, -osifbor, -osu, -preneglat, -preneglon, 
             -preswcd, -rdativ, -rdatm, -rdatunl, -typstativ, 
             -zerocd)


# select relevant columns 
species_data <- species_issues_clean$data[ , c("species", "decimalLongitude", 
                                               "decimalLatitude", "issues",
                                               "countryCode", "individualCount", 
                                               "occurrenceStatus", 
                                               "coordinateUncertaintyInMeters", 
                                               "institutionCode", "gbifID", 
                                               "references", "basisOfRecord", 
                                               "year", "month", "day", 
                                               "eventDate", "geodeticDatum", 
                                               "datasetName", "catalogNumber")]


# remove invalid rows 
species_data <- species_data%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))%>%
  filter(!is.na(countryCode))%>%
  filter(!is.na(occurrenceStatus))%>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(!is.na(institutionCode))%>%
  filter(!is.na(geodeticDatum))%>%
  filter(occurrenceStatus != "ABSENT")%>% 
  filter(countryCode != "XK") 


# identify and remove flagged records 
species_data_clean <- species_data%>%
  cc_val()%>% # invalid lat/lon coordinates
  cc_equ()%>% # identical lat/lon
  cc_cap()%>% # coordinates in vicinity of country capitals
  cc_cen()%>% # coordinates in vicinity of country or province centroids
  cc_gbif()%>% # coordinates assigned to GBIF headquarters
  cc_inst()%>% # coordinates in the vicinity of biodiversity institutions (slow)
  cc_zero()%>% # coordinates that are zero
  cc_sea()%>% # coordinates at sea
  cc_outl()%>% # outliers (slow)
  cc_dupl() # duplicate species


# save into data/cleaned_occurrence folder
n <- nrow(species_data_clean) 
destination_folder <- file.path("data", "cleaned_occurrence")
if (!file.exists(destination_folder)) {
  dir.create(destination_folder, recursive = TRUE)
}
filename <- gsub(" ", "_", species_name)
output_file <- file.path(destination_folder, paste0(filename, ".csv"))
write.csv(species_data_clean, file = output_file, row.names = FALSE)

# update summary csv
summary <- read.csv("data/cleaned_occurrence/summary.csv")
summary$Final_Occ[summary$Taxon == species_name] <- n
summary$Completed[summary$Taxon == species_name] <- "Y"
write.csv(summary, file = file.path(destination_folder, "summary.csv"), row.names = FALSE)