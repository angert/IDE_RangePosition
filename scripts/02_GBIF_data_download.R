#### PROJECT: International Drought Experiment x Range Position
#### PURPOSE: download species occurrence data from GBIF
#### AUTHOR: Ryan Ng (ryanng3001@gmail.com)/ Allison Louthan 
#### DATE LAST MODIFIED: 250223
rm(list = ls())
# setwd() as IDE_RangePosition
# load packages & data-----
#
# make vector of packages needed
packages_needed <- c("rgbif", "dplyr", "CoordinateCleaner", "sp", "tidyverse",
                     "raster", "ggplot2", "maps","stringr")

# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}


# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

species_site_consecyears <- read.csv("data/01_species_site_consecyears.csv")
sn_cap <- species_site_consecyears$Taxon

# downloading raw data for all species on GBIF & getting citations for those raw data---
downloadkeys <- as.character(rep(NA, length(sn_cap)))
citations <- rep(NA, length(sn_cap))
  for (i in 1:length(sn_cap)) {
  species_name <- stringr::str_to_sentence(sn_cap[i] )
  
   # obtain species record from gbif 
  taxonKey_i <- rgbif::name_backbone(species_name) 

  species_search <- rgbif::occ_download_wait( # wait for an occurrence download to be done, because they must be done in sequence
    rgbif::occ_download(# this command improves on Jenn's student's code because 
                            # it gets unlimited records AND generates a citation
                            # and, all of these downloads will show up on AML's GBIF account online
    rgbif::pred("hasGeospatialIssue", FALSE),user= "amlouthan", pwd= "mYwhu3-wukfaz-qahpaf", email= "amlouthan@ksu.edu", 
    rgbif::pred("hasCoordinate", TRUE), 
    rgbif::pred("occurrenceStatus","PRESENT"),
    rgbif::pred_not(rgbif::pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))), 
    rgbif::pred("taxonKey", taxonKey_i$usageKey)  ,
    rgbif::pred_lte("year", 2000), # year is less than or equal to 2000
    rgbif::pred_gte("year", 1970),  # year is greater than or equal to 1970
    rgbif::pred_or(rgbif::pred_not(
        rgbif::pred_in("establishmentMeans",c("MANAGED","INTRODUCED"))),
        rgbif::pred_isnull("establishmentMeans")),  #establishmentMeans column does not contain managed or introduced species, but can be left blank.
     #NB that the above command tries (as best it can) to remove all introduced/ invasive/ non-native records
    # but it almost certainly will allow some introduced/ invasive/ non-native records!
    rgbif::pred_or(rgbif::pred_lt("coordinateUncertaintyInMeters",1000),rgbif::pred_isnull(
      "coordinateUncertaintyInMeters")), #coordinateUncertaintyInMeters is less 1000 meter or is left blank.
    format = "SIMPLE_CSV"
    )
  )
  downloadkeys[i] <- species_search$key
  citations[i] <-  rgbif::gbif_citation(species_search$key)$download
  write.csv(downloadkeys, file= "data/02_GBIF_download_keys.csv") # putting these write commands inside the loop because this loop fails all the time
  write.csv(citations, file= "data/02_GBIF_citations_for_pub.csv")
  }
write.csv(sn_cap, file= "data/02_GBIF_species_names.csv")
