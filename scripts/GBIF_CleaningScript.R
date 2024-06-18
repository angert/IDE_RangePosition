#### PROJECT: International Drought Experiment x Range Position
#### PURPOSE: Cleaning script for extracting species occurrence from GBIF
#### AUTHOR: Ryan Ng (ryanng3001@gmail.com)
#### DATE LAST MODIFIED: 20240424

# setwd() as IDE_RangePosition

# make vector of packages needed
packages_needed <- c("rgbif", "dplyr", "CoordinateCleaner", "sp", "raster", "ggplot2", "maps")

# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}


# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

sn2 <- c("LYSIMACHIA EUROPAEA")
sn <- c("Lysimachia europaea")


for (i in c(1)){
    
  # obtain species record from gbif 
  species_name2 <- sn2[i] # TODO 
  species_name <- sn[i] # TODO
  species_search <- occ_data(scientificName = species_name, 
                             country = "US;CA;MX",
                             hasCoordinate = TRUE, 
                             year = "1970,2000",
                             hasGeospatialIssue = FALSE,
                             coordinateUncertaintyInMeters = "0,1000",
                             limit = 99999) 
  
  # remove rows with occurrence issues 
  # reference: https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html
  species_search <- species_search %>%
    occ_issues(-bri, -ccm, -cdiv, -cdout, -cdpi, -cdrepf, -cdreps,
               -cdumi, -conti, -cucdmis, -cum, -gdativ, -geodi, 
               -geodu, -iccos, -iddativ, -iddatunl, -indci, -mdativ, 
               -mdatunl, -osifbor, -osu, -preneglat, -preneglon, 
               -preswcd, -rdativ, -rdatm, -rdatunl, -zerocd)
  
  # remove potential occurrences from Hawaii 
  species_search <- species_search$data %>% 
      filter(stateProvince != "Hawaii")
  
  # select relevant columns 
  species_data <- species_search[ , c("species", "decimalLongitude", 
                                                 "decimalLatitude", "issues",
                                                 "countryCode",
                                                 "occurrenceStatus", 
                                                 "coordinateUncertaintyInMeters", 
                                                 "institutionCode", "gbifID", 
                                                 "references", "basisOfRecord", 
                                                 "year", "month", "day", 
                                                 "eventDate", "geodeticDatum", 
                                                 "catalogNumber")]
  
  
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
  filename <- gsub(" ", "_", species_name2)
  output_file <- file.path(destination_folder, paste0(filename, ".csv"))
  write.csv(species_data_clean, file = output_file, row.names = FALSE)
  
  # update summary csv
  summary <- read.csv("data/cleaned_occurrence/summary.csv")
  summary$Final_Occ[summary$Taxon == species_name2] <- n
  summary$Completed[summary$Taxon == species_name2] <- "Y"
  write.csv(summary, file = file.path(destination_folder, "summary.csv"), row.names = FALSE)
  
  # create shapefile 
  species_data_lonlat <- species_data_clean[ , c("decimalLongitude", "decimalLatitude")]
  ch <- chull(species_data_lonlat)  
  coords <- species_data_lonlat[c(ch, ch[1]), ] 
  sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)), proj4string=CRS("+proj=longlat +datum=WGS84"))
  sp_poly_df <- SpatialPolygonsDataFrame(sp_poly, data=data.frame(ID=1))
  destination_folder2 <- file.path("data", "shapefiles")
  if (!file.exists(destination_folder2)) {
    dir.create(destination_folder2, recursive = TRUE)
  }
  output_file2 <- file.path(destination_folder2, paste0(filename, ".shp"))
  shapefile(x = sp_poly_df, file = output_file2, overwrite=TRUE)
  
  
  ###### sanity checks ##########################################
  ###############################################################
  nrow(species_data_clean)
  # Load map data for specific countries
  wm_filtered <- map_data("world") %>%
    filter(region %in% c("USA", "Canada", "Mexico"))
  # Plot the cropped map with points
  ggplot() + 
    geom_polygon(data = wm_filtered, 
                 aes(x = long, y = lat, group = group),
                 colour = "gray70", fill = "gray70") +
    geom_point(data = species_data_clean, 
               aes(x = decimalLongitude, y = decimalLatitude),
               colour = "black", bg = "olivedrab2", pch = 21, size = 0.75) +
    labs(x = "Longitude", 
         y = "Latitude") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
    coord_fixed(xlim = c(-125, -65), ylim = c(20, 70))
  ###############################################################
  ###############################################################
}

