# project: IDE (international drought experiment)
# code summary: calculate posiition in range or niche for IDE sites
# written by: A Louthan
# modified: 250410

# libraries and functions -----
rm(list = ls())
library(raster)
library(sp)
c.fun<-function(df, center, scale) {
  return((df-center)/scale )
} 

# load in data----
# load terrestial map area from: https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-land/
terrestrial_sf <- rnaturalearth::ne_download(scale= 10, type= "land", "physical")

extract_coords1 <- function(x,ii){x@coords[,ii]}
S.pole.pt <- SpatialPoints(coords= cbind(x= 0, y= -90), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) # making a Spatial points object of the N pole
N.pole.pt <- SpatialPoints(coords= cbind(x= 0, y= 90), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) # making a Spatial points object of the N pole


r <- raster::stack( # I am re-loading this b/c if I load it as input data I get a weird usage error
  geodata::worldclim_global("bio", 2.5, path=tempdir()))

# load pca & species-site combos
load("data/04_input_data.RData")
species_sites <- species_sites[which(species_sites$species %in% gsub("_", " ",gsub(".csv", "",list.files("data/cleaned_occurrence")))),]
species_sites <- species_sites[which(!is.na(species_sites$Lat) & !is.na(species_sites$Long)), ]
GBIF_species_names_good_and_bad <- read.csv("data/02_GBIF_species_names.csv")[,2]
good_species <- as.vector(read.csv("data/03_good_species.csv")[1,-1], mode="logical") # first column in this data frame is not true/ false
if (length(GBIF_species_names_good_and_bad) != length(good_species)) {stop("good_species and GBIF name list do not have the same number of entries-- something is wrong")}
GBIF_species_names <- GBIF_species_names_good_and_bad[good_species]

# initialize the species-specific loop---- 
alpharange <- 15 # for the alpha shapes, alpha = 15 is the default for geographic range
alphaniche <- 10 # for the alpha shapes, alpha = 10 is the default for niche
max.alpha <- 400 # the max value-- at which you use the MCP-- is 400
alpha.step <- 5 # and this is the increment of increase

m1 <- as.data.frame(matrix(ncol=24, nrow= nrow(species_sites)))
names(m1) <- c( "rangearea", 
                "rangepcaarea", 
                "cat.descriptor", 
                "cat.descriptor.eqsp", 
                "rangeposition.quantile", 
                "rangeposition.quantile.eqsp", 
                "cat.pca.descriptor",
                "rangepcaposition.quantile", 
                
                "rangearea.95", 
                "rangepcaarea.95",  
                "cat.descriptor.95",  
                "cat.descriptor.eqsp.95", 
                "rangeposition.quantile.95",  
                "rangeposition.quantile.eqsp.95", 
                "cat.pca.descriptor.95", 
                "rangepcaposition.quantile.95", 
                
                "rangearea.alphahull", 
                "rangepcaarea.alphahull",  
                "cat.descriptor.alphahull",  
                "cat.descriptor.eqsp.alphahull", 
                "rangeposition.quantile.alphahull",  
                "rangeposition.quantile.eqsp.alphahull", 
                "cat.pca.descriptor.alphahull", 
                "rangepcaposition.quantile.alphahull")
species_sites_positions <- cbind(species_sites, m1)
sf::sf_use_s2(FALSE)
options(warn= 2)
cat.quantile <- 0.2 # what percent of the species' range is considered the "edge" (0.2 = outer 20% considered the edge, etc.)



# the loop-----

for (i in 1:length(GBIF_species_names)){ # cycling through all the unique species in the list
  species_name_i <- GBIF_species_names[i]
  MCP.range <- MCP.niche <-  "N"  # this term indicates whether the alpha shape is actually the MCP; default is no
  species_csv_i <- read.csv(paste("data/cleaned_occurrence/", sub(" ", "_", species_name_i), ".csv", sep= "")) # this gives the csv locations for the species of interest
  names(species_csv_i) <- c("decimalLongitude", "decimalLatitude")
  species_csv_i <- species_csv_i[complete.cases(cbind(as.numeric(species_csv_i$decimalLongitude), as.numeric(species_csv_i$decimalLatitude))), ] # removes any rows that do not have BOTH lat and long
  
  
  species.indices <- which(species_sites_positions$species == species_name_i)
  demo.speciesi <- species_sites_positions[species.indices,] # these two lines subsets large data file to include only species i
  results <- matrix(NA, nrow= dim(demo.speciesi)[1], ncol= ncol(m1)) # sets up a location to store the results
  
  
  # construct range maps for convex hull, alphahull, 95% convex hull----
  
  # first, construct a SpatialPolygons object using herbarium records and demo sites
  unique_demo_lat_longs <- !duplicated(cbind(demo.speciesi$Long, demo.speciesi$Lat)) # remove duplicate lat/longs in demo dataset
  x <- as.numeric(c(species_csv_i$decimalLongitude, demo.speciesi$Long[unique_demo_lat_longs])) # x= longitude of demo points
  y <- as.numeric(c(species_csv_i$decimalLatitude, demo.speciesi$Lat[unique_demo_lat_longs])) # y is latitude of demo points
  mymat <- cbind(x,y)
  speciespoints_i <- SpatialPoints(coords = mymat, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  polygon <- chull(x, y)
  polygon <- c(polygon, polygon[1])
  speciesrange_i <- SpatialPolygons(list(sp::Polygons(list( sp::Polygon(mymat[polygon,] 
                                                                        , hole= FALSE)), ID = "A")), 
                                    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) # the location data, when imported, is currently in Lat/ Long
  
  # next, calc the convex hull
  terrestrial_sf_projected <- sf::st_transform(terrestrial_sf, crs= crs(speciesrange_i)@projargs)
  overlap <- sf::st_intersection(sf::st_geometry(terrestrial_sf_projected), sf::st_geometry(sf::st_as_sf(speciesrange_i))) # confine convex hull to only terrestrial area
  rangearea <- as.numeric(sum( sf::st_area(overlap)) )# in m^2
  speciesrange_i_absyrange <- abs(range(as.vector(sapply(speciesrange_i@polygons[[1]]@Polygons, extract_coords1, ii=2, simplify= "array"))))
  # gets span of species latitude, assuming reflection over equator for Southern hemisphere 
  speciesrange_i_xrange <- (range(as.vector(sapply(speciesrange_i@polygons[[1]]@Polygons, extract_coords1, ii=1, simplify= "array"))))# gets span of species longitude, no correction for meridian
  speciesrange_i_yrange <- (range(as.vector(sapply(speciesrange_i@polygons[[1]]@Polygons, extract_coords1, ii=2, simplify= "array"))) )# gets span of species latitude, with negatives and positives
  
  #construct an alpha hull (which allows for irregular shapes & holes)
  speciespoints_i_foralphahull <- remove.duplicates(speciespoints_i, zero = 1e-02, remove.second = TRUE, memcmp = TRUE) # alpha hull does not like duplicates
  options(warn=0)
  my.alpha.shape <- try(alphahull::ahull(x = speciespoints_i_foralphahull@coords[,1], y = speciespoints_i_foralphahull@coords[,2], alpha=alpharange), 
                        silent=TRUE)
  rangearea.alphahull <-  try(alphahull::areaahull(my.alpha.shape, alpha=alpharange),  silent=TRUE) # then calculate the area of the alpha shape
  alpha.used <- alpharange
  options(warn=2)
  # if construction of the alpha shape or calc of alpha shape area fails when running the (finnicky) ahull command
  if (class(my.alpha.shape)== "try-error" | class(rangearea.alphahull) == "try-error") { 
    #then first try the max.alpha, because weird-shaped ranges tend to never have an alpha value that works
    options(warn=0)
    my.alpha.shape <- try(alphahull::ahull(x = speciespoints_i_foralphahull@coords[,1], y = speciespoints_i_foralphahull@coords[,2], alpha=max.alpha), 
                          silent=TRUE)
    rangearea.alphahull <-  try(alphahull::areaahull(my.alpha.shape, alpha=max.alpha),  silent=TRUE) # then calculate the area of the alpha shape
    options(warn=2)
    if (class(my.alpha.shape)== "try-error" | class(rangearea.alphahull) == "try-error") {
      MCP.range <- "Y" # if the max.alpha value does not result in a valid alpha shape, then the alpha shape should be the MCP
      
      } else if (
        class(my.alpha.shape)!= "try-error" & class(rangearea.alphahull) != "try-error") {
# if the max.alpha value does result in a valid alpha shape, find the smallest alpha value that results in a valid alpha shape
  for (kk in seq(from = alpharange+ alpha.step, to = max.alpha, by= alpha.step)){
       options(warn=0)
       my.alpha.shape <- try(alphahull::ahull(x = speciespoints_i_foralphahull@coords[,1], y = speciespoints_i_foralphahull@coords[,2], alpha=kk), 
                             silent=TRUE)
       rangearea.alphahull <-  try(alphahull::areaahull(my.alpha.shape, alpha=kk),  silent=TRUE) # then calculate the area of the alpha shape
       alpha.used <- kk
       options(warn=2)
       if (class(my.alpha.shape) != "try-error" & class(rangearea.alphahull) != "try-error") {break} # if the tested alpha value results in a valid alpha shape, stop the loio    
       } # end kk loop
    } 
    }  
  
  if (MCP.range !=  "Y"){
# if neither the ahull nor the areaahull command fails, then get the "track" of the alpha shape
      alpha.shape.track <- alphahull::ahull_track(x = speciespoints_i_foralphahull@coords[,1], y = speciespoints_i_foralphahull@coords[,2], nps= 1000, alpha=alpha.used)
      all.alpha.coords <- NULL
      for (no_segments in 1:length(alpha.shape.track)){ 
        all.alpha.coords <-  rbind(all.alpha.coords,alpha.shape.track[[no_segments]]$data)}
      speciesrange_i_absyrange.alphahull <- abs(range(all.alpha.coords$y)) # getting the "track" of the alpha shape allows you to calc the x and y range of the shape
      speciesrange_i_xrange.alphahull <- range(all.alpha.coords$x)
      speciesrange_i_yrange.alphahull <-  range(all.alpha.coords$y)
      # NB we do not confine this range area to terrestrial b/c alphahulls allow for holes in the species range
}    

  
  #get the 95% convex hull
 if (length(speciespoints_i) >= 5){
   range.95 <- adehabitatHR::mcp(speciespoints_i, percent= 95, unin= "m", unout= "m2") 
  range.95 <- sf::st_intersection(sf::st_geometry(terrestrial_sf_projected), sf::st_geometry(sf::st_as_sf(range.95))) # calculate the area of the hull that intersects with terrestrial area
  rangearea.95 <- as.numeric(sum( sf::st_area(range.95))) # units on this range area are m^2
  speciesrange_i_absyrange.95 <- abs(c(as.numeric(attributes(range.95)$bbox$ymin), as.numeric(attributes(range.95)$bbox$ymax)))
  speciesrange_i_xrange.95 <-  c(as.numeric(attributes(range.95)$bbox$xmin), as.numeric(attributes(range.95)$bbox$xmax))
  speciesrange_i_yrange.95 <-  c(as.numeric(attributes(range.95)$bbox$ymin), as.numeric(attributes(range.95)$bbox$ymax)) } 
  
  # for species whose ranges overlap the equator, make a north and a south hemisphere range map----- 
  
  # N and S hemisphere range map for species whose ranges overlap the equator for convex hull and alpha hull ; note they share range extents
  if (min(speciesrange_i_yrange) <0 & max(speciesrange_i_yrange)>0){ # if the species' range overlaps the equator
    

    # divide up the species' range into a North hemisphere part & a South hemisphere part
    # note that we have to generate hemisphere-specific range polygons, rather than just using the extremes of points in the 
    # North and South hemisphere,
    # because the longitudinal range along the equator could be greater than the longitudinal range in one hemisphere
    
    # North first
    species_csv_i.North <- species_csv_i[which(species_csv_i$decimalLatitude >=0),]
    demo.speciesi.North <- demo.speciesi[unique_demo_lat_longs,]
    demo.speciesi.North <- demo.speciesi.North[which(demo.speciesi.North$Lat >=0),]
    
    x.North <- as.numeric(c(species_csv_i.North$decimalLongitude, demo.speciesi.North$Long)) # x= Long 
    x.North <- c(x.North, mean(x.North)) # adds an additional point to the distribution that is ON the equator
    y.North <- as.numeric(c(species_csv_i.North$decimalLatitude, demo.speciesi.North$Lat)) # y is Lat
    y.North <- c(y.North, 0) # adds an additional point to the distribution that is ON the equator
    
    if (length(y.North)==2) {stop("there is only one point in the Northern hemisphere; check range mpas")}
    # make a spatial polygon using only North hemisphere data
    mymat.North <- cbind(x.North,y.North)
    polygon.North <- chull(x.North, y.North)
    polygon.North <- c(polygon.North, polygon.North[1])
    speciesrange_i.North <- SpatialPolygons(list(sp::Polygons(list( 
      sp::Polygon(mymat.North[polygon.North,] 
                  , hole= FALSE)), ID = "A")), 
      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
    
    # calc N hemisphere-specific range extent
    speciesrange_i_absyrange.North <- abs(range(as.vector(sapply(speciesrange_i.North@polygons[[1]]@Polygons, extract_coords1, ii=2, simplify= "array"))))
    speciesrange_i_xrange.North <- abs(range(as.vector(sapply(speciesrange_i.North@polygons[[1]]@Polygons, extract_coords1, ii=1, simplify= "array"))))
    speciesrange_i_yrange.North <- (range(as.vector(sapply(speciesrange_i.North@polygons[[1]]@Polygons, extract_coords1, ii=2, simplify= "array"))))
    
    # South next
    species_csv_i.South <- species_csv_i[which(species_csv_i$decimalLatitude <0),]
    demo.speciesi.South <- demo.speciesi[unique_demo_lat_longs,]
    demo.speciesi.South <- demo.speciesi.South[which(demo.speciesi.South$Lat <0),]
    x.South <- as.numeric(c(species_csv_i.South$decimalLongitude, demo.speciesi.South$Long)) 
    x.South <- c(x.South, mean(x.South)) 
    y.South <- as.numeric(c(species_csv_i.South$decimalLatitude, demo.speciesi.South$Lat))
    y.South <- c(y.South, 0) 
  
    if (length(y.South)==2) {stop("there is only one point in the Southern hemisphere; check range mpas")}
    
    # make a spatial polygon using only South hemisphere data 
    mymat.South <- cbind(x.South,y.South)
    polygon.South <- chull(x.South, y.South)
    polygon.South <- c(polygon.South, polygon.South[1])
    speciesrange_i.South <- SpatialPolygons(list(sp::Polygons(list( 
      sp::Polygon(mymat.South[polygon.South,] 
                  , hole= FALSE)), ID = "A")), 
      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
    speciesrange_i_absyrange.South <- abs(range(as.vector(sapply(speciesrange_i.South@polygons[[1]]@Polygons, extract_coords1, ii=2, simplify= "array"))))
    speciesrange_i_xrange.South <- (range(as.vector(sapply(speciesrange_i.South@polygons[[1]]@Polygons, extract_coords1, ii=1, simplify= "array"))))
    speciesrange_i_yrange.South <- (range(as.vector(sapply(speciesrange_i.South@polygons[[1]]@Polygons, extract_coords1, ii=2, simplify= "array"))) )
    
  }
  
  # N and S hemisphere range map for species whose ranges overlap the equator for 95% convex hull
  if (length(speciespoints_i) >= 5){
  if (min(speciesrange_i_yrange.95) <0 & max(speciesrange_i_yrange.95)>0) {# if the 95% range spans the equator
    
    polygon.coords.95 <- as.data.frame(raster::geom(as(range.95, "Spatial")))
    polygon.coords.95 <- polygon.coords.95[which(polygon.coords.95$hole==0), ]
    speciesrange_i_xrange.North.95 <- range(polygon.coords.95$x[which(polygon.coords.95$y >=0)])
    speciesrange_i_yrange.North.95 <- c(0, max(polygon.coords.95$y[which(polygon.coords.95$y >= 0)]))
    speciesrange_i_absyrange.North.95 <- abs(speciesrange_i_yrange.North.95)
    speciesrange_i_xrange.South.95 <- range(polygon.coords.95$x[which(polygon.coords.95$y <0)])
    speciesrange_i_yrange.South.95 <- c(0, min(polygon.coords.95$y[which(polygon.coords.95$y < 0)]))
    speciesrange_i_absyrange.South.95 <- abs(speciesrange_i_yrange.South.95)
    
    if (speciesrange_i_xrange.North.95[1]==speciesrange_i_xrange.North.95[2]){ stop("only one point in N hemisphere for 95%; check range for this species")}
    if (speciesrange_i_xrange.South.95[1]==speciesrange_i_xrange.South.95[2]){ stop("only one point in S hemisphere for 95%; check range for this species")}
  } }
  
  if (MCP.range== "N"){ 
    if (min(speciesrange_i_yrange.alphahull) <0 & max(speciesrange_i_yrange.alphahull)>0) {# if the range spans the equator
      speciesrange_i_xrange.North.alphahull <- range(all.alpha.coords$x[which(all.alpha.coords$y >=0)])
      speciesrange_i_xrange.South.alphahull <- range(all.alpha.coords$x[which(all.alpha.coords$y <0)])
      speciesrange_i_yrange.North.alphahull <- c(0, max(all.alpha.coords$y[which(all.alpha.coords$y >= 0)]))
      speciesrange_i_absyrange.North.alphahull <- abs(speciesrange_i_yrange.North.alphahull)
      speciesrange_i_yrange.South.alphahull <- c(0, min(all.alpha.coords$y[which(all.alpha.coords$y < 0)]))
      speciesrange_i_absyrange.South.alphahull <- abs(speciesrange_i_yrange.South.alphahull)}
  }   
  
  # make range maps in PCA space-----
  # get Bioclim data for the whole species range
  range.temp.precip <- raster::extract(r, speciespoints_i, na.omit=TRUE) # NA's here are in the ocean, so the na.omit takes those out
  # get PCA scores for demography points
  centeredNewData <-apply(range.temp.precip, MARGIN=1, FUN=c.fun, res$center, res$scale  ) # first step in getting the PCA score for the points within the species range
  species_pca_i <-t(t(res$rotation) %*% centeredNewData)[,1:2] # second step; this matrix represent the scores for PCA axes 1 & 2 for each grid cell across the range.temp.precip
  
  # then construct the range map in pca space-----
  
  # first, the convex hull around occurrences of this species in pca space
  x <- c(species_pca_i[,1]) # x= pca score 1
  y <- c(species_pca_i[,2]) # y= pca score 2
  mymat <- cbind(x,y)[which(complete.cases(cbind(x,y))),] # making a matrix of points that are actually complete cases; I am not sure why you need the complete.cases
  # get a SpatialPolygon around the occurrences
  polygon <- chull(mymat[,1], mymat[,2]) 
  polygon <- c(polygon, polygon[1]) 
  speciespcarange_i <- SpatialPolygons(list(sp::Polygons(list(
    sp::Polygon(rbind(mymat[polygon,]) , hole= FALSE)), ID = "A"))) 
  rangepcaarea <- sf::st_area(sf::st_as_sf(speciespcarange_i), which= "Euclidean") # area of species' pca range (in pca space)
  
  speciespcarange_i_yrange <-(range(as.vector(sapply(speciespcarange_i@polygons[[1]]@Polygons, extract_coords1, ii=2, simplify= "array")))) # gettng the y range of the pca polygon (y= pca 2)
  speciespcarange_i_xrange <- (range(as.vector(sapply(speciespcarange_i@polygons[[1]]@Polygons, extract_coords1, ii=1, simplify= "array")))) # getting the x range of the pca polygon (x= pca 1)
  
  # next, construct the alpha shape around occurrences of this species in pca space 
  options(warn=0)
  mymat.df <- as.data.frame(mymat[!duplicated(mymat),])
  niche.alphahull <- try(alphahull::ahull(x = mymat.df$x, y = mymat.df$y,alpha=alphaniche), silent=TRUE)
  rangepcaarea.alphahull <-  try(alphahull::areaahull(niche.alphahull), silent=TRUE)
  options(warn=2)
  # if construction of the alpha hull fails in the ahull command or the area command... 
  if (class(niche.alphahull)== "try-error" | class(rangepcaarea.alphahull) == "try-error") { 
    #then first try the max.alpha, because weird-shaped ranges tend to never have an alpha value that works
    options(warn=0)
    niche.alphahull <- try(alphahull::ahull(x = mymat.df$x, y = mymat.df$y, alpha=max.alpha), 
                           silent=TRUE)
    rangepcaarea.alphahull <-  try(alphahull::areaahull(niche.alphahull, alpha=max.alpha),  silent=TRUE) # then calculate the area of the alpha shape
    options(warn=2)
    if (class(niche.alphahull)== "try-error" | class(rangepcaarea.alphahull) == "try-error") {
      MCP.range <- "Y" # if the max.alpha value does not result in a valid alpha shape, then the alpha shape should be the MCP
    } else if (class(niche.alphahull)!= "try-error" & class(rangepcaarea.alphahull) != "try-error") {
      # if the max.alpha value does result in a valid alpha shape, find the smallest alpha value that results in a valid alpha shape
      for (kk in seq(from = alphaniche+ alpha.step, to = max.alpha, by= alpha.step)){
        options(warn=0)
        niche.alphahull <- try(alphahull::ahull(x = mymat.df$x, y = mymat.df$y, alpha=kk), 
                               silent=TRUE)
        rangepcaarea.alphahull <-  try(alphahull::areaahull(niche.alphahull, alpha=kk),  silent=TRUE) # then calculate the area of the alpha shape
        options(warn=2)
        alpha.used <- kk
        if (class(niche.alphahull) != "try-error" & class(rangepcaarea.alphahull) != "try-error") {break} # if the tested alpha value results in a valid alpha shape, stop the loio    
      } # end kk loop
    } 
  }  
  
  if (MCP.niche !=  "Y"){
    # if neither the ahull nor the areaahull command fails, then get the "track" of the alpha shape
    niche.alpha.shape.track <- alphahull::ahull_track(x = mymat.df$x, y = mymat.df$y, nps= 1000, alpha=alpha.used)
    niche.all.alpha.coords <- NULL
    for (no_segments in 1:length(niche.alpha.shape.track)){ 
      niche.all.alpha.coords <-  rbind(niche.all.alpha.coords,niche.alpha.shape.track[[no_segments]]$data)}
    speciespcarange_i_xrange.alphahull <- range(niche.all.alpha.coords$x)
    speciespcarange_i_yrange.alphahull <-  range(niche.all.alpha.coords$y)
    # NB we do not confine this range area to terrestrial b/c alphahulls allow for holes in the species range
  }    
  
  # finally, construct the 95% convex hull around occurrences of this species in pca space 
  if (length(speciespoints_i) >= 5){
  rangepca.95 <- adehabitatHR::mcp(SpatialPoints(mymat), percent= 95, unin= "m", unout= "m2") # units on this are m^2
  speciespcarange_i_xrange.95 <-   as.numeric(attributes(rangepca.95)$bbox["x",])
  speciespcarange_i_yrange.95 <-   as.numeric(attributes(rangepca.95)$bbox["y",])
  rangepcaarea.95 <- sum(rangepca.95@data$area)
  }
  # final checks on potential range issues for this species
  if (rangearea== 0 | is.na(rangearea)){stop("range area is zero or na!")} 
  if(MCP.range== "N") if (identical(rangearea.alphahull, 0)){stop("alphahull range area is zero or na!")} 
  if (length(speciespoints_i) >= 5){  if (rangearea.95== 0 | is.na(rangearea.95)){stop("95% range area is zero or na!")} }
  if (rangepcaarea== 0 | is.na(rangepcaarea)){stop("range pca area is zero")}
  if(MCP.niche== "N") if (identical(rangepcaarea.alphahull,0 )){stop("range pca area is zero for alphahull")}
  if (length(speciespoints_i) >= 5){ if (rangepcaarea.95 == 0 | is.na(rangepcaarea.95)){stop("range pca area is zero for 95")}}
  
  
  for (ii in 1:dim(demo.speciesi)[1]){ # cycles through the demographic data points for a given species 
    
    yes <- demo.speciesi$Lat[ii]	# y of the demo datapoint in lat long
    xes <- demo.speciesi$Long[ii] # x of the demo datapoint in lat long
    
    demo.pt  <- SpatialPoints(coords= cbind(x=xes,y=yes), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) # making a Spatial points object of the location of a particualr demographic location ii for species i
    if (sp::over(demo.pt, speciesrange_i) !=1) {stop("demo point is not in species range")} # double- check that the demo.pt is inside the geographic range
    
    # finding categorical metrics of position within range-----
    # NB these categorical definitions of position within range are meaningless for species whose ranges overlap the equator; the code fixes that
    # in the section called "finding alternate metrics of position within range or niche for species whose ranges overlap the poles or equator"
    
    # categorical position within range for convex hull 
    yrange.lowerquantile <- min(speciesrange_i_absyrange)+ (max(speciesrange_i_absyrange)- min(speciesrange_i_absyrange))*cat.quantile
    yrange.upperquantile <- min(speciesrange_i_absyrange)+ (max(speciesrange_i_absyrange)- min(speciesrange_i_absyrange))*(1-cat.quantile)
    xrange.lowerquantile <- min(speciesrange_i_xrange)+ (max(speciesrange_i_xrange)- min(speciesrange_i_xrange))*cat.quantile
    xrange.upperquantile <- min(speciesrange_i_xrange)+ (max(speciesrange_i_xrange)- min(speciesrange_i_xrange))*(1-cat.quantile)
    
    if (abs(demo.pt@coords[2]) > yrange.upperquantile) {cat.descriptor <- "poleward.edge"} 
    if (abs(demo.pt@coords[2]) < yrange.lowerquantile) {cat.descriptor <- "equatorward.edge"}
    if (abs(demo.pt@coords[2]) <= yrange.upperquantile & abs(demo.pt@coords[2]) >= yrange.lowerquantile) {
      if (demo.pt@coords[1] < xrange.lowerquantile) {cat.descriptor <- "central.edge"} else if (
        demo.pt@coords[1] > xrange.upperquantile) {cat.descriptor <- "central.edge"} else if (
          demo.pt@coords[1] <= xrange.upperquantile & demo.pt@coords[1] >= xrange.lowerquantile) {cat.descriptor <- "central"} else {
          cat.descriptor <- NA}
    }
    
    # categorical position within range for 95% hull 
    if (length(speciespoints_i) >= 5){
    yrange.lowerquantile.95 <- min(speciesrange_i_absyrange.95)+ (max(speciesrange_i_absyrange.95)- min(speciesrange_i_absyrange.95))*cat.quantile
    yrange.upperquantile.95 <- min(speciesrange_i_absyrange.95)+ (max(speciesrange_i_absyrange.95)- min(speciesrange_i_absyrange.95))*(1-cat.quantile)
    xrange.lowerquantile.95 <- min(speciesrange_i_xrange.95)+ (max(speciesrange_i_xrange.95)- min(speciesrange_i_xrange.95))*cat.quantile
    xrange.upperquantile.95 <- min(speciesrange_i_xrange.95)+ (max(speciesrange_i_xrange.95)- min(speciesrange_i_xrange.95))*(1-cat.quantile)
    
    if (abs(demo.pt@coords[2]) > yrange.upperquantile.95) {cat.descriptor.95 <-  "poleward.edge"} 
    if (abs(demo.pt@coords[2]) < yrange.lowerquantile.95) {cat.descriptor.95 <- "equatorward.edge"}
    if (abs(demo.pt@coords[2]) <= yrange.upperquantile.95 & abs(demo.pt@coords[2]) >= yrange.lowerquantile.95) {
      if (demo.pt@coords[1] < xrange.lowerquantile.95) {cat.descriptor.95 <- "central.edge"} else if (
        demo.pt@coords[1] > xrange.upperquantile.95) {cat.descriptor.95 <- "central.edge"} else if (
          demo.pt@coords[1] <= xrange.upperquantile.95 & demo.pt@coords[1] >= xrange.lowerquantile.95) {cat.descriptor.95 <- "central"} else {
          cat.descriptor.95 <- NA}
    }}
    
    # categorial poisition within range for alphahull
    if (MCP.range == "N"){ # assuming we were able to calculate an alpha hull.... 
      yrange.lowerquantile.alphahull <- min(speciesrange_i_absyrange.alphahull)+ (max(speciesrange_i_absyrange.alphahull)- min(speciesrange_i_absyrange.alphahull))*cat.quantile
      yrange.upperquantile.alphahull <- min(speciesrange_i_absyrange.alphahull)+ (max(speciesrange_i_absyrange.alphahull)- min(speciesrange_i_absyrange.alphahull))*(1-cat.quantile)
      xrange.lowerquantile.alphahull <- min(speciesrange_i_xrange.alphahull)+ (max(speciesrange_i_xrange.alphahull)- min(speciesrange_i_xrange.alphahull))*cat.quantile
      xrange.upperquantile.alphahull <- min(speciesrange_i_xrange.alphahull)+ (max(speciesrange_i_xrange.alphahull)- min(speciesrange_i_xrange.alphahull))*(1-cat.quantile)
      if (abs(demo.pt@coords[2]) > yrange.upperquantile.alphahull) {cat.descriptor.alphahull <-  "poleward.edge"} 
      if (abs(demo.pt@coords[2]) < yrange.lowerquantile.alphahull) {cat.descriptor.alphahull <- "equatorward.edge"}
      if (abs(demo.pt@coords[2]) <= yrange.upperquantile.alphahull & abs(demo.pt@coords[2]) >= yrange.lowerquantile.alphahull) {
        if (demo.pt@coords[1] < xrange.lowerquantile.alphahull) {cat.descriptor.alphahull <- "central.edge"} else if (
          demo.pt@coords[1] > xrange.upperquantile.alphahull) {cat.descriptor.alphahull <- "central.edge"} else if (
            demo.pt@coords[1] <= xrange.upperquantile.alphahull & demo.pt@coords[1] >= xrange.lowerquantile.alphahull) {cat.descriptor.alphahull <- "central"} else {
            cat.descriptor.alphahull <- NA}
      } }
    
    # finding alternate metrics of position within range for species whose ranges overlap the poles or equator-----
    # checking whether the species overlaps the poles or the equator & overwriting the cat.descriptor if so!!
    if (!is.na(sp::over(S.pole.pt, speciesrange_i))) {cat.descriptor <- "overlaps S pole"
    stop("this species' range map overlaps S pole")} 
    if (!is.na(sp::over(N.pole.pt, speciesrange_i))) {cat.descriptor <- "overlaps N pole"
    stop("this species' range map overlaps N pole")} 
    # as of now, no species overlap N or S pole, but if they did, you'd need to add in catches for alphahull and 0.95 too
    
    if (min(speciesrange_i_yrange)<0 & max(speciesrange_i_yrange)>0) {cat.descriptor <- "overlaps equator"}
    if ( MCP.range== "N"){  if (min(speciesrange_i_yrange.alphahull)<0 & max(speciesrange_i_yrange.alphahull)>0) {cat.descriptor.alphahull <- "overlaps equator"}}
    if (length(speciespoints_i) >= 5){  if (min(speciesrange_i_yrange.95)<0 & max(speciesrange_i_yrange.95)>0) {cat.descriptor.95 <- "overlaps equator"}}
    
    # finding alternate metrics of position within range for species whose ranges overlap the poles or equator (convex hull)
    if (cat.descriptor== "overlaps equator"){
      if(demo.pt@coords[2]>=0 ){ # if the demo pop is in the Northern hemisphere or on the equator
        # calculate N hemisphere-specific quantiles
        yrange.lowerquantile.North <- 0+ (max(speciesrange_i_absyrange.North)- 0)*cat.quantile
        yrange.upperquantile.North <- 0+ (max(speciesrange_i_absyrange.North)-  0)*(1-cat.quantile)
        xrange.lowerquantile.North <- min(speciesrange_i_xrange.North)+ (max(speciesrange_i_xrange.North)-  min(speciesrange_i_xrange.North))*cat.quantile
        xrange.upperquantile.North <- min(speciesrange_i_xrange.North)+ (max(speciesrange_i_xrange.North)- min(speciesrange_i_xrange.North))*(1-cat.quantile)
        # and use the quantiles to get an alternate version of categorical range position, called cat.descriptor.eqsp
        # this descriptor assumes that equatorward edge populations are NOT on the edge, but near the equator
        if (abs(demo.pt@coords[2]) > yrange.upperquantile.North) {cat.descriptor.eqsp <- "poleward.edge"} 
        if (abs(demo.pt@coords[2]) < yrange.lowerquantile.North) {cat.descriptor.eqsp <- "equatorward.edge"}
        if (abs(demo.pt@coords[2]) <= yrange.upperquantile.North & abs(demo.pt@coords[2]) >= yrange.lowerquantile.North) {
          if (demo.pt@coords[1] < xrange.lowerquantile.North) {cat.descriptor.eqsp <- "central.edge"} else if (
            demo.pt@coords[1] > xrange.upperquantile.North) {cat.descriptor.eqsp <- "central.edge"} else if (
              demo.pt@coords[1] <= xrange.upperquantile.North & demo.pt@coords[1] >= xrange.lowerquantile.North) {
                cat.descriptor.eqsp <- "central"} else {
                  cat.descriptor.eqsp <- NA}}
      }
      
      if(demo.pt@coords[2]<0){ # if the demo pop is in the Southern hemisphere 
        # calculate S hemisphere-specific quantiles
        yrange.lowerquantile.South <- 0+ (max(speciesrange_i_absyrange.South)- 0)*cat.quantile
        yrange.upperquantile.South <- 0+ (max(speciesrange_i_absyrange.South)- 0)*(1-cat.quantile)
        xrange.lowerquantile.South <-  min(speciesrange_i_xrange.South)+ (max(speciesrange_i_xrange.South)- min(speciesrange_i_xrange.South))*cat.quantile
        xrange.upperquantile.South <- min(speciesrange_i_xrange.South)+ (max(speciesrange_i_xrange.South)- min(speciesrange_i_xrange.South))*(1-cat.quantile)
        # and use the quantiles to get an alternate version of categorical range position, called cat.descriptor.eqsp
        # this descriptor assumes that equatorward edge populations are NOT on the edge, but near the equator
        if (abs(demo.pt@coords[2]) > yrange.upperquantile.South) {cat.descriptor.eqsp <- "poleward.edge"} 
        if (abs(demo.pt@coords[2]) < yrange.lowerquantile.South) {cat.descriptor.eqsp <-  "equatorward.edge"}
        if (abs(demo.pt@coords[2]) <= yrange.upperquantile.South & abs(demo.pt@coords[2]) >= yrange.lowerquantile.South) {
          if (demo.pt@coords[1] < xrange.lowerquantile.South) {cat.descriptor.eqsp <- "central.edge"} else if (
            demo.pt@coords[1] > xrange.upperquantile.South) {cat.descriptor.eqsp <-  "central.edge"} else if (
              demo.pt@coords[1] <= xrange.upperquantile.South & demo.pt@coords[1] >= xrange.lowerquantile.South) {cat.descriptor.eqsp <-  "central"} else {
              cat.descriptor.eqsp <- NA}}
      }
    } else {cat.descriptor.eqsp  <- cat.descriptor}
    
    # finding alternate metrics of position within range for species whose ranges overlap the poles or equator (95% hull)
    if (length(speciespoints_i) >= 5){ if (cat.descriptor.95 == "overlaps equator"){
      
      
      if(demo.pt@coords[2]>=0 ){ # if the demo pop is in the NOrthern hemisphere or on the equator
        # calculate N hemisphere-specific quantiles
        yrange.lowerquantile.North.95 <- 0+ (max(speciesrange_i_absyrange.North.95)- 0)*cat.quantile
        yrange.upperquantile.North.95 <- 0+ (max(speciesrange_i_absyrange.North.95)-0)*(1-cat.quantile)
        xrange.lowerquantile.North.95 <- min(speciesrange_i_xrange.North.95)+ (max(speciesrange_i_xrange.North.95)-  min(speciesrange_i_xrange.North.95))*cat.quantile
        xrange.upperquantile.North.95 <- min(speciesrange_i_xrange.North.95)+ (max(speciesrange_i_xrange.North.95)- min(speciesrange_i_xrange.North.95))*(1-cat.quantile)
        # and use the quantiles to get an alternate version of categorical range position, called cat.descriptor.eqsp
        # this descriptor assumes that equatorward edge populations are NOT on the edge, but near the equator
        if (abs(demo.pt@coords[2]) > yrange.upperquantile.North.95) {cat.descriptor.eqsp.95 <- "poleward.edge"} 
        if (abs(demo.pt@coords[2]) < yrange.lowerquantile.North.95) {cat.descriptor.eqsp.95 <- "equatorward.edge"}
        if (abs(demo.pt@coords[2]) <= yrange.upperquantile.North.95 & abs(demo.pt@coords[2]) >= yrange.lowerquantile.North.95) {
          if (demo.pt@coords[1] < xrange.lowerquantile.North.95) {cat.descriptor.eqsp.95 <- "central.edge"} else if (
            demo.pt@coords[1] > xrange.upperquantile.North.95) {cat.descriptor.eqsp.95 <- "central.edge"} else if (
              demo.pt@coords[1] <= xrange.upperquantile.North.95 & demo.pt@coords[1] >= xrange.lowerquantile.North.95) {cat.descriptor.eqsp.95 <- "central"} else {
              cat.descriptor.eqsp.95 <- NA}}
      }
      
      if(demo.pt@coords[2]<0){ # if the demo pop is in the Southern hemisphere 
        # calculate S hemisphere-specific quantiles
        yrange.lowerquantile.South.95 <-0+ (max(speciesrange_i_absyrange.South.95)- 0)*cat.quantile
        yrange.upperquantile.South.95 <- 0+ (max(speciesrange_i_absyrange.South.95)- 0)*(1-cat.quantile)
        xrange.lowerquantile.South.95 <- min(speciesrange_i_xrange.South.95)+ (max(speciesrange_i_xrange.South.95)- min(speciesrange_i_xrange.South.95))*cat.quantile
        xrange.upperquantile.South.95 <- min(speciesrange_i_xrange.South.95)+ (max(speciesrange_i_xrange.South.95)- min(speciesrange_i_xrange.South.95))*(1-cat.quantile)
        # and use the quantiles to get an alternate version of categorical range position, called cat.descriptor.eqsp
        # this descriptor assumes that equatorward edge populations are NOT on the edge, but near the equator
        if (abs(demo.pt@coords[2]) > yrange.upperquantile.South.95) {cat.descriptor.eqsp.95 <- "poleward.edge"} 
        if (abs(demo.pt@coords[2]) < yrange.lowerquantile.South.95) {cat.descriptor.eqsp.95 <- "equatorward.edge"}
        if (abs(demo.pt@coords[2]) <= yrange.upperquantile.South.95 & abs(demo.pt@coords[2]) >= yrange.lowerquantile.South.95) {
          if (demo.pt@coords[1] < xrange.lowerquantile.South.95) {cat.descriptor.eqsp.95 <- "central.edge"} else if (
            demo.pt@coords[1] > xrange.upperquantile.South.95) {cat.descriptor.eqsp.95 <- "central.edge"} else if (
              demo.pt@coords[1] <= xrange.upperquantile.South.95 & demo.pt@coords[1] >= xrange.lowerquantile.South.95) {cat.descriptor.eqsp.95 <- "central"} else {
              cat.descriptor.eqsp.95 <- NA}}
      }
    } else {cat.descriptor.eqsp.95 <- cat.descriptor.95}}
    
    # finding alternate metrics of position within range for species whose ranges overlap the poles or equator (alpha shape)
    if ( MCP.range== "N"){
      if (cat.descriptor.alphahull == "overlaps equator"){
        
        if(demo.pt@coords[2]>=0 ){ # if the demo pop is in the NOrthern hemisphere or on the equator
          # calculate N hemisphere-specific quantiles
          yrange.lowerquantile.North.alphahull <- 0+ (max(speciesrange_i_absyrange.North.alphahull)- 0)*cat.quantile
          yrange.upperquantile.North.alphahull <- 0+ (max(speciesrange_i_absyrange.North.alphahull)- 0)*(1-cat.quantile)
          xrange.lowerquantile.North.alphahull <- min(speciesrange_i_xrange.North.alphahull)+ (max(speciesrange_i_xrange.North.alphahull)-  min(speciesrange_i_xrange.North.alphahull))*cat.quantile
          xrange.upperquantile.North.alphahull <- min(speciesrange_i_xrange.North.alphahull)+ (max(speciesrange_i_xrange.North.alphahull)- min(speciesrange_i_xrange.North.alphahull))*(1-cat.quantile)
          # and use the quantiles to get an alternate version of categorical range position, called cat.descriptor.eqsp
          # this descriptor assumes that equatorward edge populations are NOT on the edge, but near the equator
          if (abs(demo.pt@coords[2]) > yrange.upperquantile.North.alphahull) {cat.descriptor.eqsp.alphahull <- "poleward.edge"} 
          if (abs(demo.pt@coords[2]) < yrange.lowerquantile.North.alphahull) {cat.descriptor.eqsp.alphahull <- "equatorward.edge"}
          if (abs(demo.pt@coords[2]) <= yrange.upperquantile.North.alphahull & abs(demo.pt@coords[2]) >= yrange.lowerquantile.North.alphahull) {
            if (demo.pt@coords[1] < xrange.lowerquantile.North.alphahull) {cat.descriptor.eqsp.alphahull <- "central.edge"} else if (
              demo.pt@coords[1] > xrange.upperquantile.North.alphahull) {cat.descriptor.eqsp.alphahull <- "central.edge"} else if (
                demo.pt@coords[1] <= xrange.upperquantile.North.alphahull & demo.pt@coords[1] >= xrange.lowerquantile.North.alphahull) {cat.descriptor.eqsp.alphahull <- "central"} else {
                cat.descriptor.eqsp.alphahull <- NA}}
        }
        
        if(demo.pt@coords[2]<0){ # if the demo pop is in the Southern hemisphere 
          # calculate S hemisphere-specific quantiles
          yrange.lowerquantile.South.alphahull <-0+ (max(speciesrange_i_absyrange.South.alphahull)-0)*cat.quantile
          yrange.upperquantile.South.alphahull <- 0+ (max(speciesrange_i_absyrange.South.alphahull)- 0)*(1-cat.quantile)
          xrange.lowerquantile.South.alphahull <- min(speciesrange_i_xrange.South.alphahull)+ (max(speciesrange_i_xrange.South.alphahull)-  min(speciesrange_i_xrange.South.alphahull))*cat.quantile
          xrange.upperquantile.South.alphahull <- min(speciesrange_i_xrange.South.alphahull)+ (max(speciesrange_i_xrange.South.alphahull)- min(speciesrange_i_xrange.South.alphahull))*(1-cat.quantile)
          # and use the quantiles to get an alternate version of categorical range position, called cat.descriptor.eqsp
          # this descriptor assumes that equatorward edge populations are NOT on the edge, but near the equator
          if (abs(demo.pt@coords[2]) > yrange.upperquantile.South.alphahull) {cat.descriptor.eqsp.alphahull <- "poleward.edge"} 
          if (abs(demo.pt@coords[2]) < yrange.lowerquantile.South.alphahull) {cat.descriptor.eqsp.alphahull <- "equatorward.edge"}
          if (abs(demo.pt@coords[2]) <= yrange.upperquantile.South.alphahull & abs(demo.pt@coords[2]) >= yrange.lowerquantile.South.alphahull) {
            if (demo.pt@coords[1] < xrange.lowerquantile.South.alphahull) {cat.descriptor.eqsp.alphahull <- "central.edge"} else if (
             demo.pt@coords[1] > xrange.upperquantile.South.alphahull) {cat.descriptor.eqsp.alphahull <- "central.edge"} else if (
                demo.pt@coords[1] <= xrange.upperquantile.South.alphahull & demo.pt@coords[1] >= xrange.lowerquantile.South.alphahull) {cat.descriptor.eqsp.alphahull <- "central"} else {
                cat.descriptor.eqsp.alphahull <- NA}}
        }
        
      } else {cat.descriptor.eqsp.alphahull <- cat.descriptor.alphahull}} else {cat.descriptor.eqsp.alphahull <- cat.descriptor.eqsp }
    
    
    
    
    # continuous position within range ----
    
    # rangeposition.quantile is: what fraction of the abs latitudinal span does this point encompass
    
    # first, for convex hull
    rangeposition.quantile <- ifelse(cat.descriptor %in% c("overlaps S pole", "overlaps N pole", "overlaps equator"), NA, 
                                     (abs(demo.pt@coords[2])-min(speciesrange_i_absyrange))/(max(speciesrange_i_absyrange)-min(speciesrange_i_absyrange))) 
    rangeposition.quantile.eqsp <-  ifelse(cat.descriptor == "overlaps equator" & demo.pt@coords[2]>=0, #if a species overlaps the equator, and the demo point is on equator or in N hemisphere
                                           (abs(demo.pt@coords[2])-min(speciesrange_i_absyrange.North))/(max(speciesrange_i_absyrange.North)-min(speciesrange_i_absyrange.North)), 
                                           ifelse(cat.descriptor == "overlaps equator" & demo.pt@coords[2]<0, 
                                                  (abs(demo.pt@coords[2])-min(speciesrange_i_absyrange.South))/(max(speciesrange_i_absyrange.South)-min(speciesrange_i_absyrange.South))   , 
                                                  rangeposition.quantile))
    
    # next, for 95% hull 
    if (length(speciespoints_i) >= 5){
    rangeposition.quantile.95 <- ifelse(cat.descriptor.95 %in% c("overlaps S pole", "overlaps N pole", "overlaps equator"), NA, 
                                        (abs(demo.pt@coords[2])-min(speciesrange_i_absyrange.95))/(max(speciesrange_i_absyrange.95)-min(speciesrange_i_absyrange.95))) 
    rangeposition.quantile.eqsp.95  <- ifelse(cat.descriptor.95 == "overlaps equator" & demo.pt@coords[2]>=0, #if a species overlaps the equator, and the demo point is on equator/ in N hemisphere
                                              (abs(demo.pt@coords[2])-0)/(max(speciesrange_i_absyrange.North.95)-0), 
                                              ifelse(cat.descriptor.95 == "overlaps equator" & demo.pt@coords[2]<0, 
                                                     (abs(demo.pt@coords[2])-0)/(max(speciesrange_i_absyrange.South.95)-0)   , 
                                                     rangeposition.quantile.95))
    
    if(speciesrange_i_yrange.95[1]>0 & speciesrange_i_yrange.95[2]>0  & demo.pt@coords[2] <0){ # if the demo point is on the other side of the equator from the 95% range... 
      rangeposition.quantile.95 <-rangeposition.quantile.eqsp.95 <-cat.descriptor.eqsp.95 <- cat.descriptor.95 <-   NA }
    if(speciesrange_i_yrange.95[1]<0 & speciesrange_i_yrange.95[2]<0  & demo.pt@coords[2]>0){ 
      rangeposition.quantile.95 <-rangeposition.quantile.eqsp.95 <-cat.descriptor.eqsp.95 <- cat.descriptor.95 <-   NA }
    }
    # finally, for alpha shape
    if (MCP.range== "N") {
      rangeposition.quantile.alphahull <- ifelse(cat.descriptor.alphahull %in% c("overlaps S pole", "overlaps N pole", "overlaps equator"), NA, 
                                                 (abs(demo.pt@coords[2])-min(speciesrange_i_absyrange.alphahull))/
                                                   (max(speciesrange_i_absyrange.alphahull)-min(speciesrange_i_absyrange.alphahull)))
      
      rangeposition.quantile.eqsp.alphahull  <- 
        ifelse(cat.descriptor.alphahull == "overlaps equator" & demo.pt@coords[2]>=0, #if a species overlaps the equator, and the demo point is on equator/ in N hemisphere
               (abs(demo.pt@coords[2])-0)/(max(speciesrange_i_absyrange.North.alphahull)-0), 
               ifelse(cat.descriptor.alphahull == "overlaps equator" & demo.pt@coords[2]<0, 
                      (abs(demo.pt@coords[2])-0)/(max(speciesrange_i_absyrange.South.alphahull)-0)  , 
                      rangeposition.quantile.alphahull))
      
      if(speciesrange_i_yrange.alphahull[1]>0 & speciesrange_i_yrange.alphahull[2]>0  & demo.pt@coords[2] <0){ # if the demo point is on the other side of the equator from the alphahull
        rangeposition.quantile.alphahull <-rangeposition.quantile.eqsp.alphahull <-cat.descriptor.eqsp.alphahull <- cat.descriptor.alphahull <-   NA }
      if(speciesrange_i_yrange.alphahull[1]<0 & speciesrange_i_yrange.alphahull[2]<0  & demo.pt@coords[2]>0){
        rangeposition.quantile.alphahull <-rangeposition.quantile.eqsp.alphahull <-cat.descriptor.eqsp.alphahull <- cat.descriptor.alphahull <-   NA }
    } 
    
    # getting categorical position in niche ----- 
    
    # get the PCA score for the demographic point in question
    demo.pt.temp.precip <- raster::extract(r, demo.pt) 
    if (length(which(is.na(demo.pt.temp.precip)))>0){stop("the bioclim raster is returning NA for the demo pt; is it in the ocean??!! ")}  
    centeredNewData <-apply(demo.pt.temp.precip, MARGIN=1, FUN=c.fun, res$center, res$scale  ) # using PCA on the whole world
    demo_pca_i <-t(t(res$rotation) %*% centeredNewData)[,1:2]
    demopca.pt <- SpatialPoints(coords= t(as.matrix(demo_pca_i))) # putting the demographic point in the correct format (a Spatil Points object)
    # I do not double check that the demo point is in the pca range, because if it looks like it's out of the range,
    # it's actually on the edge (there is some rounding error that can happen in the pca calculation). 
    # note that the demo point was included in the estimates of convex hull so it should be in there
    
    # getting categorical position in niche for convex hull
    yrange.lowerquantile <- min(speciespcarange_i_yrange)+ (max(speciespcarange_i_yrange)- min(speciespcarange_i_yrange))*cat.quantile # this is pca 2
    yrange.upperquantile <- min(speciespcarange_i_yrange)+ (max(speciespcarange_i_yrange)- min(speciespcarange_i_yrange))*(1-cat.quantile)
    xrange.lowerquantile <- min(speciespcarange_i_xrange)+ (max(speciespcarange_i_xrange)- min(speciespcarange_i_xrange))*cat.quantile # this is pca 1
    xrange.upperquantile <- min(speciespcarange_i_xrange)+ (max(speciespcarange_i_xrange)- min(speciespcarange_i_xrange))*(1-cat.quantile)
    
    if (demopca.pt@coords[1] < xrange.lowerquantile) {cat.pca.descriptor  <- "low.PCA1.score"}
    if (demopca.pt@coords[1] > xrange.upperquantile) {cat.pca.descriptor <- "high.PCA1.score"}
    if (demopca.pt@coords[1] <= xrange.upperquantile & demopca.pt@coords[1] >= xrange.lowerquantile) {
      if (demopca.pt@coords[2] < yrange.lowerquantile) {cat.pca.descriptor <-  "extreme.PCA2.score"} else if (
        demopca.pt@coords[2] > yrange.upperquantile) {cat.pca.descriptor <- "extreme.PCA2.score"} else if (
          demopca.pt@coords[2] <= yrange.upperquantile & demopca.pt@coords[2] >= yrange.lowerquantile) {cat.pca.descriptor <-  "central.PCA.score"} else {
          cat.pca.descriptor <- NA}
      
    }
    
    # getting categorical position in niche for 95% hull
    if (length(speciespoints_i) >= 5){
    yrange.lowerquantile.95 <- min(speciespcarange_i_yrange.95)+ (max(speciespcarange_i_yrange.95)- min(speciespcarange_i_yrange.95))*cat.quantile # this is pca 2
    yrange.upperquantile.95 <- min(speciespcarange_i_yrange.95)+ (max(speciespcarange_i_yrange.95)- min(speciespcarange_i_yrange.95))*(1-cat.quantile)
    xrange.lowerquantile.95 <- min(speciespcarange_i_xrange.95)+ (max(speciespcarange_i_xrange.95)- min(speciespcarange_i_xrange.95))*cat.quantile # this is pca 1
    xrange.upperquantile.95 <- min(speciespcarange_i_xrange.95)+ (max(speciespcarange_i_xrange.95)- min(speciespcarange_i_xrange.95))*(1-cat.quantile)
    
    if (demopca.pt@coords[1] < xrange.lowerquantile.95) {cat.pca.descriptor.95 <- "low.PCA1.score"}
    if (demopca.pt@coords[1] > xrange.upperquantile.95) {cat.pca.descriptor.95 <- "high.PCA1.score"}
    if (demopca.pt@coords[1] <= xrange.upperquantile.95 & demopca.pt@coords[1] >= xrange.lowerquantile.95) {
      if (demopca.pt@coords[2] < yrange.lowerquantile.95) {cat.pca.descriptor.95<- "extreme.PCA2.score"} else if (
        demopca.pt@coords[2] > yrange.upperquantile.95) {cat.pca.descriptor.95<- "extreme.PCA2.score"} else if (
          demopca.pt@coords[2] <= yrange.upperquantile.95 & demopca.pt@coords[2] >= yrange.lowerquantile.95) {
            cat.pca.descriptor.95 <- "central.PCA.score"} else {
              cat.pca.descriptor.95<- NA} }}
    
    # getting categorical position in niche for alpha shape  
    if (MCP.niche== "N"){
      yrange.lowerquantile.alphahull <- min(speciespcarange_i_yrange.alphahull)+ (max(speciespcarange_i_yrange.alphahull)- min(speciespcarange_i_yrange.alphahull))*cat.quantile # this is pca 2
      yrange.upperquantile.alphahull <- min(speciespcarange_i_yrange.alphahull)+ (max(speciespcarange_i_yrange.alphahull)- min(speciespcarange_i_yrange.alphahull))*(1-cat.quantile)
      xrange.lowerquantile.alphahull <- min(speciespcarange_i_xrange.alphahull)+ (max(speciespcarange_i_xrange.alphahull)- min(speciespcarange_i_xrange.alphahull))*cat.quantile # this is pca 1
      xrange.upperquantile.alphahull <- min(speciespcarange_i_xrange.alphahull)+ (max(speciespcarange_i_xrange.alphahull)- min(speciespcarange_i_xrange.alphahull))*(1-cat.quantile)
      
      if (demopca.pt@coords[1] < xrange.lowerquantile.alphahull) {cat.pca.descriptor.alphahull  <- "low.PCA1.score"}
      if (demopca.pt@coords[1] > xrange.upperquantile.alphahull) {cat.pca.descriptor.alphahull <- "high.PCA1.score"}
      if (demopca.pt@coords[1] <= xrange.upperquantile.alphahull & demopca.pt@coords[1] >= xrange.lowerquantile.alphahull) {
        if (demopca.pt@coords[2] < yrange.lowerquantile.alphahull) {cat.pca.descriptor.alphahull <-  "extreme.PCA2.score"} else if (
          demopca.pt@coords[2] > yrange.upperquantile.alphahull) {cat.pca.descriptor.alphahull <- "extreme.PCA2.score"} else if (
            demopca.pt@coords[2] <= yrange.upperquantile.alphahull & demopca.pt@coords[2] >= yrange.lowerquantile.alphahull) {
              cat.pca.descriptor.alphahull <-  "central.PCA.score"} else {
                cat.pca.descriptor.alphahull <- NA}
        
      }} else  {cat.pca.descriptor.alphahull <- cat.pca.descriptor}
    
    
    
    # continuous position in niche----
    
    rangepcaposition.quantile <- (demopca.pt@coords[1]- min(speciespcarange_i_xrange))/ (max(speciespcarange_i_xrange)- min(speciespcarange_i_xrange))
    if (length(speciespoints_i) >= 5){rangepcaposition.quantile.95 <- (demopca.pt@coords[1]- min(speciespcarange_i_xrange.95))/ (max(speciespcarange_i_xrange.95)- min(speciespcarange_i_xrange.95))}
    if (MCP.niche == "N") {rangepcaposition.quantile.alphahull <- (demopca.pt@coords[1]- min(speciespcarange_i_xrange.alphahull))/ (max(speciespcarange_i_xrange.alphahull)- min(speciespcarange_i_xrange.alphahull))
    } else {rangepcaposition.quantile.alphahull <- rangepcaposition.quantile}
    
    # fixing variables if you cannot fit an alpha shape---
    if (MCP.range == "Y") { # some of these in the below are repeated from the code above
      rangearea.alphahull <- rangearea
      cat.descriptor.alphahull <- cat.descriptor
      cat.descriptor.eqsp.alphahull <- cat.descriptor.eqsp
      rangeposition.quantile.alphahull <- rangeposition.quantile
      rangeposition.quantile.eqsp.alphahull <- rangeposition.quantile.eqsp}
    
    if (MCP.niche == "Y") { # some of these in the below are repeated from the code above
      rangepcaarea.alphahull <- rangepcaarea
      cat.pca.descriptor.alphahull <- cat.pca.descriptor
      rangepcaposition.quantile.alphahull <- rangepcaposition.quantile}
    if (length(speciespoints_i) < 5){rangearea.95 <- rangepcaarea.95 <- cat.descriptor.95 <- 
      cat.descriptor.eqsp.95 <- rangeposition.quantile.95 <- rangeposition.quantile.eqsp.95<- cat.pca.descriptor.95 <- 
      rangepcaposition.quantile.95 <- NA}
    
    results[ii,] <- c(rangearea, 
                      rangepcaarea, 
                      cat.descriptor, 
                      cat.descriptor.eqsp, 
                      rangeposition.quantile, 
                      rangeposition.quantile.eqsp, 
                      cat.pca.descriptor,
                      rangepcaposition.quantile, 
                      
                      rangearea.95, 
                      rangepcaarea.95,  
                      cat.descriptor.95,  
                      cat.descriptor.eqsp.95, 
                      rangeposition.quantile.95,  
                      rangeposition.quantile.eqsp.95, 
                      cat.pca.descriptor.95, 
                      rangepcaposition.quantile.95, 
                      
                      rangearea.alphahull, 
                      rangepcaarea.alphahull,  
                      cat.descriptor.alphahull,  
                      cat.descriptor.eqsp.alphahull, 
                      rangeposition.quantile.alphahull,  
                      rangeposition.quantile.eqsp.alphahull, 
                      cat.pca.descriptor.alphahull, 
                      rangepcaposition.quantile.alphahull)
    
    rm(cat.descriptor, 
       cat.descriptor.eqsp, 
       rangeposition.quantile, 
       rangeposition.quantile.eqsp, 
       cat.pca.descriptor,
       rangepcaposition.quantile, 
       cat.descriptor.95,  
       cat.descriptor.eqsp.95, 
       rangeposition.quantile.95,  
       rangeposition.quantile.eqsp.95, 
       cat.pca.descriptor.95, 
       rangepcaposition.quantile.95, 
       cat.descriptor.alphahull,  
       cat.descriptor.eqsp.alphahull, 
       rangeposition.quantile.alphahull,  
       rangeposition.quantile.eqsp.alphahull, 
       cat.pca.descriptor.alphahull, 
       rangepcaposition.quantile.alphahull)
    
    
    
    
  } # end of ii loop
  species_sites_positions[species.indices,which(names(species_sites_positions)== "rangearea"): length(names(species_sites_positions))] <- results
  # get rid of most things in order to speed up the loop
  gdata::keep(alpharange, alphaniche, species_sites_positions, species_sites,
              i, m1, extract_coords1, S.pole.pt, N.pole.pt, 
              c.fun, r, res, terrestrial_sf, cat.quantile, sure= TRUE)    
  print(i)}

unique(species_sites_positions$cat.descriptor) # should be, um, more than just NA's
save.image(file= paste("data/02_positions.", format(Sys.Date(), format="%y%m%d"), ".cat.quantile",cat.quantile,".RData", sep= "" ))
