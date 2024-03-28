#### PROJECT: International Drought Experiment x Range Position
#### PURPOSE: Summarize how many species are observed across multiple sites and years
#### AUTHOR: Amy Angert
#### DATE LAST MODIFIED: 20240328

# Easy code for installing packages in R (if not installed) and calling their libraries
# From: https://gist.github.com/DrK-Lo/a945a29d6606b899022d0f03109b9483

# make vector of packages needed
packages_needed <- c("tidyverse")

# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}


# Read in data
dat <- read.csv("data/cover_ppt_2024-02-12.csv") # ~56k rows

# Exclude introduced species
table(dat$local_provenance)
dat <- dat %>% 
  filter(local_provenance=="NAT") %>% droplevels() # ~44k rows; removed ~12k observations whose provenance is introduced or unknown
  
# Exclude woody species, bryophytes, and lichens
table(dat$functional_group)
dat <- dat %>% 
  filter(functional_group!="WOODY") %>% 
  filter(functional_group!="BRYOPHYTE") %>% 
  filter(functional_group!="LICHEN") %>% droplevels() # ~38k rows; removed ~5k rows of woody species + ~1k rows of bryophytes and lichens

# Quick counts
n_sites <- n_distinct(dat$site_code) # 106 sites
n_years <- n_distinct(dat$year) # 12 years
n_species <- n_distinct(dat$Taxon) # 1697 species


# Summarize sampling by species
dat_summary <- dat %>% 
  group_by(Taxon, site_code) %>% 
  summarise(n_years = n_distinct(year)) %>% 
  pivot_wider(names_from=site_code, values_from=n_years) %>% # rearrange to species as rows, sites as columns
  ungroup() %>% 
  mutate(n_sites = rowSums(!is.na(across(-Taxon)))) %>% # number of sites per species
  mutate(ave_years = rowMeans(across(-c(Taxon,n_sites)),na.rm=TRUE)) %>% # mean years per site per species
  group_by(Taxon) %>% 
  mutate(max_years = max(across(-c(n_sites,ave_years)),na.rm=TRUE)) # max years in any site per sites

summary(dat_summary$n_sites)
summary(dat_summary$ave_years)
summary(dat_summary$max_years)

# Filter to species with >1 site and >1 year
dat_summary_cull <- dat_summary %>% 
  filter(n_sites>1) %>% # drops 1233 species
  filter(max_years>1) # drops 27 species

summary_cull_simple <- dat_summary_cull %>% select(Taxon, n_sites, max_years, ave_years) 

dat_cull <- left_join(summary_cull_simple, dat) # lose ~14k rows compared to full set

# Runs of consecutive years by species per treatment per site
dat_consec <- dat_cull %>% 
  group_by(Taxon, site_code, trt) %>% 
  arrange(year, .by_group = TRUE) %>% 
  mutate(diff = c(0, diff(year))) %>% # diff=1 every time row switches to next year
  filter(diff<2) %>% # get rid of observations across non-consecutive years
  mutate(consec = cumsum(diff)) %>% 
  summarise(max_consec = max(consec)) %>% 
  ungroup()  #2430 rows

dat_consec_cull <- dat_consec %>% 
  filter(max_consec>0) %>%  # get rid of rows where either treatment has 0 consecutive years (1745 rows remaining)
  group_by(Taxon, site_code) %>% 
  summarise(paired = length(max_consec)) %>% 
  filter(paired>1)  %>% droplevels() # get rid of sites that no longer have matched droughts and controls
  # 731 rows
  
dat_consec_cull <- left_join(dat_consec, dat_consec_cull) %>% 
  filter(paired==2) %>% #1462 rows (= 731 sites x 2 treatment levels)
  mutate(site.trt = paste(site_code, "_", trt, sep="")) %>% 
  select(-c(site_code, trt, paired)) %>% 
  pivot_wider(names_from=site.trt, values_from=max_consec) %>% 
  mutate(n_sitetrt = rowSums(!is.na(across(-Taxon)))) %>% # number of site x treatment combos per species
  filter(n_sitetrt > 2) # get rid of species that have only 1 site remaining --> 199 species remaining

summary(dat_consec_cull$n_sitetrt)  

write_csv(dat_consec_cull, "data/species_site_consecyears.csv")
  
