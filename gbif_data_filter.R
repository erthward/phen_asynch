

# curating GBIF dataset

######################################################################

# search parameters

# enough observations (include)
  # phylum: Tracheophyta
      # class: 
          #Magnoliopsida
          #Liliopsida
              # species: Selaginella selaginoides 
          #Polypodiopsida
              # species: Azolla filiculoides Lam., Osmunda regalis L., Botrychium lanceolatum, 
              # Botrychium matricariifolium, Sceptridium multifidum, Ophioglossum vulgatum L., 
              # Botrychium lunaria, Equisetum hyemale L., Equisetum variegatum Schleich., 
              # Equisetum sylvaticum L., Equisetum telmateia Ehrh., Equisetum fluviatile L.,
              # Equisetum palustre L., Equisetum arvense L.,
              # LEFT OFF ON Polypodiales 
          #Pinopsida
              # species: Picea sitchensis, Larix decidua Mill., Picea abies, Taxus baccata, Pinus sylvestris, Juniperus communis
          #Lycopodiopsida
              # species: Lycopodiella inundata, Lycopodium clavatum, Spinulum annotinum, Huperzia selago

# species with non-annual phenology
  # genus
    # ficus 

non_ann_pheno <- c()

# list of genuses that are domesticated
domesticated_genus <- c("Pyrus", "Malus", "Cydonia" )

# list of domesticated species
domesticated_species <- c(
  "Eriobotrya japonica", "Mespilus germanica", "Citrus medica", "Citrus paradisi",
  "Citrus limon", "Citrus maxima" 
)

# curating gBIF dataset 

library(rgbif)

dat <- occ_data(hasCoordinate = T, hasGeospatialIssue = T, establishmentMeans = "NATIVE", limit = 20) # 270,552

# get all plant occurrences with coordinates
plantae_key <- name_backbone(name = "Plantae", rank = "kingdom")$usageKey
dat <- occ_data(hasCoordinate = T, hasGeospatialIssue = T, kingdomKey = plantae_key) # 1,530,552

# extract their column names for reference 
gbif_dat_col_names <- colnames(dat$data)

# get the data frame of all plant occurrences with coordinates
gbif_dat <- dat$data
# create a filter to exclude all observations without coordinate uncertainty estimates 
coord_filter1 <- !is.na(gbif_dat$coordinateUncertaintyInMeters)
# filter out all observations without coordinate uncertainty estimates 
gbif_dat <- gbif_dat[coord_filter1, ]
# create a second filter to exclude all observations with coordinate uncertainty estimates < 10,000 meters 
coord_filter2 <- coord_filtered_dat1 <= 10000
# filter out all observations with coordinate uncertainty estimates < 10,000 meters 
gbif_dat <- gbif_dat[coord_filter2, ]


#sum(is.na(gbif_dat$lifeStage))

# rgbif doesn't contain all occurrences of plants
# next steps: download data directly from gBIF and manipulate in R

###################################################

# DOWNLOADED DATA DIRECTLY FROM gBIF
library(magrittr)
library(dplyr)
library(ggplot2)
# load a subset of the data: Tracheophyta > Magnoliopsida > Ranunculales, Brassicales, Rosales, Malpighiales, Gentianales
setwd("~/Downloads")
gbif_dat_subset1 <- read.csv("gBIFplant_subset1.csv", header = T, sep = "\t")

# function returns a filtered data set of gBIF data 
gBIF_data_filter = function(data, numSpeciesMin, coordUncertainty){
  # data = gBIF dataset
  # numSpeciesMin = minimum sample size for all the observations of a particular species 
  # coordUncertainty = exclude all observations with coordinate uncertainty of X meters
  
  # create a filter to remove all observations without longitudinal coordinates
  coord_filter1 <- !is.na(data$decimalLongitude)
  # filter out all observations without longitudinal coordinates 
  gbif_dat <- data[coord_filter1, ]
  # create a filter to remove all observations without latitudinal coordinates
  coord_filter2 <- !is.na(gbif_dat$decimalLatitude)
  # filter out all observations without latitudinal coordinates 
  gbif_dat1 <- gbif_dat[coord_filter2, ]
  # create a filter to exclude all observations without coordinate uncertainty estimates 
  coord_filter3 <- !is.na(gbif_dat1$coordinateUncertaintyInMeters)
  # filter out all observations without coordinate uncertainty estimates 
  gbif_dat2 <- gbif_dat1[coord_filter3, ]
  # create a second filter to exclude all observations with coordinate uncertainty estimates < coordUncertainty meters 
  coord_filter4 <- gbif_dat2$coordinateUncertaintyInMeters <= coordUncertainty
  # filter out all observations with coordinate uncertainty estimates < coordUncertainty meters 
  gbif_dat3 <- gbif_dat2[coord_filter4, ] 
  # find the frequency of occurrence for all the species in the table 
  species_counts <- gbif_dat3 %>% count(as.factor(scientificName))
  # create a mask to filter out all species with fewer than numSpeciesMin observations 
  min_species_filter  <- which(species_counts$n >= numSpeciesMin) 
  filtered_species_names <- species_counts[min_species_filter, 1]
  # filter out all species with fewer than numSpeciesMin observations
  all_species = gbif_dat3$scientificName
  num_species_filter = filtered_species_names
  species_filter <-  all_species %in% num_species_filter
  gbif_dat4 <- gbif_dat3[species_filter, ]
}

# compile the entire data set

# load a subset of the data: Tracheophyta > Magnoliopsida > Ranunculales, Brassicales, Rosales, Malpighiales, Gentianales
filtered_gbif_data1 <- gBIF_data_filter(data = gbif_dat_subset1, numSpeciesMin = 30, coordUncertainty = 10000)
#nrow(filtered_gbif_data1) # 413917
#length(unique(filtered_gbif_data1$scientificName)) # 679

# load a subset of the data: Tracheophyta > Magnoliopsida > Asterales, Caryophyllales
gbif_dat_subset2 <- read.csv("gBIFplant_subset2.csv", header = T, sep = "\t", quote = "")

filtered_gbif_data2 <- gBIF_data_filter(data = gbif_dat_subset2, numSpeciesMin = 30, coordUncertainty = 10000)
#nrow(filtered_gbif_data2) #   465388
#length(unique(filtered_gbif_data2$scientificName)) # 766

# load a subset of the data: Tracheophyta > Magnoliopsida > Lamiales, Fabales, Ericales, Apiales 
gbif_dat_subset3 <- read.csv("gBIFplant_subset3.csv", header = T, sep = "\t")

filtered_gbif_data3 <- gBIF_data_filter(data = gbif_dat_subset3, numSpeciesMin = 30, coordUncertainty = 10000)
#nrow(filtered_gbif_data3) #  273920
#length(unique(filtered_gbif_data3$scientificName)) # 549

# load a subset of everything else except Tracheophyta > Liliopsida
gbif_dat_subset4 <- read.csv("gBIFplant_subset4.csv", header = T, sep = "\t")

filtered_gbif_data4 <- gBIF_data_filter(data = gbif_dat_subset4, numSpeciesMin = 30, coordUncertainty = 10000)
#nrow(filtered_gbif_data4) # 232270
#length(unique(filtered_gbif_data4$scientificName)) # 477

# load a dataset only including Tracheophyta > Liliopsida
gbif_dat_subset5 <- read.csv("gBIFplant_subset5.csv", header = T, sep = "\t")

filtered_gbif_data5 <- gBIF_data_filter(data = gbif_dat_subset5, numSpeciesMin = 30, coordUncertainty = 10000)
#nrow(filtered_gbif_data5) # 71045
#length(unique(filtered_gbif_data5$scientificName)) # 395

# bind all the 5 subsets by their rows and create one massive dataset
gbif_plant <- bind_rows(filtered_gbif_data1, filtered_gbif_data2, 
                        filtered_gbif_data3, filtered_gbif_data4, filtered_gbif_data5)
nrow(gbif_plant) # 1,504,411
length(unique(gbif_plant$scientificName)) # 2,955

############################

## TEST CODE FOR FUNCTION ## 
# create a filter to exclude all observations without coordinate uncertainty estimates 
coord_filter1 <- !is.na(gbif_dat_subset1$coordinateUncertaintyInMeters)
# filter out all observations without coordinate uncertainty estimates 
gbif_dat <- gbif_dat_subset1[coord_filter1, ]
# create a second filter to exclude all observations with coordinate uncertainty estimates < 10,000 meters 
coord_filter2 <- gbif_dat$coordinateUncertaintyInMeters <= 10000
# filter out all observations with coordinate uncertainty estimates < 10,000 meters 
gbif_dat1 <- gbif_dat_subset1[coord_filter2, ] # 467103 results
# find the frequency of occurence for all the species in the table 
species_freq <- table(gbif_dat1$scientificName)
# create a mask to filter out all species with fewer than 30 observations 
species_30plus <- species_freq > 29
species_30plus <- rownames(species_30plus)
# filter out all species with fewer than 30 observations
all_species = gbif_dat1$scientificName
species_30 = species_30plus
species_filter <-  all_species %in% species_30
gbif_dat1 <- gbif_dat1[species_filter, ]

species_counts <- gbif_plant %>% count(as.factor(scientificName))
species_counts30 <- which(species_counts$n >= 30) # 1456540 total observations, 2866 species
species_30names <- species_counts[species_counts30, 1]

all_species = gbif_plant$scientificName
species_30 = species_30names
species_filter <-  all_species %in% species_30
gbif_dat1 <- gbif_plant[species_filter, ]

nrow(gbif_dat1) # 467103
length(unique(gbif_dat1$scientificName)) # 9023
############################

#### MAKE THE GRAPHICS FOR THE DATASET ###

# first: compile a table with the variance of longitude and latitude for each species
species_sd <- gbif_plant %>% 
  group_by(scientificName) %>% 
  summarise(species_sd_lon = sd(decimalLongitude), 
            species_sd_lat = sd(decimalLatitude)) 

# second: find the mean latitude for each species and add column to the table above 
  # mean latitude intended to be used for color-coding the points on the data visualization
species_lat_avg <- gbif_plant %>% 
  group_by(scientificName) %>% 
  summarise(species_lat_mean = mean(decimalLatitude)) 

species_data <- cbind(species_sd, species_lat_avg[ ,2])

nrow(species_data) #2,955 --> this checks out
sum(is.na(species_data$species_sd_lon)) # 0 
sum(is.na(species_data$species_sd_lat)) # 0 
sum(is.na(species_data$species_lat_mean)) # 0 

# plot the data
  # color-code the species data point in the plot as: red = equator, blue = poles 
ggplot(species_data, aes(x = species_sd_lon, y = species_sd_lat)) + 
  geom_point(aes(color = species_lat_mean), alpha = 0.7) +
  labs(title = "Species Latitudinal Standard Deviation vs Species Longitudinal Standard Deviation") + 
  xlab(label = "Species Latitudinal Standard Deviation") + 
  ylab(label = "Species Longitudinal Standard Deviation") + 
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_color_gradient2(low = "red", high = "blue", mid = "purple", 
                        midpoint = 42.25, name = "Mean Species Latitude") 
  

#The plot is intended to examine the geography variability of all observations 
# Note: low variance = data points are very close to one another.

# pick several species with large variance

# get the species by those with low latitudinal variance but high longitudinal variance
# first get species with large longitudinal variance
high_species_lon_sd <- sort(species_data$species_sd_lon, decreasing = T, index.return = TRUE)
high_lon_sd_name_filter <- high_species_lon_sd$ix[1:10]
high_lon_sd_species <- species_data$scientificName[high_lon_sd_name_filter]
# make sure these species also don't have large latitudinal sd  
lat_sd_for_high_lon <- species_data$species_sd_lat[high_lon_sd_name_filter]
lat_sd_for_high_lon_filter <- lat_sd_for_high_lon < 31

#sort the species by those that have high within species latitudinal variance
high_lat_species_sd <- sort(species_data$species_sd_lat, decreasing = T, index.return = TRUE)
# gather the names of the first 5 of these species
high_lat_sd_name_filter <- high_lat_species_sd$ix[1:5]
# filter out all species except for the selected 5  
high_lat_sd_species_names <- species_data$scientificName[high_lat_sd_name_filter]
# compile the dataframe of all information on these 5 species  
species_names <- c(high_lon_sd_species[1:6], high_lat_sd_species_names)
subset_species_data <- gbif_plant %>% filter(scientificName == species_names)

# species in table: "Tridax procumbens L." ,"Argemone mexicana L.", "Calotropis procera (Aiton) Aiton fil.",
#"Tribulus terrestris L.", "Reichardia tingitana (L.) Roth", "Mesembryanthemum nodiflorum L.", 
#"Carpobrotus edulis (L.) N.E.Br.", "Leycesteria formosa Wall.", "Oxalis articulata Savigny", "


#head(subset_species_data [1:5, 10:13])
#length(unique(subset_species_data$scientificName)) # 5




##########################################################
# other species used 
Boraginales 48,502
Myrtales 45,103
Solanales 35,104
Dipsacales 34,541
Geraniales 33,777
Saxifragales 27,043
Malvales 24,800
Piperales 9,305
Sapindales 8,700
Celastrales 6,947
Oxalidales 6,798
Nymphaeales 6,125
Fagales 5,679
Cucurbitales 5,012
Laurales 3,618
Cornales 2,980
Santalales 2,563
Magnoliales 2,158
Vitales 794
Aquifoliales 748
Proteales 
Zygophyllales 
Chloranthales 
Dilleniales 
Picramniales 
Canellales 
Crossosomatales 
Buxales 
Icacinales
Huerteales 
Gunnerales
Escalloniales
Metteniusales
Ceratophyllales
Bruniales
Garryales 
Austrobaileyales 

Tracheophyta
#Liliopsida 70,566,378
  Polypodiopsida 10,538,886
  Pinopsida 4,575,925
  Lycopodiopsida 1,016,562
  Gnetopsida 61,017
  Cycadopsida 56,165
  Ginkgoopsida 28,222

phylum 
  Bryophyta 10,924,267
  Marchantiophyta 2,472,518
  Rhodophyta 1,811,833
  Chlorophyta 1,339,179
  Charophyta 326,552
  Anthocerotophyta 121,995
  Glaucophyta 2,159
