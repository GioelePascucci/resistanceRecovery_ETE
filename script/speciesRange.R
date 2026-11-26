# ==============================================================================
# Script title: SpeciesRange.R
# Description:  The script integrates spatial data from multiple sources 
#               (IUCN range maps, GBIF occurrences, and WorldClim 
#               Digital Elevation Models) to calculate species distribution range.
# Manuscript title: Correlated evolution of thermal tolerance and egg size 
#                links resistance and recovery to extreme heat in amphibians.
# Manuscript authors: Gioele Pascucci, Miguel Tejedo, Lumír Gvoždík
# Script author: Gioele Pascucci
# Contact:  gioelepascucci@live.it
# Date: August 30, 2026
# R Version:  4.5.3
# ==============================================================================

rm(list = ls()) # Clean workspace

library(rgbif)             # Access GBIF data via DOI
library(CoordinateCleaner) # Clean occurrence coordinates
library(dplyr)
library(tidyr)
library(sf)                
library(terra)             
library(stringr)      
library(readxl)      
library(spData)           # Provides the 'us_states' dataset for spatial cropping

# ------------------------------------------------------------------------------
# Load and prepare IUCN shapefiles and trait dataset
# ------------------------------------------------------------------------------

# Load trait dataset, adjust to your path
main_data <- read_excel("../data/datasetTraits.xlsx")
# Extract and standardize species names from the database
species_original <- str_to_sentence(str_squish(main_data$Species))

# Unzip the IUCN shapefiles into /data directory (Download link: https://www.iucnredlist.org/resources/spatial-data-download#amphibians)
unzip("../data/AMPHIBIANS.zip", exdir = "../data/IUCN_temp")
# Read the parts of the IUCN database
iucn_part1 <- st_read("../data/IUCN_temp/spatial_data_1.shp", quiet = TRUE)
iucn_part2 <- st_read("../data/IUCN_temp/spatial_data_2.shp", quiet = TRUE)
iucn_part3 <- st_read("../data/IUCN_temp/spatial_data_3.shp", quiet = TRUE)
iucn_part4 <- st_read("../data/IUCN_temp/spatial_data_4.shp", quiet = TRUE)
iucn_part5 <- st_read("../data/IUCN_temp/spatial_data_5.shp", quiet = TRUE)
iucn_part6 <- st_read("../data/IUCN_temp/spatial_data_6.shp", quiet = TRUE)
iucn_part7 <- st_read("../data/IUCN_temp/spatial_data_7.shp", quiet = TRUE)
iucn_part8 <- st_read("../data/IUCN_temp/spatial_data_8.shp", quiet = TRUE)
# Merge into a single spatial dataframe
iucn_master <- bind_rows(iucn_part1, iucn_part2,iucn_part3,iucn_part4,iucn_part5,iucn_part6,iucn_part7,iucn_part8)
rm(iucn_part1, iucn_part2,iucn_part3,iucn_part4,iucn_part5,iucn_part6,iucn_part7,iucn_part8)

# Split P. glutinosus complex distribution range (Based on literure; AmphibiaWeb, 2026)
print("Splitting P. glutinosus and P. grobmani ranges...")

# Load US states and transform to match IUCN CRS
data("us_states", package = "spData")
us_states <- st_transform(us_states, st_crs(iucn_master))

# Define P. grobmani range: FL, GA, AL (Southeastern Coastal Plain, south of 32.5° N latitude)
al_ga_fl <- us_states %>% filter(NAME %in% c("Florida", "Georgia", "Alabama"))
southern_mask <- st_crop(al_ga_fl, xmin = -90, xmax = -75, ymin = 24, ymax = 32.5)

# extract original P. glutinosus polygon
lumped_poly <- iucn_master %>% filter(sci_name == "Plethodon glutinosus")

# Extract P. grobmani, overlap of the original polygon and the southern mask
poly_grobmani <- st_intersection(lumped_poly, st_union(southern_mask)) %>% 
  mutate(sci_name = "Plethodon grobmani")

# Extract new P. glutinosus, spatial difference between the original polygon and the southern mask
poly_glutinosus <- st_difference(lumped_poly, st_union(southern_mask)) %>% 
  mutate(sci_name = "Plethodon glutinosus")

# Remove the old P. glutinosus polygon and append the new ones
iucn_master <- iucn_master %>%
  filter(sci_name != "Plethodon glutinosus") %>%
  bind_rows(poly_glutinosus)

# Align taxonomy to Frost (2026)
print("Updating taxonomy to Frost (2026)...")

iucn_master <- iucn_master %>%
  mutate(sci_name = case_when(
    sci_name == "Aquarana catesbeianus"  ~ "Lithobates catesbeianus",
    sci_name == "Lithobates pipiens"     ~ "Rana pipiens",
    sci_name == "Aquarana clamitans"     ~ "Rana clamitans",
    sci_name == "Aquarana virgatipes"    ~ "Lithobates virgatipes",
    sci_name == "Mesotriton alpestris"   ~ "Ichthyosaura alpestris",
    sci_name == "Litoria chloris"        ~ "Chlorohyla chloris",
    sci_name == "Litoria aurea"          ~ "Ranoidea aurea",
    sci_name == "Litoria lesueurii"      ~ "Rhyaconastes lesueurii",
    TRUE ~ sci_name # Keep all other species names unchanged
  ))

print("Filtering IUCN polygons to include only native ranges and species within our dataset...")

iucn_filtered <- iucn_master %>%
  # Keep only species listed in the database
  filter(sci_name %in% species_original) %>%
  # Keep only native (origin == 1 or 2) and extant (presence == 1) polygons
  filter(presence == 1 & origin %in% c(1, 2))

# Ensure a unique polygon for each species
IUCN <- iucn_filtered %>%
  group_by(sci_name) %>%
  summarize(geometry = st_union(geometry)) %>%
  ungroup() %>%
  # Repair topology faults created during unions or splitting
  st_make_valid() 

# Rename columns
names(IUCN)[names(IUCN) == "sci_name"] <- "species"
names(IUCN)[names(IUCN) == "geometry"] <- "geom"
st_geometry(IUCN) <- "geom"

# Export the final file
st_write(IUCN, "../data/IUCN_final.gpkg", delete_dsn = TRUE, quiet = TRUE)
print("IUCN polygons loaded and correctly elaborated. Final output saved as 'IUCN_final.gpkg'.")




# ------------------------------------------------------------------------------
# Load GBIF and WorldClim DEM
# ------------------------------------------------------------------------------

# Load WorldClim elevation raster, adjust to your path (Download link: https://www.worldclim.org/data/worldclim21.html)
unzip("../data/wc2.1_30s_elev.zip", exdir="../data")
elevation <- rast("../data/wc2.1_30s_elev.tif")

# Download and import GBIF data directly using DOIs
gbif_dataset_1 <- occ_download_get("0008706-260423192947929") %>% occ_download_import()
gbif_dataset_2 <- occ_download_get("0071961-260519110011954") %>% occ_download_import()

# Combine both datasets and perform initial filtering
gbif_data <- bind_rows(gbif_dataset_1, gbif_dataset_2) %>%
  select(gbifID, species, decimalLatitude, decimalLongitude, countryCode, 
         coordinateUncertaintyInMeters, genus, order) %>%
  
  # Drop rows with missing coordinates or species
  drop_na(decimalLongitude, decimalLatitude, species) %>%
    # Filter for valid species names and acceptable uncertainty
  filter(
    species != "",
    coordinateUncertaintyInMeters <= 1000 | is.na(coordinateUncertaintyInMeters)
  ) %>%
    # Remove occurence duplicates (due to previous merge)
  distinct(gbifID, .keep_all = TRUE)

# Free up RAM by removing raw downloads
rm(gbif_dataset_1, gbif_dataset_2)
gc() 


# ------------------------------------------------------------------------------
# Match GBIF occurrences to IUCN species native range
# ------------------------------------------------------------------------------

# Align taxonomy to Frost (2026)
gbif_data <- gbif_data %>%
mutate(species = case_when(
  species == "Ranoidea chloris"        ~ "Chlorohyla chloris",
  species == "Lithobates clamitans"    ~ "Rana clamitans",
  species == "Ranoidea lesueurii"      ~ "Rhyaconastes lesueurii",
  TRUE ~ species # Keep all other names as they are
))
gbif_filtered_list <- list()
all_species_initial <- unique(gbif_data$species)

for(sp in all_species_initial) {
  
  # Subset points and polygon for the specific species
  sp_points <- gbif_data %>% filter(species == sp)
  sp_poly <- IUCN %>% filter(species == sp)
  
  # Only proceed if we have both points and a polygon for the species
  if(nrow(sp_poly) > 0 && nrow(sp_points) > 0) {
    
    # Convert GBIF dataframe to sf object using IUCN's CRS to ensure perfect match
    sp_points_sf <- st_as_sf(sp_points, 
                             coords = c("decimalLongitude", "decimalLatitude"), 
                             crs = st_crs(IUCN), 
                             remove = FALSE) # keep original lon/lat columns
    
    # Keep only points inside the IUCN polygon
    sp_points_in_range <- st_filter(sp_points_sf, sp_poly)
    
    # Drop spatial geometry to return to a standard dataframe
    gbif_filtered_list[[sp]] <- st_drop_geometry(sp_points_in_range)
  }
}

# Combine all spatially filtered species back into a single dataframe
gbif_filtered <- bind_rows(gbif_filtered_list)



# ------------------------------------------------------------------------------
# Clean GBIF coordinates (Capitals, centroids, sea...)
# ------------------------------------------------------------------------------

#Split the occurrences in batches to avoid RAM shortages,
all_species <- unique(gbif_filtered$species)
num_batches <- min(50, length(all_species)) # Prevent empty batches if species count is low
species_batches <- split(all_species, cut(seq_along(all_species), num_batches, labels = FALSE))

results_list <- list()

for(i in seq_along(species_batches)) {
  
  message(sprintf("Processing batch %d of %d...", i, num_batches))
  
  batch_data <- gbif_filtered %>% filter(species %in% species_batches[[i]])
  
  batch_clean <- tryCatch({
    clean_coordinates(
      x = batch_data,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "species", 
      countries = if("countryCode" %in% names(batch_data)) "countryCode" else NULL,
      tests = c("capitals", "centroids", "equal", "gbif", "institutions","sea", "zeros", "outliers"),
      value = "clean", 
      verbose = FALSE 
    )
  }, error = function(e) {
    message("Error in batch ", i, ": ", e$message)
    return(NULL)
  })
  
  if(!is.null(batch_clean)) {
    results_list[[i]] <- batch_clean
  }
  
  gc() # Force R to flush memory between batches
}

gbif_final_clean <- bind_rows(results_list)


# ------------------------------------------------------------------------------
# Intersection between GBIF and WorldClim elevation raster (Species elevation range extraction)
# ------------------------------------------------------------------------------

# Convert cleaned dataframe into a SpatVector for terra extraction
occ_vect <- vect(gbif_final_clean, 
                 geom = c("decimalLongitude", "decimalLatitude"), 
                 crs = "EPSG:4326")

# Extract values from the WorldClim elevation raster
elev_extracted <- terra::extract(elevation, occ_vect)

# Append extracted values back to the cleaned dataframe
gbif_final_clean$elev_value <- elev_extracted[, 2]

# Group by species and calculate elevation metrics
elevation_summary <- gbif_final_clean %>%
  group_by(species) %>%
  summarize(
    # Select 5th and 95th quantile, as an additional control for outliers
    perc_05_elevation = quantile(elev_value, probs = 0.05, na.rm = TRUE),
    perc_95_elevation = quantile(elev_value, probs = 0.95, na.rm = TRUE),
    num_occurrences = n() # Total verified points for each species
  )


# ------------------------------------------------------------------------------
# Extract species range coordinates from IUCN polygons
# ------------------------------------------------------------------------------

# Merge polygons by species
iucn_metrics <- IUCN %>%
  group_by(species) %>%
  summarize(geom = st_union(geom))

# Reproject to Mollweide equal-area to prevent poleward bias
iucn_ea <- iucn_metrics %>%
  mutate(geom_ea = st_transform(geom, crs = "+proj=moll +datum=WGS84"))

# Compute centroid in equal-area space, then back-project to WGS84 for coordinates
centroids_ea    <- st_centroid(iucn_ea$geom_ea)
centroids_wgs84 <- st_transform(centroids_ea, crs = 4326)
centroid_coords <- st_coordinates(centroids_wgs84)

iucn_metrics$centroid_lat  <- centroid_coords[, "Y"]

# Calculate absolute latitude
iucn_metrics$abs_lat <- abs(iucn_metrics$centroid_lat)

# Drop heavy geometry before merging
iucn_metrics_clean <- iucn_metrics %>% st_drop_geometry()



# ------------------------------------------------------------------------------
# Merge and Export Final Data
# ------------------------------------------------------------------------------

# Join the IUCN metrics with the elevation summary
final_spatial_metrics <- iucn_metrics_clean %>%
  left_join(elevation_summary, by = "species")

# Export to a single unified CSV
write.csv(final_spatial_metrics, "../data/Species_Spatial.csv", row.names = FALSE)

cat("\nProcess completed, your spatial data are available to data/Species_Spatial.csv\n")
