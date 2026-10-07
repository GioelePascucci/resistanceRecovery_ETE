# ==============================================================================
# Script title: heatTolerance.R
# Description:  The script standardizes CTmax values from different studies to the
#               same acclimation temperature, acclimation time and endpoint.
# Manuscript title: Correlated evolution of thermal tolerance and egg size 
#                links resistance and recovery to extreme heat in amphibians.
# Manuscript authors: Gioele Pascucci, Miguel Tejedo, Lumír Gvoždík
# Script author: Gioele Pascucci
# Contact:  gioelepascucci@live.it
# Date: August 30, 2026
# R Version:  4.5.3
# ==============================================================================

rm(list = ls()) # Clean workspace

library(lme4)
library(emmeans)
library(dplyr)
library(lmerTest)
library(pbkrtest)
library(readxl)
library(tidyr)
library(writexl)

# ------------------------------------------------------------------------------
# Load and prepare datasets: AmphiTherm v1.01 (Pottier et al., 2025) & datasetTrait.xlsl
# ------------------------------------------------------------------------------

# Load AmphiTherm dataset, adjust to your path (Download link: https://doi.org/10.5281/zenodo.16751427)
# Unzip the file into /data directory
unzip("../data/AmphiTherm-v1.0.1.zip", exdir = "../data")
amphiTherm<-read.csv("../data/p-pottier-AmphiTherm-234a12d/data/Cleaned_and_curated_data_Amphitherm.csv")

amphiTherm<-amphiTherm %>% 
  select(species, population_ID, cohort_ID, acclimation_temp, acclimation_time, metric, endpoint, mean_trait, n_trait, strain, peer_reviewed, acclimated,life_stage_acclimated,life_stage_tested, chemical, hormone, infected) %>%
  # Keep only observation for CTmax and adult specimens, from peer_reviewed sources, without chemicals and/or infections influences
  filter(metric == "CTmax") %>%
  filter(peer_reviewed == "peer-reviewed") %>% 
  filter(acclimated == "acclimated") %>% 
  filter(life_stage_acclimated == "adults") %>% 
  filter(life_stage_tested == "adults") %>% 
  filter(
    (chemical %in% c("control", "NA") | is.na(chemical)) &
    (hormone %in% c("control", "NA") | is.na(hormone)) &
    (infected == "NA" | is.na(infected))
  ) %>% 
  # Align taxonomy to Frost (2026)
  mutate(
    species = case_when(
      species == "Desmognathus fuscus" & strain == "auriculatus" ~ "Desmognathus auriculatus",
      species == "Desmognathus fuscus" & strain == "brimleyorum" ~ "Desmognathus brimleyorum",
      species == "Desmognathus ochrophaeus" & strain == "carolinensis" ~ "Desmognathus carolinensis",
      species == "Plethodon cinereus" & strain == "cinereus" ~ "Plethodon cinereus",
      species == "Plethodon wehrlei" & strain == "dixi" ~ "Plethodon dixi",
      species == "Anaxyrus woodhousii" & strain == "fowleri" ~ "Anaxyrus fowleri",
      species == "Eurycea longicauda" & strain == "guttolineata" ~ "Eurycea guttolineata",
      species == "Plethodon jordani" & strain == "jordani" ~ "Plethodon jordani",
      species == "Eurycea lucifuga" & strain == "lucifuga" ~ "Eurycea lucifuga",
      species == "Eurycea longicauda" & strain == "melanopleura" ~ "Eurycea melanopleura",
      species == "Plethodon jordani" & strain == "metcalfi" ~ "Plethodon metcalfi",
      species == "Litoria chloris" ~ "Chlorohyla chloris",
      species == "Litoria aurea" ~ "Ranoidea aurea",
      species == "Litoria lesueurii" ~ "Rhyaconastes lesueurii",
      species == "Hyla andersonii" ~ "Dryophytes andersonii",
      species == "Hyla chrysoscelis" ~ "Dryophytes chrysoscelis",
      species == "Hyla cinerea" ~ "Dryophytes cinereus",
      species == "Hyla versicolor" ~ "Dryophytes versicolor",
      species == "Rana catesbeiana" ~ "Lithobates catesbeianus",
      species == "Rana sylvatica" ~ "Lithobates sylvaticus",
      species == "Rana virgatipes" ~ "Lithobates virgatipes",
      species == "Rhinella schneideri" ~ "Rhinella diptycha",
      
      # keep all the other species name
      TRUE ~ species 
    )
  ) %>% 
  select(species, population_ID, cohort_ID, acclimation_temp, acclimation_time, endpoint, mean_trait, n_trait) %>% 
  # Check for missing values
  mutate(across(where(is.character), ~ na_if(.x, "NA"))) %>%
  drop_na()

# Load trait dataset, adjust to your path
main_data <- read_excel("../data/datasetTraits.xlsx")

# Keep only species present in the main dataset
amphiCTmax <- amphiTherm %>%
  filter(species %in% main_data$Species)
# Set column type
amphiCTmax<-amphiCTmax %>% 
  mutate(
    species = as.factor(species),
    population_ID = as.factor(population_ID),
    cohort_ID = as.factor(cohort_ID),
    acclimation_temp=as.numeric(acclimation_temp),
    acclimation_time=as.numeric(acclimation_time),
    endpoint=as.factor(endpoint),
    mean_trait=as.numeric(mean_trait),
    n_trait=as.numeric(n_trait),
  )



# ------------------------------------------------------------------------------
# Check which species are suitable for standardization (at least 2 acclimation temperature)
# ------------------------------------------------------------------------------
# Count distinct acclimation temperatures per species
ctmax_counts <- amphiCTmax %>%
  group_by(species) %>%
  summarise(n_temps = n_distinct(acclimation_temp, na.rm = TRUE))

# Identify valid species and singletons (species tested only at one acclimation temperature)
valid_species <- ctmax_counts %>% filter(n_temps >= 2) %>% pull(species)
singleton_species <- ctmax_counts %>% filter(n_temps < 2) %>% pull(species)

# Split the dataset
valid_data_ctmax <- amphiCTmax %>% filter(species %in% valid_species)
singleton_data_ctmax <- amphiCTmax %>% filter(species %in% singleton_species)



# ------------------------------------------------------------------------------
# Fit Linear Mixed-Effects Model & extract standardized CTmax
# ------------------------------------------------------------------------------

# Modeling CTmax using acclimation temperature and species interaction,
# accounting for experimental endpoint and acclimation time,
# with random effects for nested data within population.
# weighting observation by sample size (n_trait)
model_ctmax <- lmer(mean_trait ~ acclimation_temp * species + endpoint + acclimation_time + (1 | population_ID/cohort_ID), 
                    weights = n_trait, 
                    data = valid_data_ctmax)

# Calculate Estimated Marginal Means (emmeans) at standard conditions:
# acclimation temp = 20, acclimation time = 7 days, endpoint = "LRR"
emms_ctmax <- emmeans(model_ctmax, ~ species, 
                      at = list(acclimation_temp = 20, 
                                acclimation_time = 7, 
                                endpoint = "LRR"))

# Convert to a clean dataframe
results_lmm <- as.data.frame(emms_ctmax) %>%
  rename(CTmax = emmean) %>%
  mutate(
    CTmax = round(CTmax, 2),
    acclimation_temp = 20,
    acclimation_time = 7,
    endpoint = "LRR",
    data_quality = "LMM_Standardized"
  ) %>%
  select(species, acclimation_temp, acclimation_time, endpoint, CTmax, data_quality)



# ------------------------------------------------------------------------------
# Rescue singletons acclimated at 20°C, but at different acclimation time and endpoint
# ------------------------------------------------------------------------------

# Rescue species excluded from the main model,
# Since they were already acclimated at exactly 20°C,
# we just need to standardize acclimation time and endpoint,
# don't need a minimum number of observations to calculate "acclimation_temp * species"

# Extract coefficients from previous model
coefs <- as.data.frame(fixef(model_ctmax))
colnames(coefs) <- c("Effect")
coefs$Term <- rownames(coefs)

# Set the same standardization target
target_endpoint <- "LRR"
target_time <- 7
target_endpoint_term <- paste0("endpoint", target_endpoint)
target_effect_lss <- coalesce(coefs$Effect[match(target_endpoint_term, coefs$Term)], 0)
time_slope <- coalesce(coefs$Effect[match("acclimation_time", coefs$Term)], 0)

rescued_singletons <- singleton_data_ctmax %>%
  filter(acclimation_temp == 20) %>%
  mutate(
    # Pull the mathematical effect of the observed endpoint
    obs_endpoint_term = paste0("endpoint", endpoint),
    obs_effect = coalesce(coefs$Effect[match(obs_endpoint_term, coefs$Term)], 0),
    
    # Calculate the standardized CTmax
    Adj_CTmax = mean_trait - obs_effect + target_effect_lss + ((target_time - acclimation_time) * time_slope)
  ) %>%
  group_by(species) %>%
  summarise(
    CTmax = round(weighted.mean(Adj_CTmax, w = n_trait, na.rm = TRUE), 2),
    acclimation_temp = 20,
    acclimation_time = 7,
    endpoint = "LRR",
    .groups = "drop"
  ) %>%
  mutate(data_quality = "Rescued_Singleton_20C") %>%
  select(species, acclimation_temp, acclimation_time, endpoint, CTmax, data_quality)



# ------------------------------------------------------------------------------
# Process remaining singletons
# ------------------------------------------------------------------------------

remaining_singletons <- singleton_data_ctmax %>%
  filter(!species %in% rescued_singletons$species) %>%
  group_by(species) %>%
  # If present keep the observations with LRR as endpoint
  filter(if (any(endpoint == "LRR", na.rm = TRUE)) endpoint == "LRR" else TRUE) %>%
  # Calculate weighted mean if multiple observations for same experimental conditions are available
  summarise(
    CTmax = round(weighted.mean(mean_trait, w = n_trait, na.rm = TRUE), 2),
    acclimation_temp = mean(acclimation_temp, na.rm = TRUE),
    acclimation_time = mean(acclimation_time, na.rm = TRUE),
    endpoint = first(endpoint), 
    .groups = "drop"
  ) %>%
  mutate(data_quality = "No_standardization") %>%
  select(species, acclimation_temp, acclimation_time, endpoint, CTmax, data_quality)

# Merge the outputs
final_ctmax <- bind_rows(results_lmm, rescued_singletons, remaining_singletons)

# Print the summary of standardization process
cat("\n--- STANDARDIZATION SUMMARY (Acclimation time=7 days; acclimation temperature=20 °C; endpoint=LRR ---\n")
cat("Full Model Predictions (Species with multiple acclimation temperatures available): ", nrow(results_lmm), "\n")
cat("Species already (and only) acclimated at 20 °C: ", nrow(rescued_singletons), "\n")
cat("Not-standardized species:     ", nrow(remaining_singletons), "\n")
cat("Total Species:", nrow(final_ctmax), "\n\n")

# Export to a single Excel file
write_xlsx(final_ctmax, "../data/CTmax_Standardized.xlsx")
