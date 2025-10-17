#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

library(auk)
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(sf)
library(magrittr)
library(tidyr)
library(purrr)
library(units)


#########################
### SPECIES SELECTION ###
#########################

# Tally counts of all species across, for each Atlas period
species_counts_summary <- observations_comp %>% 
  count(common_name, atlas_period, name = "n_obs") %>% 
  group_by(common_name) %>% 
  summarise(
    obs_Atlas1 = sum(n_obs[atlas_period == "Atlas1"], na.rm = TRUE),
    obs_Atlas2 =sum(n_obs[atlas_period == "Atlas2"], na.rm = TRUE),
    obs_Total = sum(n_obs,na.rm = TRUE),
    .groups = "drop"
  )

# Tally counts of breeding species across, for each Atlas period
breeders_counts_summary <- observations_comp %>% 
  filter(
    breeding_category %in% c("C2", "C3", "C4")
  ) %>% 
  count(common_name, atlas_period, name = "n_obs") %>% 
  group_by(common_name) %>% 
  summarise(
    obs_Atlas1 = sum(n_obs[atlas_period == "Atlas1"], na.rm = TRUE),
    obs_Atlas2 =sum(n_obs[atlas_period == "Atlas2"], na.rm = TRUE),
    obs_Total = sum(n_obs,na.rm = TRUE),
    .groups = "drop"
  )


# Develop minimum # of observations to be usable in modeling - WIP





####################
### ZERO-FILLING ###
####################

# Create detection/non-detection data


### --- ALL SPECIES PROTOCOL --- ###

# Extract all unique species names
species_list <- unique(observations_comp$common_name)

# Function to automate zero-fill
Zerofill_species <- function(sp_name, checklists) {
  
  obs_sp <- observations_comp %>%
    filter(
      common_name == sp_name,
      breeding_category %in% c("C2", "C3", "C4")
    )
  
  zf_sp <- auk_zerofill(obs_sp, checklists, collapse = TRUE)
  
  dets_summary <- zf_sp %>%
    group_by(atlas_block, atlas_period) %>%
    summarise(species_observed = any(species_observed), .groups = "drop") %>%
    pivot_wider(
      names_from = atlas_period,
      values_from = species_observed,
      names_prefix = "det_",
      values_fill = FALSE
    ) %>%
    # ensure both columns exist, even if species missing in one atlas
    mutate(
      det_Atlas1 = if (!"det_Atlas1" %in% names(.)) FALSE else det_Atlas1,
      det_Atlas2 = if (!"det_Atlas2" %in% names(.)) FALSE else det_Atlas2
    ) %>%
    mutate(
      det_Atlas1 = as.integer(det_Atlas1),
      det_Atlas2 = as.integer(det_Atlas2),
      transition_state = case_when(
        det_Atlas1 == 0 & det_Atlas2 == 0 ~ "Absence",
        det_Atlas1 == 0 & det_Atlas2 == 1 ~ "Colonization",
        det_Atlas1 == 1 & det_Atlas2 == 0 ~ "Extinction",
        det_Atlas1 == 1 & det_Atlas2 == 1 ~ "Persistence",
        TRUE ~ NA_character_
      ),
      common_name = sp_name
    )
  
  return(dets_summary)
}

# Apply function
breeders_zf_summary <- map_dfr(species_list, Zerofill_species, checklists = checklists_comp)



### --- SINGLE SPECIES PROTOCOL --- ###

# Zero-fill all comp checklists with species det/nondet (TRUE/FALSE)
obs_rcki <- observations_comp %>% # filter obs to single species
  filter(
    common_name == "Ruby-crowned Kinglet",
    breeding_category %in% c("C2", "C3", "C4") # only breeding individuals
  ) %>%
  mutate(all_species_reported = TRUE) # force incomplete checklists to complete (all Atlas1 data coded "incomplete")

lists_rcki <- checklists_comp %>% # include all comp checklists
  mutate(all_species_reported = TRUE)

zf_rcki <- auk_zerofill(obs_rcki, lists_rcki, collapse = TRUE) # zf df

# Summarize detections by block, Atlas period
rcki_dets_summary <- zf_rcki %>%
  group_by(atlas_block, atlas_period) %>%
  summarise(species_observed = any(species_observed), .groups = "drop") %>%
  pivot_wider(
    names_from = atlas_period,
    values_from = species_observed,
    names_prefix = "det_",
    values_fill = FALSE
  ) %>%
  mutate(
    det_Atlas1 = as.integer(det_Atlas1), #convert T/F to 1/0
    det_Atlas2 = as.integer(det_Atlas2)
  )

# Define transition states
rcki_dets_summary <- rcki_dets_summary %>%
  mutate(
    transition_state = case_when(
      det_Atlas1 == 0 & det_Atlas2 == 1 ~ "Colonization",
      det_Atlas1 == 1 & det_Atlas2 == 0 ~ "Extinction",
      det_Atlas1 == 1 & det_Atlas2 == 1 ~ "Persistence",
      det_Atlas1 == 0 & det_Atlas2 == 0 ~ "Absence",
      TRUE ~ NA_character_
    )
  )



#################
### VISUALIZE ###
#################

# Visualize state transitions between Atlas periods (single species)
blocks_comp_rcki <- blocks_comp_shp %>% left_join(rcki_dets_summary, by = "atlas_block") # join data, shp 

ggplot() + # map dets across comp blocks
  geom_sf(data = blocks_comp_rcki, aes(fill = transition_state), size = 0.1) +
  scale_fill_manual(
    values = c(
      "Colonization" = "#FFD700",   # gold/yellow
      "Extinction"   = "#1E90FF",   # blue
      "Persistence"  = "#32CD32",   # green
      "Absence"      = "gray80"     # light gray
    ),
    na.value = "white"
  ) +
  labs(
    title = "Ruby-crowned Kinglet Atlas Block Transitions",
    fill = "Transition State"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )