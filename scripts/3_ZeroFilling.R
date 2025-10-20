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



#################
### VISUALIZE ###
#################

### --- SPECIES DATA --- ###

# Reproducible for any species 
species_to_plot <- "Bobolink"

# Join species data to blocks shp data by block ID
blocks_species <- blocks_comp_shp %>%
  left_join(
    breeders_zf_summary %>%
      filter(common_name == species_to_plot),
    by = "atlas_block"
  )


### --- PLOT --- ###

vir_colors <- viridis::viridis(3) # make color blind-friendly palette

ggplot(blocks_species) +
  geom_sf(aes(fill = transition_state), size = 0.1) +
  scale_fill_manual( # hacking manual fill with viridis palette
    values = c(
      "Colonization" = vir_colors[3], # assign accessible colors to states of interest
      "Persistence"  = vir_colors[2],
      "Extinction"   = vir_colors[1],
      "Absence"      = "white" # force absence to white so it doesn't distract
    ),
    breaks = c("Colonization", "Persistence", "Extinction"), # show only states of relevance
    drop = FALSE
  ) +
  labs(
    fill = "Transition State",
    caption = "Figure 1. Map of state transitions for the Bobolink (Dolichonyx oryzivorus) across comparable survey blocks using data from the Wisconsin Breeding Bird Atlas. \nTransitions reflect detection/non-detection data aggregated at the block level (~25 km^2) across two survey periods: 1995-2000 (Atlas 1) and 2015-2019 (Atlas 2). \nWhite blocks denote those grid cells part of the comparable pool but in which the species was not detected in either Atlas period."
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13, margin = margin(b = 10)),
    legend.text = element_text(size = 11),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0.5, size = 9, margin = margin(t = 20))
  )