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



####################
### ZERO-FILLING ###
####################

## --- CREATE DETECTION/NON-DETECTION DATA --- ##

# Single species protocol

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


## --- VISUALIZE DETECTION-BASED DATA --- ##

# Visualize species transitions between Atlas periods
blocks_comp_rcki <- blocks_comp_shp %>% left_join(rcki_dets_summary, by = "atlas_block") # join data, shp 

ggplot() + # map dets across comp blocks
  geom_sf(data = blocks_comp_rcki, aes(fill = transition_state), color = "black", size = 0.1) +
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