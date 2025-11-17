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
library(here)



########################
### PREP EBIRD FILES ###
########################

## --- HELPER FUNCTIONS --- ##

# Assign atlas period based on observation year
atlas_periods <- function(date) {
  case_when(
    year(date) %in% 1995:2000 ~ "Atlas1",
    year(date) %in% 2015:2019 ~ "Atlas2",
    TRUE ~ NA_character_
  )
}

# Project filters for atlas data
filter_atlas_data <- function(df) {
  df %>%
    filter(
      year(observation_date) %in% c(1995:2000, 2015:2019),
      project_names == "Wisconsin Breeding Bird Atlas"
    )
}

## --- UPLOAD & FILTER EBD, SED FILES --- ## 

# Set auk path
auk_set_ebd_path("C:/Users/lablue/Desktop/MS Project/WI Atlas/WIBBA Data/data/ebird/ebird_raw", overwrite = TRUE)

# Set relative file paths
wibba_ebd_base_raw <- "data/ebird/ebird_raw/ebd_US-WI_smp_relMar-2025.txt"
wibba_sed_raw <- "data/ebird/ebird_raw/ebd_US-WI_smp_relMar-2025_sampling.txt"
wibba_ebd_sensitive_raw <- "data/ebird/ebird_raw/ebd_sensitive_wibba_relSep2025.txt"

# Auk filtering full ebd, sed
auk_ebd(wibba_ebd_base_raw, file_sampling = wibba_sed_raw) %>% # joint filtering ebd, sed files; all wibba1 data are coded "incomplete"
  auk_date(date = c("*-05-01", "*-08-31")) %>% # May-August, any year
  auk_filter(file = "data/ebird/ebird_filtered/wibba_ebd_base_filtered.txt", file_sampling = "data/ebird/ebird_filtered/wibba_sed_filtered.txt", overwrite = TRUE) # rerun using , overwrite = TRUE

# Auk filtering sensitive ebd
auk_ebd(wibba_ebd_sensitive_raw, file_sampling = wibba_sed_raw) %>% # joint filtering sensitive ebd, sed files
  auk_date(date = c("*-05-01", "*-08-31")) %>% # June-August, any year
  auk_filter(file = "data/ebird/ebird_filtered/wibba_ebd_sensitive_filtered.txt", file_sampling = "data/ebird/ebird_filtered/wibba_sed_filtered.txt", overwrite = TRUE) # rerun using , overwrite = TRUE

# Additional obs filtering
observations_base <- read_ebd("data/ebird/ebird_filtered/wibba_ebd_base_filtered.txt") %>% filter_atlas_data()

observations_sensitive <- read_ebd("data/ebird/ebird_filtered/wibba_ebd_sensitive_filtered.txt") %>% filter_atlas_data()

observations_filt <- bind_rows(observations_base, observations_sensitive)

# Additional checklist filtering
checklists_filt <- read_sampling("data/ebird/ebird_filtered/wibba_sed_filtered.txt") %>%
  filter_atlas_data()

# Drop unmatched obs
observations_filt <- semi_join(observations_filt, checklists_filt, by = "checklist_id")

# Add atlas period to obs, lists dfs
checklists_filt <- checklists_filt %>% mutate(atlas_period = atlas_periods(observation_date))
observations_filt <- observations_filt %>% mutate(atlas_period = atlas_periods(observation_date))


##############################
### SUMMARIZING ATLAS DATA ###
##############################

## --- Summarize Checklists --- ##

checklist_summary <- checklists_filt %>%
  count(atlas_block, atlas_period, name = "n_checklists") %>%
  pivot_wider(
    names_from = atlas_period,
    values_from = n_checklists,
    names_prefix = "checklists_"
  ) %>%
  mutate(across(starts_with("checklists_"), ~replace_na(., 0))) %>%
  mutate(checklists_total = checklists_Atlas1 + checklists_Atlas2)

  
## --- Summarize Temporal Effort (mins) --- ##

minutes_summary <- checklists_filt %>%
  group_by(atlas_block, atlas_period) %>%
  summarise(effort_minutes = sum(duration_minutes, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = atlas_period,
    values_from = effort_minutes,
    names_prefix = "mins_"
  ) %>%
  mutate(across(starts_with("mins_"), ~replace_na(., 0))) %>%
  mutate(mins_total = mins_Atlas1 + mins_Atlas2)

supp_effort <- read_csv("data/effort/wibba1_effort.csv") %>%
  filter(!is.na(quad_prefix), !is.na(quad_suffix)) %>%
  mutate(
    quad_prefix = as.character(quad_prefix),
    atlas_block = paste0(quad_prefix, quad_suffix),
    effort_min = effort_hr * 60
  ) %>%
  group_by(atlas_block) %>%
  summarise(mins_Atlas1_supp = sum(effort_min, na.rm = TRUE), .groups = "drop")


## --- Summarize Species Richness --- ##

sr_summary <- observations_filt %>%
  group_by(atlas_block, atlas_period) %>%
  summarise(species = list(unique(common_name)), .groups = "drop") %>%
  complete(atlas_block, atlas_period = c("Atlas1", "Atlas2")) %>%
  pivot_wider(names_from = atlas_period, values_from = species) %>%
  mutate(
    sr_Atlas1 = lengths(Atlas1),
    sr_Atlas2 = lengths(Atlas2),
    sr_total = lengths(map2(Atlas1, Atlas2, union))
  ) %>%
  left_join(select(checklist_summary, atlas_block, checklists_Atlas1, checklists_Atlas2), by = "atlas_block") %>%
  mutate(
    sr_Atlas1 = ifelse(checklists_Atlas1 == 0, NA, replace_na(sr_Atlas1, 0)),
    sr_Atlas2 = ifelse(checklists_Atlas2 == 0, NA, replace_na(sr_Atlas2, 0)),
    sr_total = ifelse(is.na(sr_Atlas1) & is.na(sr_Atlas2), NA, replace_na(sr_total, 0)),
    sr_diff = sr_Atlas2 - sr_Atlas1,
    sr_diff_z = scale(sr_diff)[, 1]
  ) %>%
  select(atlas_block, sr_total, sr_Atlas1, sr_Atlas2, sr_diff, sr_diff_z)


## --- ATLAS BLOCK SUMMARY DF --- ##

blocks_all <- checklists_filt %>%
  distinct(atlas_block) %>%
  arrange(atlas_block)

wibba_summary_full <- blocks_all %>%
  left_join(checklist_summary, by = "atlas_block") %>%
  left_join(minutes_summary, by = "atlas_block") %>%
  left_join(supp_effort, by = "atlas_block") %>%
  mutate(
    mins_Atlas1 = coalesce(mins_Atlas1_supp, mins_Atlas1),
    mins_total = rowSums(across(c(mins_Atlas1, mins_Atlas2)), na.rm = TRUE),
    mins_diff = mins_Atlas2 - mins_Atlas1,
    mins_diff_z = scale(mins_diff)[, 1]
  ) %>%
  select(-mins_Atlas1_supp)

wibba_summary_full <- wibba_summary %>%
  left_join(sr_summary, by = "atlas_block")



#######################
### COMPARABLE DATA ###
#######################

## --- FILTER COMPARABLE BLOCKS (RLL) --- ##

# Remove blocks missing checklists in either atlas period
# AND blocks with 0 sr, which are not "empty" sed lists but lists that had invalidated bird records that are not in the ebd

wibba_summary_comp <- wibba_summary_full %>%
  mutate(across(
    c(checklists_Atlas1, checklists_Atlas2, sr_Atlas1, sr_Atlas2), 
    ~replace_na(., 0)
  )) %>%
  filter(
    checklists_Atlas1 > 0,
    checklists_Atlas2 > 0,
    sr_Atlas1 > 0,
    sr_Atlas2 > 0
  )

sr_summary_comp <- sr_summary %>%
  filter(atlas_block %in% blocks_comp)

blocks_comp <- wibba_summary_comp$atlas_block

checklists_comp <- checklists_filt %>%
  filter(atlas_block %in% blocks_comp)

observations_comp <- observations_filt %>%
  filter(atlas_block %in% blocks_comp)


## -- MAPS/SHP FILES -- ##

# Load, create maps
blocks_shp <- st_read("data/maps/wibba blocks/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # load full block map
  rename(atlas_block = BLOCK_ID)

blocks_comp_shp <- blocks_shp %>% # create sf object of comparable blocks
  filter(atlas_block %in% blocks_comp)

# Visualize all Atlas blocks
ggplot(blocks_shp) + # all blocks
  geom_sf(fill = "orange", color = "white") +
  ggtitle("All Atlas Blocks")

# Visualize surveyed v un-surveyed blocks
blocks_surveyed <- checklist_summary %>% 
  select(atlas_block, checklists_Atlas1, checklists_Atlas2) %>%
  pivot_longer(
    cols = starts_with("checklists_"),
    names_to = "atlas_period",
    names_prefix = "checklists_",
    values_to = "n_checklists"
  ) %>%
  mutate(surveyed = n_checklists > 0)

blocks_surveyed_shp <- blocks_shp %>%
  left_join(blocks_surveyed, by = "atlas_block")

ggplot(blocks_surveyed_shp) +
  geom_sf(aes(fill = surveyed), color = NA) +
  facet_wrap(~atlas_period) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgray")) +
  labs(
    title = "Surveyed Atlas Blocks by Atlas Period",
    fill = "Surveyed"
  ) +
  theme_minimal()

# Visualize comparable Atlas blocks (surveyed in both Atlas periods)
ggplot(blocks_comp_shp) + 
  geom_sf(fill = "orange", color = "white") +
  ggtitle("Comparable Atlas Blocks")





### --- COMPARABLE BLOCKS (DNR) --- ###

wibba_summary_dnrcomp <- dnr_compblocks %>%
  left_join(wibba_summary_comp, by = "atlas_block")

