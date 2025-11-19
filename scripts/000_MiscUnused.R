########################
### MISC CODE CHUNKS ###
########################

### --- ELEVATION DATA (NED) --- ###
# Mean elevation per block done in ArcPro (1/3 arc-sec DEMs,~10m)

elevation_summary <- read_csv("maps/ned/ned_wibba_table.csv") # read in data calc'd in ArcPro

elevation_summary <- elevation_summary %>%
  rename(
    atlas_block = BLOCK_ID, #rename variables to match
    elev_mean = MEAN
  ) %>%
  mutate(atlas_block = as.character(atlas_block)) %>%
  mutate(
    elev_mean_z = as.numeric(scale(elev_mean)) # z-standard
  )

wibba_summary <- wibba_summary %>%
  left_join(elevation_summary %>% dplyr::select(atlas_block, elev_mean_z), by = "atlas_block") # join to summary df






### --- GEOGRAPHIC FILTER --- ###

### Locates checklists within Atlas block grid by lat, lon coordinates and 
# assigns missing block IDs to checklists, observations. 
# Uneeded for WIBBA.


# geographic filtering, clean-up
wibba_grid <- read_sf("maps/wibba blocks/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # read in atlas grid shape file
  st_transform(4269) # convert to NAD83, crs of shp file

checklists_sf <- checklists %>%
  select(checklist_id, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% # convert checklist locations to points
  st_transform(4269) # match atlas block crs

checklists_joined <- st_join(checklists_sf, wibba_grid[, c("BLOCK_ID")], left = TRUE) # create df with checklist id, lat/lon, and BLOCK_ID

in_region <- checklists_joined %>%
  filter(!is.na(BLOCK_ID)) # keep only checklists in joined atlas block shp

checklists_wip <- checklists %>%
  semi_join(st_drop_geometry(in_region), by = "checklist_id") %>% # remove checklists outside atlas blocks
  left_join(st_drop_geometry(in_region[, c("checklist_id", "BLOCK_ID")]), by = "checklist_id") %>% # adds BLOCK_ID from shp to checklist data w/o overwriting
  mutate(
    atlas_block = if_else(is.na(atlas_block), BLOCK_ID, atlas_block) # fills in missing atlas blocks w/o overwriting
  ) %>%
  select(-BLOCK_ID) #remove BLOCK_ID column

observations_wip <- observations %>%
  semi_join(checklists, by = "checklist_id") %>%  # include only matching checklist_id from checklists df
  left_join(
    checklists %>% 
      select(checklist_id, atlas_block) %>%  # join atlas_block to observations using checklist_id
      rename(joined_block = atlas_block),
    by = "checklist_id"
  ) %>%
  mutate(
    old_block = atlas_block, # preserve existing
    atlas_block = joined_block # overwrite with cleaned value
  ) %>%
  select(-joined_block) # drop temp join column



### --- OLD MODELING APPROACHES --- ###



######################### OLD MODEL APPROACHES #################################


# GROUPED DATA # 
# 'Colonization-Persistence' df
data_colper <- spp_mod_data_full %>%
  filter(transition_state %in% c("Colonization", "Persistence"))

data_colper_z <- data_colper %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%    
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))



# 'Extinction-Absence' df
data_extabs <- spp_mod_data_full %>%
  filter(transition_state %in% c("Extinction", "Absence"))

data_extabs_z <- data_extabs %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))


# Binomial state IDs
data_colper_z <- data_colper_z %>%
  mutate(
    col_per = ifelse(transition_state %in% c("Colonization"), 1, 0),
  )

data_extabs_z <- data_extabs_z %>%
  mutate(
    ext_abs = ifelse(transition_state %in% c("Extinction"), 1, 0)
  )


# Col-Per
rbwo_mod_colper <- glm(col_per ~  pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_colper_z, family = binomial)
summary(rbwo_mod_colper)
autoplot(rbwo_mod_colper)

# Ext-Abs
rbwo_mod_extabs <- glm(ext_abs ~ pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_extabs_z, family = binomial)

summary(rbwo_mod_extabs)
autoplot(rbwo_mod_extabs)


# Visualize

# COL-PER
facet_colper <- c(
  pa_percent_z = "Protected Area",
  developed_lower_base_z = "Open + Lower Development",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest ***",
  forest_deciduous_base_z = "Deciduous Forest .",
  pasture_crop_base_z = "Pasture/Cropland",
  tmax_38yr_z = "Max Temp ***",
  prcp_38yr_z = "Precipitation",
  developed_total_diff_z = "Difference in Total Developed Land",
  forest_total_diff_z = "Difference in Total Forest"
)

rbwo_predplot_colper_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_colper, 
                                   data = data_colper_z)


ggplot(rbwo_predplot_colper_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_colper)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across Blocks"
  )


# EXT-ABS
facet_extabs <- c(
  pa_percent_z = "Protected Area",
  developed_lower_base_z = "Open + Lower Development",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest **",
  forest_deciduous_base_z = "Deciduous Forest",
  pasture_crop_base_z = "Pasture/Cropland",
  tmax_38yr_z = "Max Temp ***",
  prcp_38yr_z = "Precipitation",
  developed_total_diff_z = "Difference in Total Developed Land",
  forest_total_diff_z = "Difference in Total Forest"
)

rbwo_predplot_ext_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                model = rbwo_mod_ext, 
                                data = spp_mod_data_ext_z)


ggplot(rbwo_predplot_ext_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_extabs)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across Blocks"
  )
