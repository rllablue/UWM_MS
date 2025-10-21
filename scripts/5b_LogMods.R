#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

library(dplyr)
library(ggplot2)
library(lubridate)
library(magrittr)
library(tidyr)
library(purrr)
library(units)



#########################
### 3 LOGISTIC MODELS ###
#########################

Model_logistic <- function(species_name, 
                           zf_summary = breeders_zf_summary,
                           covariates = wibba_modeling_comp) {
  
  species_dets <- zf_summary %>% # filter full species zf data
    filter(common_name == species_name)
  
  species_mod_df <- species_dets %>% # combine zf, covariate data ### WIP WIP WIP ###
    left_join(
      covariates %>%
        dplyr::select(
          atlas_block,
          pa_z,
          lat_z,
          lon_z,
          srA1_resid,
          srA2_resid,
          developed_total_z,
          forest_total_z,
          grass_pasture_crop_z,
          wetlands_total_z
        ),
      by = "atlas_block"
    )
  
  species_mod_df <- species_mod_df %>% # create state identifiers
    mutate(
      is_colonization = ifelse(transition_state == "Colonization", 1, 0),
      is_extinction   = ifelse(transition_state == "Extinction", 1, 0),
      is_persistence  = ifelse(transition_state == "Presence", 1, 0)
    )
  
  # --- Fit logistic models --- # WIP WIP WIP ###
  mod_colonization <- glm(is_colonization ~ pa_z + lat_z + lon_z + srA1_resid + srA2_resid +
                            developed_total_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z,
                          data = species_mod_df, family = binomial)
  
  mod_extinction <- glm(is_extinction ~ pa_z + lat_z + lon_z + srA1_resid + srA2_resid +
                          developed_total_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z,
                        data = species_mod_df, family = binomial)
  
  mod_persistence <- glm(is_persistence ~ pa_z + lat_z + lon_z + srA1_resid + srA2_resid +
                           developed_total_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z,
                         data = species_mod_df, family = binomial)
  
  # --- Predict probabilities --- #
  species_mod_df <- species_mod_df %>%
    mutate(
      p_colonization = predict(mod_colonization, type = "response"),
      p_extinction   = predict(mod_extinction, type = "response"),
      p_persistence  = predict(mod_persistence, type = "response")
    )
  
  # --- Output --- #
  return(list(
    models = list(
      colonization = mod_colonization,
      extinction = mod_extinction,
      persistence = mod_persistence
    ),
    data = species_mod_df,
    pred_probs = species_mod_df %>%
      select(atlas_block, p_colonization, p_extinction, p_persistence)
  ))
}

### --- APPLY, INSPECT --- ###

rbwo_log_results <- Model_three_logistic("Red-bellied Woodpecker")

summary(rbwo_log_results$models$colonization)
summary(rbwo_log_results$models$extinction)
summary(rbwo_log_results$models$persistence)

head(rbwo_log_results$data)










































# Create binary response variable, subset by transition state variables

# --- RCKI ---
rcki_wibba_summary <- rcki_wibba_summary %>%
  mutate(
    state = case_when(
      det_Atlas1 == 0 & det_Atlas2 == 0 ~ "A",
      det_Atlas1 == 0 & det_Atlas2 == 1 ~ "C",
      det_Atlas1 == 1 & det_Atlas2 == 0 ~ "E",
      det_Atlas1 == 1 & det_Atlas2 == 1 ~ "P"
    )
  )

rcki_colo_data <- rcki_wibba_summary %>%
  filter(state %in% c("A", "C")) %>%
  mutate(colo_rcki = as.numeric(state == "C"))

rcki_ext_data <- rcki_wibba_summary %>%
  filter(state %in% c("P", "E")) %>%
  mutate(
    ext_rcki = as.numeric(state == "E"),
    perst_rcki = as.numeric(state == "P")  # used also for persistence model
  )

# --- RBWO ---
rbwo_wibba_summary <- rbwo_wibba_summary %>%
  mutate(
    state = case_when(
      det_Atlas1 == 0 & det_Atlas2 == 0 ~ "A",
      det_Atlas1 == 0 & det_Atlas2 == 1 ~ "C",
      det_Atlas1 == 1 & det_Atlas2 == 0 ~ "E",
      det_Atlas1 == 1 & det_Atlas2 == 1 ~ "P"
    )
  )

rbwo_colo_data <- rbwo_wibba_summary %>%
  filter(state %in% c("A", "C")) %>%
  mutate(colo_rbwo = as.numeric(state == "C"))

rbwo_ext_data <- rbwo_wibba_summary %>%
  filter(state %in% c("P", "E")) %>%
  mutate(
    ext_rbwo = as.numeric(state == "E"),
    perst_rbwo = as.numeric(state == "P")
  )


## LOGISTIC REGRESSIONS ##
# (temporarily not including effort/sr_diff)

# --- RCKI Models ---

# Colonization
colonization_mod_rcki <- glm(
  colo_rcki ~ pa_percent_z + summer_tmax_diff_z + water_open_z + developed_open_z + 
    developed_low_z + developed_med_z + developed_high_z + barren_land_z + 
    forest_deciduous_z + forest_evergreen_z + forest_mixed_z + shrub_scrub_z + 
    grassland_z + pasture_z + cropland_z + wetlands_woody_z + wetlands_herb_z,
  data = rcki_colo_data,
  family = binomial(link = "logit")
)

# Extinction
extinction_mod_rcki <- glm(
  ext_rcki ~ pa_percent_z + summer_tmax_diff_z + water_open_z + developed_open_z + 
    developed_low_z + developed_med_z + developed_high_z + barren_land_z + 
    forest_deciduous_z + forest_evergreen_z + forest_mixed_z + shrub_scrub_z + 
    grassland_z + pasture_z + cropland_z + wetlands_woody_z + wetlands_herb_z,
  data = rcki_ext_data,
  family = binomial(link = "logit")
)

# Persistence
persistence_mod_rcki <- glm(
  perst_rcki ~ pa_percent_z + summer_tmax_diff_z + water_open_z + developed_open_z + 
    developed_low_z + developed_med_z + developed_high_z + barren_land_z + 
    forest_deciduous_z + forest_evergreen_z + forest_mixed_z + shrub_scrub_z + 
    grassland_z + pasture_z + cropland_z + wetlands_woody_z + wetlands_herb_z,
  data = rcki_ext_data,
  family = binomial(link = "logit")
)

# Summaries
summary(colonization_mod_rcki)
summary(extinction_mod_rcki)
summary(persistence_mod_rcki)


# --- RBWO Models ---

# Colonization
colonization_mod_rbwo <- glm(
  colo_rbwo ~ pa_percent_z + summer_tmax_diff_z + water_open_z + developed_open_z + 
    developed_low_z + developed_med_z + developed_high_z + barren_land_z + 
    forest_deciduous_z + forest_evergreen_z + forest_mixed_z + shrub_scrub_z + 
    grassland_z + pasture_z + cropland_z + wetlands_woody_z + wetlands_herb_z,
  data = rbwo_colo_data,
  family = binomial(link = "logit")
)

# Extinction
extinction_mod_rbwo <- glm(
  ext_rbwo ~ pa_percent_z + summer_tmax_diff_z + water_open_z + developed_open_z + 
    developed_low_z + developed_med_z + developed_high_z + barren_land_z + 
    forest_deciduous_z + forest_evergreen_z + forest_mixed_z + shrub_scrub_z + 
    grassland_z + pasture_z + cropland_z + wetlands_woody_z + wetlands_herb_z,
  data = rbwo_ext_data,
  family = binomial(link = "logit")
)

# Persistence
persistence_mod_rbwo <- glm(
  perst_rbwo ~ pa_percent_z + summer_tmax_diff_z + water_open_z + developed_open_z + 
    developed_low_z + developed_med_z + developed_high_z + barren_land_z + 
    forest_deciduous_z + forest_evergreen_z + forest_mixed_z + shrub_scrub_z + 
    grassland_z + pasture_z + cropland_z + wetlands_woody_z + wetlands_herb_z,
  data = rbwo_ext_data,
  family = binomial(link = "logit")
)

# Summaries
summary(colonization_mod_rbwo)
summary(extinction_mod_rbwo)
summary(persistence_mod_rbwo)


## PRELIM PLOTS ##

# CASC plot with RCKI, RBWO transition state and % pa across wibba blocks

# Add species labels and raw pa_percent to both datasets
rcki_plot_data <- rcki_wibba_summary %>%
  mutate(species = "RCKI", pa_percent = pa_percent) %>%
  filter(!is.na(state))

rbwo_plot_data <- rbwo_wibba_summary %>%
  mutate(species = "RBWO", pa_percent = pa_percent) %>%
  filter(!is.na(state))

# combine
combo_plot_data <- bind_rows(rcki_plot_data, rbwo_plot_data)

