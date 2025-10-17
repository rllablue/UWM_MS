## 4 LOGISTIC REG MODELS (CASC POSTER) ##
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

