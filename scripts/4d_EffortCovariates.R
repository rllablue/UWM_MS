#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)


###################################
### SPECIES RICHNESS / 'EFFORT' ###
###################################

### Create effort proxy by regressing SR on other ecological covariates to untangle 
# ecological from survey effort-driven signal in richness. Include SR residuals for 
# both periods: Atlas 1 has disproportionate/asymmetric pull on what state transition 
# will be assigned per block per species; to discern source of effort bias, need 
# "baseline" effort to understand 'under-' and 'over-surveying' illuminated by 
# inclusion of resids for Atlas 2:
### Blocks poorly sampled in Atlas 1 will have low residuals, which will down-weight 
# apparent “extinction"; blocks with high residuals in Atlas 2 will up-weight 
# "colonization" since high effort means lowers probability of true non-detection.

# Put raw SR into the modeling df
wibba_modeling_covars <- wibba_modeling_covars %>%
  left_join(
    sr_summary_comp %>%
      dplyr::select(atlas_block, sr_Atlas1, sr_Atlas2),
      by = "atlas_block"
  )


# Extract relevant covars for residual regression, z-scale
covars_for_sr <- wibba_modeling_covars %>%
  dplyr::select(
    atlas_block,
    pa_percent, water_open_2008, shrub_scrub_2008, grassland_2008,
    developed_total_2008, forest_total_2008, wetlands_total_2008,
    lat, lon
  )

covars_for_sr <- covars_for_sr %>%
  mutate(
    across(
      .cols = -atlas_block,   
      .fns = ~ as.numeric(scale(.x)),    
      .names = "{.col}_z"
    )
  )

covars_for_sr <- covars_for_SR %>%
  left_join(
  sr_summary_comp %>%
    dplyr::select(atlas_block, sr_Atlas1, sr_Atlas2),
  by = "atlas_block"
)

# Fit models on un-standardized SR
# Use grouped land cover covariates, climate proxy covariates for ease
mod_srA1 <- lm(sr_Atlas1 ~ pa_percent_z + water_open_2008_z + shrub_scrub_2008_z + grassland_2008_z + developed_total_2008_z + forest_total_2008_z + wetlands_total_2008_z + lat_z + lon_z, data = covars_for_sr)
mod_srA2 <- lm(sr_Atlas2 ~ pa_percent_z + water_open_2008_z + shrub_scrub_2008_z + grassland_2008_z + developed_total_2008_z + forest_total_2008_z + wetlands_total_2008_z + lat_z + lon_z, data = covars_for_sr)

# Extract residuals, add to summaries
sr_resids <- covars_for_sr %>%
  mutate(
    srA1_resid = resid(mod_srA1),
    srA2_resid = resid(mod_srA2)
  ) %>%
  dplyr::select(atlas_block, srA1_resid, srA2_resid)

sr_summary_comp <- sr_summary_comp %>%
  left_join(sr_resids, by = "atlas_block")

wibba_modeling_covars <- wibba_modeling_covars %>%
  left_join(sr_resids, by = "atlas_block") %>%
  dplyr::select(-sr_Atlas1, -sr_Atlas2)