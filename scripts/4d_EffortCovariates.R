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

# Fit models on un-standardized SR
# Use grouped land cover covariates, climate proxy covariates for ease
mod_srA1 <- lm(sr_Atlas1 ~ pa_z + water_open_z + shrub_scrub_z + grassland_z + developed_total_z + forest_total_z + wetlands_total_z + lat_z + lon_z, data = wibba_modeling_comp)
mod_srA2 <- lm(sr_Atlas2 ~ pa_z + water_open_z + shrub_scrub_z + grassland_z + developed_total_z + forest_total_z + wetlands_total_z + lat_z + lon_z, data = wibba_modeling_comp)

# Extract residuals, standardize, add to summaries
sr_resids <- wibba_modeling_comp %>%
  mutate(
    srA1_resid = resid(mod_srA1),
    srA2_resid = resid(mod_srA2),
    srA1_resid_z = as.numeric(scale(srA1_resid)),
    srA2_resid_z = as.numeric(scale(srA2_resid))
  ) %>%
  dplyr::select(atlas_block, srA1_resid, srA2_resid, srA1_resid_z, srA2_resid_z)

sr_summary_comp <- sr_summary_comp %>%
  left_join(sr_resids, by = "atlas_block")

wibba_modeling_comp <- wibba_modeling_comp %>%
  mutate(
    srA1_resid_z = sr_resids$srA1_resid_z,
    srA2_resid_z = sr_resids$srA2_resid_z
  ) %>%
  dplyr::select(-sr_Atlas1, -sr_Atlas2)

