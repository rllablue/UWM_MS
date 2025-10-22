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
wibba_modeling_comp <- wibba_modeling_comp %>%
  left_join(
    sr_summary_comp %>%
      dplyr::select(atlas_block, sr_Atlas1, sr_Atlas2),
      by = "atlas_block"
  )

# Fit models on un-standardized SR ################################## DECIDE WHICH COVARS TO PUT IN HERE ##################
mod_srA1 <- lm(sr_Atlas1 ~ pa_z + developed_lower_z_08 + developed_upper_z_08 + forest_total_z + grass_pasture_crop_z + wetlands_total_z, data = wibba_modeling_comp)
mod_srA2 <- lm(sr_Atlas2 ~ pa_z + developed_lower_z + developed_upper_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z, data = wibba_modeling_comp)

# Extract residuals, standardize, add to summary
wibba_modeling_comp <- wibba_modeling_comp %>%
  mutate(
    srA1_resid = resid(mod_srA1),
    srA2_resid = resid(mod_srA2),
    srA1_resid_z = as.numeric(scale(srA1_resid)),
    srA2_resid_z = as.numeric(scale(srA2_resid))
  )

# Add resids to sr_summary_comp