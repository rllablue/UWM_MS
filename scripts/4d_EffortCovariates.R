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

# Regressing SR on other ecological covariates to untangle ecological from survey effort-driven signal in richness as effort/completeness proxy
# Will include residual sr for Atlas 1, Atlas 2
# Change in effort between two periods only captures some of effort incongruence because...
# ...SR in Atlas 1 has disproportionate/asymmetric pull on what state transition will be assigned per block per species
# To discern source of effort bias, need to "baseline" effort to understand 'under-' and 'over-surveying'

# Blocks that were poorly sampled in Atlas 1 will have low residuals, and this will down-weight apparent “extinctions.”
# Blocks with high residuals in Atlas 2 will up-weight colonization likelihoods, since high effort means it’s unlikely you missed a true detection.

# Fit models on un-standardized SR
mod_srA1 <- lm(sr_Atlas1 ~ pa_z + developed_total_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z, data = wibba_summary_comp)
mod_srA2 <- lm(sr_Atlas2 ~ pa_z + developed_total_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z, data = wibba_summary_comp)

# Extract residuals, standardize, add to summary
wibba_summary_comp <- wibba_summary_comp %>%
  mutate(
    srA1_resid = resid(mod_srA1),
    srA2_resid = resid(mod_srA2),
    srA1_resid_z = as.numeric(scale(srA1_resid)),
    srA2_resid_z = as.numeric(scale(srA2_resid))
  )
