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

### --- GRAPHICAL ANALYSIS --- ###

# Check distribution of differences in SR between RLL and DNR comp block sets 

srdiff_hist_data <- bind_rows(
  sr_summary_comp %>%
    dplyr::select(sr_diff) %>%
    mutate(source = "rll"),
  
  sr_summary_dnrcomp %>%
    dplyr::select(sr_diff) %>%
    mutate(source = "dnr")
)


ggplot(srdiff_hist_data, aes(x = sr_diff, fill = source)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_viridis_d(option = "plasma") +
  theme_bw() +
  labs(
    title = "Distribution of SR Differences",
    x = "Species Richness Difference (sr_diff)",
    y = "Count"
  )



### --- SR DIFFERENCE --- ###

# Join to summary df
wibba_covars_rll <- wibba_covars_raw %>%
  left_join(climate_summary_covars, by = "atlas_block")





######### NEW 2.27.26] ###############

covars_raw_rll <- covars_raw_rll %>% # add spp richness effort proxy
  mutate(sr_diff = sr_Atlas2 - sr_Atlas1)

# write.csv(covars_raw_rll, "data/summaries/covars_raw_rll.csv", row.names = FALSE)

