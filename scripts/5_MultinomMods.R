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
library(broom)
library(nnet)



#################################
### 4-STATE MULTINOMIAL MODEL ###
#################################

### RCKI ###
# Combine species detection data and modeling covariates from summary data
rcki_multinom_mod <- rcki_dets_summary %>%
  left_join(
    wibba_summary_comp %>%
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

# Factorize response variable, set model baseline/ref
rcki_multinom_mod <- rcki_multinom_mod %>%
  mutate(
    transition_state = factor(transition_state),
    transition_state = relevel(transition_state, ref = "Absence")  
  )

# Set up model
mod_rcki <- multinom(
  transition_state ~ pa_z + lat_z + lon_z + srA1_resid + srA2_resid +
    developed_total_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z,
  data = rcki_multinom_mod
)

summary(mod_rcki)
coef(mod_rcki)
exp(coef(mod_rcki))

rcki_multinom_mod$pred_probs <- predict(mod_rcki, type = "probs")



### RBWO ###
# Combine species detection data and modeling covariates from summary data
rbwo_multinom_mod <- rbwo_dets_summary %>%
  left_join(
    wibba_summary_comp %>%
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

# Factorize response variable, set model baseline/ref
rbwo_multinom_mod <- rbwo_multinom_mod %>%
  mutate(
    transition_state = factor(transition_state),
    transition_state = relevel(transition_state, ref = "Absence")  
  )

# Set up model
mod_rbwo <- multinom(
  transition_state ~ pa_z + lat_z + lon_z + srA1_resid + srA2_resid +
    developed_total_z + forest_total_z + grass_pasture_crop_z + wetlands_total_z,
  data = rbwo_multinom_mod
)

summary(mod_rbwo)
coef(mod_rbwo)
exp(coef(mod_rbwo))

rbwo_multinom_mod$pred_probs <- predict(mod_rbwo, type = "probs")




################
### PLOTTING ###
################

# Tidy and exponentiate coefficients (odds ratios)
rbwo_tidy <- tidy(mod_rbwo, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(species = "Red-bellied Woodpecker")

rcki_tidy <- tidy(mod_rcki, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(species = "Ruby-crowned Kinglet")

# Combine both species' results
mod_combined <- bind_rows(rbwo_tidy, rcki_tidy)


mod_combined <- mod_combined %>%
  mutate(
    term = recode(term,
                  pa_z = "Protected Area Coverage (%)",
                  lat_z = "Latitude",
                  lon_z = "Longitude",
                  srA1_resid = "Species Richness Atlas 1 (resid.)",
                  srA2_resid = "Species Richness Atlas 2 (resid.)",
                  developed_total_z = "Developed Land (km2)",
                  forest_total_z = "Forest cover (km2)",
                  grass_pasture_crop_z = "Grass/Pasture/Crop (km2)",
                  wetlands_total_z = "Wetland cover (km2)"
    )
  )

ggplot(mod_combined, aes(x = term, y = estimate, color = y.level)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.6), width = 0.2
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_y_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10),
    labels = c("0.1", "0.5", "1", "2", "5", "10")
  ) +
  coord_flip() +
  facet_wrap(~species, ncol = 2, scales = "free_x") +
  labs(
    title = "Effects of Environmental Covariates on Species Transition States",
    subtitle = "Odds ratios (relative to Absence) from multinomial models",
    x = NULL,
    y = "Odds Ratio (log scale)",
    color = "Transition State"
  ) +
  theme_minimal(base_size = 13) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 15)
  )