#############
### SETUP ###
#############

### --- LOAD PACKAGES --- ###

library(dplyr)
library(ggplot2)
library(lubridate)
library(magrittr)
library(tidyr)
library(purrr)
library(units)
library(broom)
library(nnet)


### --- DATA PREP --- ###

# Add species alpha codes to data to support flexible naming system during modeling

breeders_alpha <- read.csv("data/species/wibba_breeders_alpha.csv", stringsAsFactors = FALSE) 

breeders_zf_summary <- breeders_zf_summary %>%
  left_join(breeders_alpha, by = "common_name")


#################################
### 3-STATE MULTINOMIAL MODEL ###
#################################

### Modeling all states except "Absence" because lesser relevance/more signal 
# dampening effect by calculating probabilities for blocks where species 
# 'aren't already are.' Set "Presence" as base-state of comparison for 
# "Colonization" and "Extinction" states.

### --- FLEXIBLE MODEL --- ###

Model_multinomial <- function(species_name, 
                              alpha_lookup = breeders_alpha,
                              zf_summary = breeders_zf_summary,
                              covariates = wibba_modeling_comp,
                              baseline_state = "Persistence",
                              min_cov_pct = 0.01,
                              min_change_pct = 0.01
                              ) {
 
  
 # Lookup species alpha code
  alpha_code <- alpha_lookup %>%
    filter(common_name == species_name) %>%
    pull(alpha_code)
  
  if (length(alpha_code) == 0) {
    stop(paste0("No alpha_code found for ", species_name, ". Check breeders_alpha."))
  } 
   

  # Filter, aggregate detection data 
  species_dets <- zf_summary %>% # filter full species zf data
    filter(common_name == species_name) %>%
    filter(transition_state != "Absence")
  
  species_mod_df <- species_dets %>% # combine zf, covariate data
    left_join(covariates, by = "atlas_block")
  
  
  # Determine which landcover covariates to keep per species by thresholds
  # Baseline landcover (z reference year: 2008)
  baseline_landcov <- grep("_z_08$", names(species_mod_df), value = TRUE)
  keep_baseline <- baseline_landcov[
    sapply(species_mod_df[baseline_landcov], function(x) mean(x, na.rm = TRUE) >= min_cov_pct)
  ]
  
  # Landcover change (z difference reference years: 2015 - 1995)
  change_landcov <- grep("_z_9515$", names(species_mod_df), value = TRUE)
  keep_change <- change_landcov[
    sapply(species_mod_df[change_landcov], function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) >= min_change_pct)
  ]
  
  # Model specifications
  ### WIP WIP WIP: add climate ###
  fixed_covars <- c("pa_z", "srA1_resid", "srA2_resid") # aggregate covariates
  final_covars <- c(fixed_covars, keep_baseline, keep_change)
  
  species_mod_df <- species_mod_df %>%
    mutate(
      transition_state = factor(transition_state), # factorize response
      transition_state = relevel(transition_state, ref = baseline_state) # set baseline state
    )

  fmla <- as.formula(
    paste("transition_state ~", paste(final_covars, collapse = " + "))   # flexible formula construction
  ) 
  
  species_mod <- nnet::multinom(fmla, data = species_mod_df) # fit model
  species_mod_df$pred_probs <- predict(species_mod, type = "probs") # add predicted probabilities
  
  # Model summary objects
  model_summary <- summry(species_mod)
  model_coef <- coef(species_mod)
  model_odds <- exp(model_coef)
  
  # Return, label all results with flexible naming convention
  result <- list(
    species = species_name,
    alpha = alpha_code,
    model = species_mod,
    data = species_mod_df,
    pred_probs = species_mod_df$pred_probs,
    used_covars = final_covars,
    summary = model_summary,
    coef = model_coef,
    odds = model_odds
  )
  
  assign(paste0(alpha_code, "_multinom"), result, envir = .GlobalEnv)
  message(paste0("Model for ", species_name, " (", alpha_code, ") completed and saved as ", alpha_code, "_multinom"))
  
  return(result)
}


### --- APPLY, INSPECT --- ###

rbwo_multinom <- Model_multinomial(
  species_name = "Red-bellied Woodpecker", # select species to run
  min_cov_pct = 0.02,    # only include landcover >=2% of blocks
  min_change_pct = 0.01  # only include change >=1% across blocks
)

summary(rbwo_multinom$model)
rbwo_multinom$used_covars
coef(rbwo_multinom$model)
exp(coef(rbwo_multinom$model))

head(rbwo_multinom$data) # view combined data

head(rbwo_multinom$pred_probs) # predicted probs per block




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