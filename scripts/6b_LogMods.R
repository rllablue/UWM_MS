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

spp_mod_data <- spp_mod_data %>% # create state identifiers
  mutate(
    is_colonization = ifelse(transition_state == "Colonization", 1, 0),
    is_extinction   = ifelse(transition_state == "Extinction", 1, 0),
    is_persistence  = ifelse(transition_state == "Presence", 1, 0)
  )


###
rbwo_mod_colonization <- glm(col ~ 
                               water_open_2008_z + barren_land_2008_z + shrub_scrub_2008_z +
                                developed_total_2008_z + forest_deciduous_2008_z + forest_evergreen_2008_z + forest_mixed_2008_z +
                                wetlands_total_2008_z + grassland_2008_z + pasture_crop_2008_z + water_open_diff_z +
                                + forest_total_diff_z + wetlands_total_diff_z + developed_total_diff_z,
                             data = spp_mod_data, family = binomial)

rbwo_mod_extinction <- glm(ext ~ 
                             water_open_2008_z + barren_land_2008_z + shrub_scrub_2008_z +
                             developed_total_2008_z + forest_deciduous_2008_z + forest_evergreen_2008_z + forest_mixed_2008_z +
                             wetlands_total_2008_z + grassland_2008_z + pasture_crop_2008_z + water_open_diff_z +
                             + forest_total_diff_z + wetlands_total_diff_z + developed_total_diff_z,
                           data = spp_mod_data, family = binomial)

rbwo_mod_persistence <- glm(per ~
                              water_open_2008_z + barren_land_2008_z + shrub_scrub_2008_z +
                              developed_total_2008_z + forest_deciduous_2008_z + forest_evergreen_2008_z + forest_mixed_2008_z +
                              wetlands_total_2008_z + grassland_2008_z + pasture_crop_2008_z + water_open_diff_z +
                              + forest_total_diff_z + wetlands_total_diff_z + developed_total_diff_z,
                            data = spp_mod_data, family = binomial)



# Helper function for tidy summary
summarize_model <- function(model, model_name) {
  tidy(model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      model = model_name,
      abs_estimate = abs(estimate)
    ) %>%
    arrange(p.value) %>%
    dplyr::select(model, term, estimate, std.error, p.value, abs_estimate)
}

# Apply to all three models
rbwo_cov_summary <- bind_rows(
  summarize_model(rbwo_mod_colonization, "Colonization"),
  summarize_model(rbwo_mod_extinction, "Extinction"),
  summarize_model(rbwo_mod_persistence, "Persistence")
)

# Rank within each model by importance
rbwo_cov_summary_ranked <- rbwo_cov_summary %>%
  group_by(model) %>%
  arrange(p.value, desc(abs_estimate)) %>%
  mutate(rank = row_number())

print(rbwo_cov_summary_ranked)


rbwo_top3_covars <- rbwo_cov_summary_ranked %>%
  filter(rank <= 3) %>%
  arrange(model, rank)

print(rbwo_top3_covars)





# FULL 
rbwo_mod_col_full <- glm(col ~ pa_percent_z + srA1_resid_z + srA2_resid_z + 
                               tmax_38yr_z + prcp_38yr_z +
                               developed_total_2008_z + forest_mixed_2008_z +
                               wetlands_total_2008_z,
                             data = spp_mod_data, family = binomial)

rbwo_mod_ext_full <- glm(ext ~ pa_percent_z + srA1_resid_z + srA2_resid_z + 
                             forest_mixed_2008_z + grassland_2008_z + forest_total_diff_z,
                           data = spp_mod_data, family = binomial)

rbwo_mod_per_full <- glm(per ~ pa_percent_z + srA1_resid_z + srA2_resid_z + 
                              water_open_2008_z + forest_mixed_2008_z + wetlands_total_2008_z,
                            data = spp_mod_data, family = binomial)

summary(rbwo_mod_col_full)
summary(rbwo_mod_ext_full)
summary(rbwo_mod_per_full)



##


land_covars_final <- c("water_open_2008_z", "barren_land_2008_z", "shrub_scrub_2008_z",
                       "developed_total_2008_z", "forest_deciduous_2008_z", "forest_evergreen_2008_z", "forest_mixed_2008_z",
                       "wetlands_total_2008_z", "grassland_2008_z", "pasture_crop_2008_z", "water_open_diff_z",
                       "forest_total_diff_z", "wetlands_total_diff_z", "developed_total_diff_z")


wip_land_covars_final <- c("developed_lower_2008_z",
                       "developed_upper_2008_z", "forest_total_2008_z",
                      "grassland_2008_z", "pasture_crop_2008_z"
                       )



clim_covars_final <- c("tmax_38yr_z", "prcp_38yr_z")



other_covars <- c("pa_percent_z", "srA1_resid_z", "srA2_resid_z")


all_covars <- c(wip_land_covars_final, clim_covars_final, other_covars)





# --- RBWO WI WIP WIP WIP WIP WIP
rbwo_mod_colonization <- glm(is_colonization ~ pa_percent_z + srA1_resid_z + srA2_resid_z +
                          developed_lower_2008_z + developed_upper_2008_z + forest_total_2008_z + pasture_crop_2008_z +
                          tmax_38yr_z + prcp_38yr_z,
                        data = spp_mod_data, family = binomial)

rbwo_mod_extinction <- glm(is_extinction ~ pa_percent_z + srA1_resid_z + srA2_resid_z +
                             developed_lower_2008_z + developed_upper_2008_z + forest_total_2008_z + pasture_crop_2008_z +
                             tmax_38yr_z + prcp_38yr_z,
                           data = spp_mod_data, family = binomial)

rbwo_mod_persistence <- glm(is_persistence ~ pa_percent_z + srA1_resid_z + srA2_resid_z +
                              developed_lower_2008_z + developed_upper_2008_z + forest_total_2008_z + pasture_crop_2008_z +
                              tmax_38yr_z + prcp_38yr_z,
                            data = spp_mod_data, family = binomial)



summary(rbwo_mod_colonization)
summary(rbwo_mod_extinction)
summary(rbwo_mod_persistence)




# Multinom 
spp_mod_data <- spp_mod_data %>%
  mutate(
    transition_state = case_when(
      col == 1 ~ "col",
      ext == 1 ~ "ext",
      per == 1 ~ "per",
      TRUE ~ NA_character_
    ),
    transition_state = factor(transition_state, levels = c("per", "col", "ext")) # per = baseline
  )

# --- Fit multinomial model ---
rbwo_multinom <- nnet::multinom(
  transition_state ~ pa_percent_z + srA1_resid_z + srA2_resid_z +
    developed_lower_2008_z + developed_upper_2008_z + forest_total_2008_z + pasture_crop_2008_z +
    tmax_38yr_z + prcp_38yr_z,
  data = spp_mod_data
)

# --- Summaries ---
summary(rbwo_multinom)

# --- Convert coefficients to odds ratios ---
exp(coef(rbwo_multinom))





############ AIC RANKING #######

# --- AIC comparison ---
aic_comparison <- tibble(
  model = c("multinom", "colonization", "extinction", "persistence"),
  AIC = c(AIC(rbwo_multinom),
          AIC(rbwo_mod_colonization),
          AIC(rbwo_mod_extinction),
          AIC(rbwo_mod_persistence))
)

aic_comparison


logLik(rbwo_multinom)
logLik(rbwo_mod_colonization)
logLik(rbwo_mod_extinction)
logLik(rbwo_mod_persistence)



# --- Extract multinomial coefficients ---
multi_coefs <- broom::tidy(rbwo_multinom, exponentiate = TRUE) %>%
  mutate(model = paste0("multinom_", y.level))

# --- Extract logit coefficients ---
logit_coefs <- bind_rows(
  broom::tidy(rbwo_mod_colonization, exponentiate = TRUE) %>% mutate(model = "col"),
  broom::tidy(rbwo_mod_extinction, exponentiate = TRUE) %>% mutate(model = "ext"),
  broom::tidy(rbwo_mod_persistence, exponentiate = TRUE) %>% mutate(model = "per")
)

# --- Combine and compare ---
coef_comparison <- bind_rows(multi_coefs, logit_coefs) %>%
  dplyr::select(model, term, estimate, std.error, p.value)

coef_comparison



library(pROC)

auc_col <- roc(spp_mod_data$col, spp_mod_data$p_col_multinom)$auc
auc_ext <- roc(spp_mod_data$ext, spp_mod_data$p_ext_multinom)$auc
auc_per <- roc(spp_mod_data$per, spp_mod_data$p_per_multinom)$auc

c(auc_col, auc_ext, auc_per)




########### BOBO WIP WIP WIP WIP WIP


bobo_mod_colonization <- glm(is_colonization ~ pa_percent_z + srA1_resid_z + srA2_resid_z +
                               developed_total_2008_z + forest_total_2008_z + pasture_crop_2008_z +
                               grassland_2008_z + tmax_38yr_z + prcp_38yr_z,
                             data = spp_mod_data, family = binomial)

bobo_mod_extinction <- glm(is_colonization ~ pa_percent_z + srA1_resid_z + srA2_resid_z +
                             developed_total_2008_z + forest_total_2008_z + pasture_crop_2008_z +
                             grassland_2008_z + tmax_38yr_z + prcp_38yr_z,
                           data = spp_mod_data, family = binomial)

bobo_mod_persistence <- glm(is_colonization ~ pa_percent_z + srA1_resid_z + srA2_resid_z +
                              developed_total_2008_z + forest_total_2008_z + pasture_crop_2008_z +
                              grassland_2008_z + tmax_38yr_z + prcp_38yr_z,
                            data = spp_mod_data, family = binomial)



summary(bobo_mod_colonization)
summary(bobo_mod_extinction)
summary(bobo_mod_persistence)






#############################################################################

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
  
  spp_mod_data <- spp_mod_data %>% # create state identifiers
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

