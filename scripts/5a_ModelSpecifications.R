#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

# Core
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(stringr)

# Diagnostics
library(nnet)
library(stats)
library(broom)
library(corrplot)
library(car)
library(performance)
library(DescTools)
library(Metrics)
library(pROC)
library(rsample)
library(glmnet)
library(MuMIn)

# Visualization
library(ggplot2)
library(ggfortify)
library(viridis)



## DFS ##

wibba_covars_raw
breeders_zf_summary


# Full, species-specific detection df
spp_name <- "Red-bellied Woodpecker"

spp_mod_data_full <- breeders_zf_summary %>%
  filter(common_name == spp_name) %>%
  left_join(
    wibba_covars_raw,
    by = "atlas_block"
  )

# Species-, state-specific detection dfs

### Separate data into bins where A2 detection is either 1 or 0; not necessarily
# in precise probabilities between colo, pres, abs, ext, but more what promotes
# 'new' v. 'continued' colonization, ie. what promotes det = 1, as opposed to 
# what promotes det = 0. 
### Combining absence with extinction will also solve problem of signal dampening 
# when absence data is included along with (both kinds) od presence data, ie. we're 
# only looking at relationship to pa (and other covars) for places where a species
# already resides.

# 'Colonization-Persistence' df
data_colper <- spp_mod_data_full %>%
  filter(transition_state %in% c("Colonization", "Persistence"))

data_colper_z <- data_colper %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))
  


# 'Extinction-Absence' df
data_extabs <- spp_mod_data_full %>%
  filter(transition_state %in% c("Extinction", "Absence"))

data_extabs_z <- data_extabs %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))



# 'Colonization-Absence' df
data_colabs <- spp_mod_data_full %>%
  filter(transition_state %in% c("Colonization", "Absence"))

data_colabs_z <- data_colabs %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))



# 'Extinction-Persistence' df
data_extper <- spp_mod_data_full %>%
  filter(transition_state %in% c("Extinction", "Persistence"))

data_extper_z <- data_extper %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))




# Assign binomial state identifiers for model
data_colper_z <- data_colper_z %>%
  mutate(
    col_per = ifelse(transition_state %in% c("Colonization"), 1, 0),
  )

data_extabs_z <- data_extabs_z %>%
  mutate(
    ext_abs = ifelse(transition_state %in% c("Extinction"), 1, 0)
  )

data_colabs_z <- data_colabs_z %>%
  mutate(
    col_abs = ifelse(transition_state %in% c("Colonization"), 1, 0)
  )

data_extper_z <- data_extper_z %>%
  mutate(
    ext_per = ifelse(transition_state %in% c("Extinction"), 1, 0)
  )





############ RED-BELLIED WOODPECKER ############################################
# Col-Per
rbwo_mod_colper <- glm(col_per ~  pa_percent_z +
                               developed_lower_base_z + developed_upper_base_z + 
                               forest_deciduous_base_z + forest_mixed_base_z +
                               pasture_crop_base_z + 
                               forest_total_diff_z + developed_total_diff_z +
                               tmax_38yr_z + prcp_38yr_z,
                             data = data_colper_z, family = binomial)
summary(rbwo_mod_colper)
autoplot(rbwo_mod_colper)

# Ext-Abs
rbwo_mod_extabs <- glm(ext_abs ~ pa_percent_z +
                      developed_lower_base_z + developed_upper_base_z + 
                      forest_deciduous_base_z + forest_mixed_base_z +
                      pasture_crop_base_z + 
                      forest_total_diff_z + developed_total_diff_z +
                      tmax_38yr_z + prcp_38yr_z,
                    data = data_extabs_z, family = binomial)

summary(rbwo_mod_extabs)
autoplot(rbwo_mod_extabs)

# Col-Abs
rbwo_mod_colabs <- glm(col_abs ~  pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_colabs_z, family = binomial)
summary(rbwo_mod_colabs)
autoplot(rbwo_mod_colabs)

# Ext-Per
rbwo_mod_extper <- glm(ext_per ~  pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_extper_z, family = binomial)
summary(rbwo_mod_extper)
autoplot(rbwo_mod_extper)











#### VISUALIZE ####

rbwo_vars_plot <- c(
  "pa_percent_z",
  "developed_lower_base_z",
  "forest_mixed_base_z",
  "forest_deciduous_base_z",
  "pasture_crop_base_z",
  "tmax_38yr_z",
  "prcp_38yr_z",
  "developed_total_diff_z"
)


make_pred_df <- function(var_name, model, data, n = 100) {
  # baseline at mean for all predictors
  base_vals <- data %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  # vary only this variable
  newdata <- base_vals[rep(1, n), ]
  newdata[[var_name]] <- seq(
    min(data[[var_name]], na.rm = TRUE),
    max(data[[var_name]], na.rm = TRUE),
    length.out = n
  )
  
  # predicted probability
  newdata$pred_prob <- predict(model, newdata = newdata, type = "response")
  newdata$variable <- var_name
  newdata$x_value <- newdata[[var_name]]
  
  newdata %>% dplyr::select(variable, x_value, pred_prob)
}


### PLOT 
facet_labels_col <- c(
  pa_percent_z = "Protected Area (%)***",
  tmax_38yr_z = "Max Temp (38-year avg)***",
  prcp_38yr_z = "Precipitation (38-year avg)",
  forest_deciduous_base_z = "Deciduous Forest***",
  forest_mixed_base_z = "Mixed Forest***",
  pasture_crop_base_z = "Pasture/Cropland*",
  developed_lower_base_z = "Open/Lower Development***",
  developed_total_diff_z = "Difference in Total Developed Land (%)"
)


facet_labels_ext <- c(
  pa_percent_z = "Protected Area (%)",
  tmax_38yr_z = "Max Temp (38-year avg)*",
  prcp_38yr_z = "Precipitation (38-year avg)",
  forest_deciduous_base_z = "Deciduous Forest",
  forest_mixed_base_z = "Mixed Forest",
  pasture_crop_base_z = "Pasture/Cropland .",
  developed_lower_base_z = "Open/Lower Development",
  developed_total_diff_z = "Difference in Total Developed Land (%)"
)


# COL-PER
rbwo_predplot_colper_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                        model = rbwo_mod_colper, 
                        data = data_colper_z)


ggplot(rbwo_predplot_colper_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_col)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across Blocks"
  )


# EXT-ABS
rbwo_predplot_ext_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                             model = rbwo_mod_ext, 
                             data = spp_mod_data_ext_z)


ggplot(rbwo_predplot_ext_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_ext)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across Blocks"
  )

# COL-ABS
rbwo_predplot_colabs_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_colabs, 
                                   data = data_colabs_z)


ggplot(rbwo_predplot_colabs_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_col)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across Blocks"
  )

# EXT-PER
rbwo_predplot_extper_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_extper, 
                                   data = data_extper_z)


ggplot(rbwo_predplot_extper_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_ext)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across Blocks"
  )





########################## BOBOLINK ###########################################

bobo_mod_col <- glm(col_per ~  pa_percent_z +
                      developed_lower_base_z + developed_upper_base_z + 
                      forest_total_base_z +
                      pasture_crop_base_z + grassland_base_z +
                      pasture_crop_diff_z + grassland_diff_z +
                      tmax_38yr_z + prcp_38yr_z,
                    data = spp_mod_data_col_z, family = binomial)
summary(bobo_mod_col)
autoplot(bobo_mod_col)


bobo_mod_ext <- glm(ext_abs ~ pa_percent_z +
                      developed_lower_base_z + developed_upper_base_z + 
                      forest_total_base_z +
                      pasture_crop_base_z + grassland_base_z +
                      pasture_crop_diff_z + grassland_diff_z +
                      tmax_38yr_z + prcp_38yr_z,
                    data = spp_mod_data_ext_z, family = binomial)

summary(bobo_mod_ext)
autoplot(bobo_mod_ext)




#### VISUALIZE ####

bobo_vars_plot <- c(
  "pa_percent_z",
  "pasture_crop_base_z",
  "grassland_diff_z",
  "prcp_38yr_z",
  "tmax_38yr_z"
)


make_pred_df <- function(var_name, model, data, n = 100) {
  # baseline at mean for all predictors
  base_vals <- data %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  # vary only this variable
  newdata <- base_vals[rep(1, n), ]
  newdata[[var_name]] <- seq(
    min(data[[var_name]], na.rm = TRUE),
    max(data[[var_name]], na.rm = TRUE),
    length.out = n
  )
  
  # predicted probability
  newdata$pred_prob <- predict(model, newdata = newdata, type = "response")
  newdata$variable <- var_name
  newdata$x_value <- newdata[[var_name]]
  
  newdata %>% dplyr::select(variable, x_value, pred_prob)
}


### PLOT 
facet_labels_col <- c(
  pa_percent_z = "Protected Area (%)",
  tmax_38yr_z = "Max Temp (38-year avg).",
  prcp_38yr_z = "Precipitation (38-year avg).",
  pasture_crop_base_z = "Pasture + Crop",
  grassland_diff_z = "Grassland (chnage)"
)


facet_labels_ext <- c(
  pa_percent_z = "Protected Area (%)***",
  tmax_38yr_z = "Max Temp (38-year avg).",
  prcp_38yr_z = "Precipitation (38-year avg).",
  pasture_crop_base_z = "Pasture + Crop**",
  grassland_diff_z = "Grassland (change)"
)


# COLO
bobo_predplot_col_df <- map_dfr(bobo_vars_plot, make_pred_df, 
                                model = bobo_mod_col, 
                                data = spp_mod_data_col_z)


ggplot(bobo_predplot_col_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_col)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across Blocks"
  )


# EXT
bobo_predplot_ext_df <- map_dfr(bobo_vars_plot, make_pred_df, 
                                model = bobo_mod_ext, 
                                data = spp_mod_data_ext_z)


ggplot(bobo_predplot_ext_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_ext)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across Blocks"
  )



############################## RCKI ############################################

rcki_mod_col <- glm(col_per ~  pa_percent_z +
                      developed_lower_base_z +
                      forest_evergreen_base_z + forest_mixed_base_z +
                      wetlands_total_base_z + 
                      tmax_38yr_z + prcp_38yr_z,
                    data = spp_mod_data_col_z, family = binomial)
summary(rcki_mod_col)
autoplot(rcki_mod_col)


rcki_mod_ext <- glm(ext_abs ~ pa_percent_z +
                      developed_lower_base_z +
                      forest_evergreen_base_z + forest_mixed_base_z +
                      wetlands_total_base_z + 
                      tmax_38yr_z + prcp_38yr_z,
                    data = spp_mod_data_ext_z, family = binomial)

summary(rcki_mod_ext)
autoplot(rcki_mod_ext)




#### VISUALIZE ####

rcki_vars_plot <- c(
  "pa_percent_z",
  "forest_mixed_base_z",
  "wetlands_total_base_z",
   "prcp_38yr_z",
  "tmax_38yr_z"
)


make_pred_df <- function(var_name, model, data, n = 100) {
  # baseline at mean for all predictors
  base_vals <- data %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  # vary only this variable
  newdata <- base_vals[rep(1, n), ]
  newdata[[var_name]] <- seq(
    min(data[[var_name]], na.rm = TRUE),
    max(data[[var_name]], na.rm = TRUE),
    length.out = n
  )
  
  # predicted probability
  newdata$pred_prob <- predict(model, newdata = newdata, type = "response")
  newdata$variable <- var_name
  newdata$x_value <- newdata[[var_name]]
  
  newdata %>% dplyr::select(variable, x_value, pred_prob)
}


### PLOT 
facet_labels_col <- c(
  pa_percent_z = "Protected Area (%)",
  tmax_38yr_z = "Max Temp (38-year avg)***",
  prcp_38yr_z = "Precip (38-year avg)",
  wetlands_total_base_z = "Wetlands (herb + woody)",
  forest_mixed_base_z = "Mixed Forest ."
)


facet_labels_ext <- c(
  pa_percent_z = "Protected Area (%)",
  tmax_38yr_z = "Max Temp (38-year avg)***",
  prcp_38yr_z = "Precip (38-year avg) .",
  wetlands_total_base_z = "Wetlands (herb + woody)***",
  forest_mixed_base_z = "Mixed Forest ."
)


# COLO
rcki_predplot_col_df <- map_dfr(rcki_vars_plot, make_pred_df, 
                                model = rcki_mod_col, 
                                data = spp_mod_data_col_z)


ggplot(rcki_predplot_col_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_col)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across Blocks"
  )


# EXT
rcki_predplot_ext_df <- map_dfr(rcki_vars_plot, make_pred_df, 
                                model = rcki_mod_ext, 
                                data = spp_mod_data_ext_z)


ggplot(rcki_predplot_ext_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_labels_ext)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across Blocks"
  )














################################################################################

### Covariate Selection
# Diagnostics: correlations, vif, aic ranking and selection




spp_mod_data_col_z <- data.frame()
spp_mod_data_ext_z <- data.frame()