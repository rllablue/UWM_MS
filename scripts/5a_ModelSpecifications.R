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
library(pscl)
library(AICcmodavg)

# Visualization
library(ggplot2)
library(ggfortify)
library(viridis)
library(gt)


#######################
### RLL COMP BLOCKS ###
#######################

## DFS ##
breeders_zf_summary
wibba_covars_raw 


# Full, species-specific detection df
spp_name <- "Red-bellied Woodpecker"

# RLL COMP BLOCKS #
spp_mod_data_rll <- breeders_zf_summary %>%
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


# 'Colonization-Absence' df
data_colabs_rll <- spp_mod_data_rll %>%
  filter(transition_state %in% c("Colonization", "Absence"))

data_colabs_rll_z <- data_colabs_rll %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))



# 'Extinction-Persistence' df
data_extper_rll <- spp_mod_data_rll %>%
  filter(transition_state %in% c("Extinction", "Persistence"))

data_extper_rll_z <- data_extper_rll %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))



# Assign binomial state identifiers for model
data_colabs_rll_z <- data_colabs_rll_z %>%
  mutate(
    col_abs = ifelse(transition_state %in% c("Colonization"), 1, 0)
  )

data_extper_rll_z <- data_extper_rll_z %>%
  mutate(
    ext_per = ifelse(transition_state %in% c("Extinction"), 1, 0)
  )



############ RED-BELLIED WOODPECKER ############################################

# Col-Abs
rbwo_mod_colabs_rll <- glm(col_abs ~  pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_colabs_rll_z, family = binomial)
summary(rbwo_mod_colabs_rll)
autoplot(rbwo_mod_colabs_rll)

# Ext-Per
rbwo_mod_extper_rll <- glm(ext_per ~  pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_extper_rll_z, family = binomial)
summary(rbwo_mod_extper_rll)
autoplot(rbwo_mod_extper_rll)



#### VISUALIZE ####

rbwo_vars_plot <- c(
  "pa_percent_z",
  "developed_lower_base_z",
  "developed_upper_base_z",
  "forest_mixed_base_z",
  "forest_deciduous_base_z",
  "pasture_crop_base_z",
  "tmax_38yr_z",
  "prcp_38yr_z",
  "developed_total_diff_z",
  "forest_total_diff_z"
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


### PLOT ###

# COL-ABS
facet_colabs <- c(
  pa_percent_z = "Protected Area ***",
  developed_lower_base_z = "Open + Lower Development ***",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest ***",
  forest_deciduous_base_z = "Deciduous Forest ***",
  pasture_crop_base_z = "Pasture/Cropland *",
  tmax_38yr_z = "Max Temp ***",
  prcp_38yr_z = "Precipitation ***",
  developed_total_diff_z = "Difference in Total Developed Land *",
  forest_total_diff_z = "Difference in Total Forest"
)

rbwo_predplot_colabs_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_colabs, 
                                   data = data_colabs_z)


ggplot(rbwo_predplot_colabs_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_colabs)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across Blocks"
  )


# EXT-PER
facet_extper <- c(
  pa_percent_z = "Protected Area",
  developed_lower_base_z = "Open + Lower Development",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest",
  forest_deciduous_base_z = "Deciduous Forest",
  pasture_crop_base_z = "Pasture/Cropland .",
  tmax_38yr_z = "Max Temp *",
  prcp_38yr_z = "Precipitation",
  developed_total_diff_z = "Difference in Total Developed Land",
  forest_total_diff_z = "Difference in Total Forest"
)

rbwo_predplot_extper_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_extper, 
                                   data = data_extper_z)


ggplot(rbwo_predplot_extper_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_extper)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across Blocks"
  )



#######################
### DNR COMP BLOCKS ###
#######################

# DFS #
zf_dnr_compblocks <- dnr_compblocks %>% # zf data for dnr fair comp blocks (N. Anich)
  left_join(breeders_zf_summary, by = "atlas_block")

covars_dnr_compblocks_raw <- dnr_compblocks %>% # modeling covariates for dnr fair comp blocks
  left_join(wibba_covars_raw, by = "atlas_block")


# Full, species-specific detection df
spp_name <- "Red-bellied Woodpecker"

spp_mod_data_dnr <- zf_dnr_compblocks %>%
  filter(common_name == spp_name) %>%
  left_join(
    covars_dnr_compblocks_raw,
    by = "atlas_block"
  )


# 'Colonization-Absence' df
data_colabs_dnr <- spp_mod_data_dnr_full %>%
  filter(transition_state %in% c("Colonization", "Absence"))

data_colabs_dnr_z <- data_colabs_dnr %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))



# 'Extinction-Persistence' df
data_extper_dnr <- spp_mod_data_dnr_full %>%
  filter(transition_state %in% c("Extinction", "Persistence"))

data_extper_dnr_z <- data_extper_dnr %>%
  mutate(
    across(
      .cols = -c(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, det_Atlas1, det_Atlas2, common_name, alpha_code, ends_with("_z"))



# Assign binomial state identifiers for model
data_colabs_dnr_z <- data_colabs_dnr_z %>%
  mutate(
    col_abs = ifelse(transition_state %in% c("Colonization"), 1, 0)
  )

data_extper_dnr_z <- data_extper_dnr_z %>%
  mutate(
    ext_per = ifelse(transition_state %in% c("Extinction"), 1, 0)
  )



############ RED-BELLIED WOODPECKER ########

# Col-Abs
rbwo_mod_colabs_dnr <- glm(col_abs ~  pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_colabs_dnr_z, family = binomial)
summary(rbwo_mod_colabs_dnr)
autoplot(rbwo_mod_colabs_dnr)

# Ext-Per
rbwo_mod_extper_dnr <- glm(ext_per ~  pa_percent_z +
                         developed_lower_base_z + developed_upper_base_z + 
                         forest_deciduous_base_z + forest_mixed_base_z +
                         pasture_crop_base_z + 
                         forest_total_diff_z + developed_total_diff_z +
                         tmax_38yr_z + prcp_38yr_z,
                       data = data_extper_dnr_z, family = binomial)
summary(rbwo_mod_extper_dnr)
autoplot(rbwo_mod_extper_dnr)



#### VISUALIZE ####

rbwo_vars_plot <- c(
  "pa_percent_z",
  "developed_lower_base_z",
  "developed_upper_base_z",
  "forest_mixed_base_z",
  "forest_deciduous_base_z",
  "pasture_crop_base_z",
  "tmax_38yr_z",
  "prcp_38yr_z",
  "developed_total_diff_z",
  "forest_total_diff_z"
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



### PLOT ###

# COL-ABS
facet_colabs_dnr <- c(
  pa_percent_z = "Protected Area *",
  developed_lower_base_z = "Open + Lower Development **",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest .",
  forest_deciduous_base_z = "Deciduous Forest *",
  pasture_crop_base_z = "Pasture/Cropland",
  tmax_38yr_z = "Max Temp ***",
  prcp_38yr_z = "Precipitation *",
  developed_total_diff_z = "Difference in Total Developed Land",
  forest_total_diff_z = "Difference in Total Forest"
)

rbwo_predplot_colabs_dnr_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_colabs_dnr, 
                                   data = data_colabs_dnr_z)


ggplot(rbwo_predplot_colabs_dnr_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_colabs_dnr)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across DNR Comparable Blocks"
  )


# EXT-PER
facet_extper_dnr <- c(
  pa_percent_z = "Protected Area",
  developed_lower_base_z = "Open + Lower Development",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest",
  forest_deciduous_base_z = "Deciduous Forest",
  pasture_crop_base_z = "Pasture/Cropland",
  tmax_38yr_z = "Max Temp *",
  prcp_38yr_z = "Precipitation",
  developed_total_diff_z = "Difference in Total Developed Land",
  forest_total_diff_z = "Difference in Total Forest *"
)

rbwo_predplot_extper_dnr_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_extper_dnr, 
                                   data = data_extper_dnr_z)


ggplot(rbwo_predplot_extper_dnr_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_extper_dnr)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across DNR Comparable Blocks"
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



########################
### MODEL COMPARISON ###
########################

# Combine all model data
tab_colabs_rll  <- tidy(rbwo_mod_colabs_rll)  %>% mutate(dataset = "RLL",  model = "Colonization")
tab_extper_rll  <- tidy(rbwo_mod_extper_rll)  %>% mutate(dataset = "RLL",  model = "Extinction")

tab_colabs_dnr  <- tidy(rbwo_mod_colabs_dnr)  %>% mutate(dataset = "DNR",  model = "Colonization")
tab_extper_dnr  <- tidy(rbwo_mod_extper_dnr)  %>% mutate(dataset = "DNR",  model = "Extinction")

tab_all <- bind_rows(tab_colabs_rll, tab_extper_rll,
                     tab_colabs_dnr, tab_extper_dnr)

# Compute odds ratios, confidence intervals
tab_all <- tab_all %>%
  mutate(odds_ratio = exp(estimate),
         OR_low = exp(estimate - 1.96*std.error),
         OR_high = exp(estimate + 1.96*std.error))

# Comparison table
compare_table <- tab_all %>%
  dplyr::select(dataset, model, term, estimate, odds_ratio, p.value) %>%
  pivot_wider(names_from = dataset,
              values_from = c(estimate, odds_ratio, p.value),
              names_sep = "_")
compare_table


# Plot coefficient effect sizes
ggplot(tab_all, aes(x = term, y = estimate, color = dataset)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = estimate - 1.96*std.error,
                    ymax = estimate + 1.96*std.error),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  coord_flip() +
  facet_wrap(~model, scales = "free_y") +
  theme_bw() +
  labs(title = "Coefficient Comparison: RLL vs DNR",
       y = "Log-odds estimate")


# Model fit: AIC, deviance, McFadden pseudo-R^2
# Helper function
extract_metrics <- function(model, model_name) {
  data.frame(
    Model = model_name,
    AIC = AIC(model),
    Deviance = deviance(model),
    McFadden_R2 = pscl::pR2(model)[["McFadden"]]
  )
}

### ---- 1. COLONIZATION: DNR vs RLL ---- ###
col_results <- bind_rows(
  extract_metrics(rbwo_mod_colabs_dnr, "Colonization — DNR"),
  extract_metrics(rbwo_mod_colabs_rll, "Colonization — RLL")
)

### ---- 2. EXTINCTION: DNR vs RLL ---- ###
ext_results <- bind_rows(
  extract_metrics(rbwo_mod_extper_dnr, "Extinction — DNR"),
  extract_metrics(rbwo_mod_extper_rll, "Extinction — RLL")
)

### ---- Print results ---- ###
cat("=== Colonization Model Comparison ===\n")
print(col_results)

cat("\n\n=== Extinction Model Comparison ===\n")
print(ext_results)

# Export tidy table 
# Combine results
combined_results <- rbind(
  data.frame(Model = "Colonization — DNR", col_results[1, -1]),
  data.frame(Model = "Colonization — RLL", col_results[2, -1]),
  data.frame(Model = "Extinction — DNR", ext_results[1, -1]),
  data.frame(Model = "Extinction — RLL", ext_results[2, -1])
)

gt_table <- combined_results %>%
  gt() %>%
  tab_header(
    title = "Model Comparison: Colonization and Extinction",
    subtitle = "AIC, Deviance, and McFadden Pseudo-R²"
  ) %>%
  fmt_number(
    columns = c(AIC, Deviance, McFadden_R2),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 18,
    heading.subtitle.font.size = 14
  )

gt_table

