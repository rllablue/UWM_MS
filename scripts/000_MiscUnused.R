########################
### MISC CODE CHUNKS ###
########################

### --- ELEVATION DATA (NED) --- ###
# Mean elevation per block done in ArcPro (1/3 arc-sec DEMs,~10m)

elevation_summary <- read_csv("maps/ned/ned_wibba_table.csv") # read in data calc'd in ArcPro

elevation_summary <- elevation_summary %>%
  rename(
    atlas_block = BLOCK_ID, #rename variables to match
    elev_mean = MEAN
  ) %>%
  mutate(atlas_block = as.character(atlas_block)) %>%
  mutate(
    elev_mean_z = as.numeric(scale(elev_mean)) # z-standard
  )

wibba_summary <- wibba_summary %>%
  left_join(elevation_summary %>% dplyr::select(atlas_block, elev_mean_z), by = "atlas_block") # join to summary df



### --- GEOGRAPHIC FILTER --- ###

### Locates checklists within Atlas block grid by lat, lon coordinates and 
# assigns missing block IDs to checklists, observations. 
# Uneeded for WIBBA.


# geographic filtering, clean-up
wibba_grid <- read_sf("maps/wibba blocks/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # read in atlas grid shape file
  st_transform(4269) # convert to NAD83, crs of shp file

checklists_sf <- checklists %>%
  select(checklist_id, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% # convert checklist locations to points
  st_transform(4269) # match atlas block crs

checklists_joined <- st_join(checklists_sf, wibba_grid[, c("BLOCK_ID")], left = TRUE) # create df with checklist id, lat/lon, and BLOCK_ID

in_region <- checklists_joined %>%
  filter(!is.na(BLOCK_ID)) # keep only checklists in joined atlas block shp

checklists_wip <- checklists %>%
  semi_join(st_drop_geometry(in_region), by = "checklist_id") %>% # remove checklists outside atlas blocks
  left_join(st_drop_geometry(in_region[, c("checklist_id", "BLOCK_ID")]), by = "checklist_id") %>% # adds BLOCK_ID from shp to checklist data w/o overwriting
  mutate(
    atlas_block = if_else(is.na(atlas_block), BLOCK_ID, atlas_block) # fills in missing atlas blocks w/o overwriting
  ) %>%
  select(-BLOCK_ID) #remove BLOCK_ID column

observations_wip <- observations %>%
  semi_join(checklists, by = "checklist_id") %>%  # include only matching checklist_id from checklists df
  left_join(
    checklists %>% 
      select(checklist_id, atlas_block) %>%  # join atlas_block to observations using checklist_id
      rename(joined_block = atlas_block),
    by = "checklist_id"
  ) %>%
  mutate(
    old_block = atlas_block, # preserve existing
    atlas_block = joined_block # overwrite with cleaned value
  ) %>%
  select(-joined_block) # drop temp join column



### --- OLD MODELING APPROACHES --- ###

# GROUPED DATA # 
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


# Binomial state IDs
data_colper_z <- data_colper_z %>%
  mutate(
    col_per = ifelse(transition_state %in% c("Colonization"), 1, 0),
  )

data_extabs_z <- data_extabs_z %>%
  mutate(
    ext_abs = ifelse(transition_state %in% c("Extinction"), 1, 0)
  )


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


# Visualize

# COL-PER
facet_colper <- c(
  pa_percent_z = "Protected Area",
  developed_lower_base_z = "Open + Lower Development",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest ***",
  forest_deciduous_base_z = "Deciduous Forest .",
  pasture_crop_base_z = "Pasture/Cropland",
  tmax_38yr_z = "Max Temp ***",
  prcp_38yr_z = "Precipitation",
  developed_total_diff_z = "Difference in Total Developed Land",
  forest_total_diff_z = "Difference in Total Forest"
)

rbwo_predplot_colper_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                   model = rbwo_mod_colper, 
                                   data = data_colper_z)


ggplot(rbwo_predplot_colper_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_colper)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Colonization",
    title = "Marginal Effects of Covariates on Colonization Probability Across Blocks"
  )


# EXT-ABS
facet_extabs <- c(
  pa_percent_z = "Protected Area",
  developed_lower_base_z = "Open + Lower Development",
  developed_upper_base_z = "Moderate + High Development",
  forest_mixed_base_z = "Mixed Forest **",
  forest_deciduous_base_z = "Deciduous Forest",
  pasture_crop_base_z = "Pasture/Cropland",
  tmax_38yr_z = "Max Temp ***",
  prcp_38yr_z = "Precipitation",
  developed_total_diff_z = "Difference in Total Developed Land",
  forest_total_diff_z = "Difference in Total Forest"
)

rbwo_predplot_ext_df <- map_dfr(rbwo_vars_plot, make_pred_df, 
                                model = rbwo_mod_ext, 
                                data = spp_mod_data_ext_z)


ggplot(rbwo_predplot_ext_df, aes(x = x_value, y = pred_prob)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, 
             scales = "free_x",
             labeller = labeller(variable = facet_extabs)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Covariate Value (z-score)",
    y = "Predicted Probability of Extinction",
    title = "Marginal Effects of Covariates on Extinction Probability Across Blocks"
  )



### --- SR RESIDUAL REGRESSION --- ###

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
wibba_covars_raw <- wibba_covars_raw %>%
  left_join(
    sr_summary_comp %>%
      dplyr::select(atlas_block, sr_Atlas1, sr_Atlas2),
    by = "atlas_block"
  )


# Extract relevant covars for residual regression, z-scale ################### Note: Refine this list? Also how to scale? (global vs response-subset block data)
covars_for_sr <- wibba_covars_raw %>%
  dplyr::select(
    atlas_block,
    pa_percent, water_open_base, shrub_scrub_base, grassland_base,
    developed_total_base, forest_total_base, wetlands_total_base,
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

covars_for_sr <- covars_for_sr %>%
  left_join(
    sr_summary_comp %>%
      dplyr::select(atlas_block, sr_Atlas1, sr_Atlas2),
    by = "atlas_block"
  ) %>%
  dplyr::select(-c(2:10))

# Fit models on un-standardized SR; logistic mod
# Use grouped land cover covariates, climate proxy covariates for ease
mod_srA1 <- lm(sr_Atlas1 ~ pa_percent_z + water_open_base_z + shrub_scrub_base_z + grassland_base_z + developed_total_base_z + forest_total_base_z + wetlands_total_base_z + lat_z + lon_z, data = covars_for_sr)
mod_srA2 <- lm(sr_Atlas2 ~ pa_percent_z + water_open_base_z + shrub_scrub_base_z + grassland_base_z + developed_total_base_z + forest_total_base_z + wetlands_total_base_z + lat_z + lon_z, data = covars_for_sr)

# Extract residuals, add to summaries
sr_resids <- covars_for_sr %>%
  mutate(
    srA1_resid = resid(mod_srA1),
    srA2_resid = resid(mod_srA2)
  ) %>%
  dplyr::select(atlas_block, srA1_resid, srA2_resid)

sr_summary_comp <- sr_summary_comp %>%
  left_join(sr_resids, by = "atlas_block")

wibba_covars_raw <- wibba_covars_raw %>%
  left_join(sr_resids, by = "atlas_block") %>%
  dplyr::select(-sr_Atlas1, -sr_Atlas2)



### --- SIMPLE MODELS --- ### 
# Models by response, block sets w/o model and covariate selection applied

######################
### BASIC MODELING ###
######################

# Modeling w/o covariate thinning, selection

# Species = RBWO

### --- RLL SUBSET --- ###

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

# Species = RBWO

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