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
library(readxl)

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
library(lme4)
library(MuMIn)
library(pscl)
library(AICcmodavg)

# Visualization
library(ggplot2)
library(ggfortify)
library(viridis)
library(gt)


### --- FLEXIBLE SPECS --- ###

spp_name <- "Red-bellied Woodpecker"


### --- DATAFRAMES --- ###

# Carry-over DataFrames #
spp_zf_rll <- read.csv("outputs/data/breeders_zf_summary.csv") 
covars_raw_rll <- read.csv("outputs/data/covars_raw_all.csv")

wibba_summary_rll <- read.csv("data/summaries/wibba_summary_rll.csv") # df
blocks_rll <- wibba_summary_rll$atlas_block # vector

blocks_dnr <- read_xlsx("data/summaries/CompBlocks_DNR2023.xlsx") # df
blocks_dnr <- blocks_dnr$atlas_block # vector


# SPECIES RICHNESS / EFFORT PROXY (WIP, RE-ORDER)
covars_raw_rll <- covars_raw_rll %>%
  mutate(sr_Diff = sr_Atlas2 - sr_Atlas1)


# Covariate Sets #
factor_covars_all <- c("atlas_block", "common_name", "alpha_code", "transition_state")

stable_covars_all <- c("lon", "lat", "sr_Diff", "pa_percent")

land_covars_all <- c("water_open_base", "barren_land_base", "shrub_scrub_base", "grassland_base",
                     "developed_open_base", "developed_low_base", "developed_med_base", "developed_high_base", 
                     "developed_lower_base", "developed_upper_base", "developed_total_base", 
                     "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                     "forest_total_base", "pasture_base", "cropland_base", "pasture_crop_base", 
                     "wetlands_woody_base", "wetlands_herb_base", "wetlands_total_base",
                     
                     "water_open_diff", "barren_land_diff", "shrub_scrub_diff", "grassland_diff",
                     "developed_open_diff", "developed_low_diff", "developed_med_diff", "developed_high_diff", 
                     "developed_lower_diff", "developed_upper_diff", "developed_total_diff", 
                     "forest_deciduous_diff", "forest_evergreen_diff", "forest_mixed_diff", 
                     "forest_total_diff", "pasture_diff", "cropland_diff", "pasture_crop_diff", 
                     "wetlands_woody_diff", "wetlands_herb_diff", "wetlands_total_diff")

land_covars_base <- c("water_open_base", "barren_land_base", "shrub_scrub_base", "grassland_base",
                      "developed_open_base", "developed_low_base", "developed_med_base", "developed_high_base", 
                      "developed_lower_base", "developed_upper_base", "developed_total_base", 
                      "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                      "forest_total_base", "pasture_base", "cropland_base", "pasture_crop_base", 
                      "wetlands_woody_base", "wetlands_herb_base", "wetlands_total_base")

land_covars_diff <- c("water_open_diff", "barren_land_diff", "shrub_scrub_diff", "grassland_diff",
                      "developed_open_diff", "developed_low_diff", "developed_med_diff", "developed_high_diff", 
                      "developed_lower_diff", "developed_upper_diff", "developed_total_diff", 
                      "forest_deciduous_diff", "forest_evergreen_diff", "forest_mixed_diff", 
                      "forest_total_diff", "pasture_diff", "cropland_diff", "pasture_crop_diff", 
                      "wetlands_woody_diff", "wetlands_herb_diff", "wetlands_total_diff")

climate_covars_all <- c("tmax_38yr", "tmin_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff")

climate_covars_base <- c("tmax_38yr", "tmin_38yr", "prcp_38yr")

climate_covars_diff <- c("tmax_diff", "tmin_diff", "prcp_diff")



# New DataFrames #
# RLL modeling df
mod_data_rll <- spp_zf_rll %>%
  filter(common_name == spp_name) %>%
  left_join(covars_raw_rll, by = "atlas_block")

write.csv(mod_data_rll, "outputs/data/mod_data_rll.csv", row.names = FALSE)

# DNR modeling df
mod_data_dnr <- spp_zf_rll %>%
  filter(atlas_block %in% blocks_dnr, 
         common_name == spp_name) %>%  
  left_join(covars_raw_rll, by = "atlas_block")

write.csv(mod_data_dnr, "outputs/data/mod_data_dnr.csv", row.names = FALSE)



### -- COVARIATE THINNING --- ###

# PRE-SCREENING #
### Screen for biologically relevant covariates (for landcover and climate)
# on a species-, state-specific basis, ie. within-group reduction; run pairwise 
# correlations, VIF to assess multi-collinearity (alt. approach: PCA &/or CLT)
# and further thin predictors

# Full covariate sets
factor_covars_all
stable_covars_all
land_covars_all
  land_covars_base
  land_covars_diff
climate_covars_all
  climate_covars_base
  climate_covars_diff


# Species-specific Thinned Covariate Sets
spp_name <- "Red-bellied Woodpecker"


factor_covars_reduced <- c("atlas_block", "transition_state")

stable_covars_all <- c("lon", "lat", "sr_Diff", "pa_percent")

land_covars_reduced <- c("shrub_scrub_base", 
                         "grassland_base", "developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "pasture_crop_base", "wetlands_total_base",
                         
                         "developed_total_diff", "forest_total_diff", "wetlands_total_diff")

covars_numeric_reduced <- c(stable_covars_all, land_covars_reduced, climate_covars_all)


# State-specific Covariate Sets
### Separate data into bins where A2 detection is either 1 or 0; not necessarily
# in precise probabilities between colo, pers, abs, ext, but more what promotes
# 'new' v. 'continued' colonization, ie. what promotes det = 1, as opposed to 
# what promotes det = 0. 
### Scaling of covars w/in subsets for relevant normalized values

spp_name <- "Red-bellied Woodpecker"

# RLL
mod_colabs_rll_z <- mod_data_rll %>%
  filter(transition_state %in% c("Colonization", "Absence")) %>% # filter by state
  dplyr::select(all_of(factor_covars_reduced), # select numeric, factor covars
                all_of(covars_numeric_reduced)) %>%
  mutate(across( # scale only numeric covars
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(col_abs = ifelse(transition_state == "Colonization", 1, 0)) %>% # binomial response variable
  dplyr::select(all_of(factor_covars_reduced), # columns to keep
         ends_with("_z"),
         col_abs)

mod_extper_rll_z <- mod_data_rll %>%
  filter(transition_state %in% c("Extinction", "Persistence")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(ext_per = ifelse(transition_state == "Extinction", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                ext_per)

# DNR
mod_colabs_dnr_z <- mod_data_dnr %>%
  filter(transition_state %in% c("Colonization", "Absence")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(col_abs = ifelse(transition_state == "Colonization", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                col_abs)

mod_extper_dnr_z <- mod_data_dnr %>%
  filter(transition_state %in% c("Extinction", "Persistence")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(ext_per = ifelse(transition_state == "Extinction", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                ext_per)


# CORRELATIONS, COLLINEARITY # 

# Pairwise Correlations
covar_cols_colabs_rll <- grep("_z$", names(mod_colabs_rll_z), value = TRUE)
covar_cols_extper_rll <- grep("_z$", names(mod_extper_rll_z), value = TRUE)

covar_cols_colabs_dnr <- grep("_z$", names(mod_colabs_dnr_z), value = TRUE)
covar_cols_extper_dnr <- grep("_z$", names(mod_extper_dnr_z), value = TRUE)


M1 <- cor(mod_colabs_rll_z[, covar_cols_colabs_rll], use = "pairwise.complete.obs")
M2 <- cor(mod_extper_rll_z[, covar_cols_extper_rll], use = "pairwise.complete.obs")

M3 <- cor(mod_colabs_dnr_z[, covar_cols_colabs_dnr], use = "pairwise.complete.obs")
M4 <- cor(mod_extper_dnr_z[, covar_cols_extper_dnr], use = "pairwise.complete.obs")


corrplot(M1, method = "color", tl.cex = 0.7, number.cex = 0.6)
corrplot(M2, method = "color", tl.cex = 0.7, number.cex = 0.6)

corrplot(M3, method = "color", tl.cex = 0.7, number.cex = 0.6)
corrplot(M4, method = "color", tl.cex = 0.7, number.cex = 0.6)


high_corr1 <- which(abs(M1) > 0.7 & abs(M1) < 1, arr.ind = TRUE)
apply(high_corr1, 1, function(i) cat(rownames(M1)[i[1]], "-", colnames(M1)[i[2]], "r =", M1[i[1],i[2]], "\n"))

high_corr2 <- which(abs(M2) > 0.7 & abs(M2) < 1, arr.ind = TRUE)
apply(high_corr2, 1, function(i) cat(rownames(M2)[i[1]], "-", colnames(M2)[i[2]], "r =", M2[i[1],i[2]], "\n"))

high_corr3 <- which(abs(M3) > 0.7 & abs(M3) < 1, arr.ind = TRUE)
apply(high_corr3, 1, function(i) cat(rownames(M3)[i[1]], "-", colnames(M3)[i[2]], "r =", M3[i[1],i[2]], "\n"))

high_corr4 <- which(abs(M4) > 0.7 & abs(M4) < 1, arr.ind = TRUE)
apply(high_corr4, 1, function(i) cat(rownames(M4)[i[1]], "-", colnames(M4)[i[2]], "r =", M4[i[1],i[2]], "\n"))


# Correlation Thinned Covariate Sets
factor_covars_reduced

land_covars_reduced <- c("shrub_scrub_base", "pasture_crop_base",
                         "grassland_base","developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_total_base","forest_total_diff")

climate_covars_reduced <- c("tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff")

stable_covars_reduced <- c("sr_Diff", "pa_percent")

covars_numeric_reduced <- c(land_covars_reduced, climate_covars_reduced, stable_covars_reduced)
covars_numeric_reduced_z <- paste0(covars_numeric_reduced, "_z")



# VIF

vif_model1 <- glm(col_abs ~ ., data = mod_colabs_rll_z[, c("col_abs", covars_numeric_reduced_z)], family = binomial)
vif(vif_model1)
alias(vif_model1)

vif_model2 <- glm(ext_per ~ ., data = mod_extper_rll_z[, c("ext_per", covars_numeric_reduced_z)], family = binomial)
vif(vif_model2)
alias(vif_model2)

vif_model3 <- glm(col_abs ~ ., data = mod_colabs_dnr_z[, c("col_abs", covars_numeric_reduced_z)], family = binomial)
vif(vif_model3)
alias(vif_model3)

vif_model4 <- glm(ext_per ~ ., data = mod_extper_dnr_z[, c("ext_per", covars_numeric_reduced_z)], family = binomial)
vif(vif_model4)
alias(vif_model4)

# VIF Thinned Covariate Sets
# Final set prior to AICc/Model Selection process

factor_covars_final <- c("atlas_block")

land_covars_final <- c("shrub_scrub_base", "grassland_base","developed_total_base",
                       "forest_deciduous_base", "forest_evergreen_base", 
                       "forest_mixed_base", "wetlands_total_base","forest_total_diff")

climate_covars_final <- climate_covars_reduced # "tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff"

stable_covars_final <- stable_covars_reduced # "sr_Diff", "pa_percent"


covars_numeric_final <- c("land_covars_final", "climate_covars_final", "stable_covars_final")
covars_numeric_final_z <- paste0(covars_numeric_final, "_z")

covars_all_final <- c(covars_numeric_reduced_z, factor_covars_final)
# w/in each modeling df: col_abs or ext_per



#######################
### MODEL SELECTION ### 
#######################

### --- PROCESS --- ###

### AICc from pre-screened, correlation-controlled, biologically plausible 
# candidates with both additive and interaction terms; include atlas_block as  
# a random effect, as both block sets are sub-samples of full Atlas block set.    ### Does this make sense? 
                                                                                  # Are those that were chosen/completed
                                                                                  # to a satisfactory degree "random" enough?

# MuMIn::dredge() automates candidate set construction of all main effect covar combos
# but user must manually define, limit interaction terms. User should also further
# restrict candidate models prior to running dredge(), otherwise models will not 
# converge and/or overload processor. 
# Also: add atlas_block as a random effect.

### Two directions and/or steps for candidate model thinning:
# 1) Partitioned model sets, ie. bin covar groups ('land', 'climate', etc.), run model subsets
# 2) Targeted, manually constructed 2-way interactions, ie. choose which interactions to run

# Covar set: covars_all_final


### APPROACH 1: Partition Covariate Models ###


# DataFrames # 
mod_colabs_rll_aicc <- mod_colabs_rll_z[, c("col_abs", covars_all_final)]
mod_extper_rll_aicc <- mod_extper_rll_z[, c("ext_per", covars_all_final)]

mod_colabs_dnr_aicc <- mod_colabs_dnr_z[, c("col_abs", covars_all_final)]
mod_extper_dnr_aicc <- mod_extper_dnr_z[, c("ext_per", covars_all_final)]


# Models #
# Model structure
### Step needed to create, set limit for interactions

mod_form_colabs <- as.formula(paste("col_abs ~ (", paste(covars_all_final, collapse = " + "), ")^2")) # ^2 = max 2-way interactions


mod_form_extper <- as.formula(paste("ext_per ~ (", paste(covars_all_final, collapse = " + "), ")^2"))

# Model data
# RLL
full_mod_colabs_rll <- glm(mod_form_colabs, data = mod_colabs_rll_aicc, family = binomial)
full_mod_extper_rll <- glm(mod_form_extper, data = mod_extper_rll_aicc, family = binomial)

# DNR
full_mod_colabs_dnr <- glm(mod_form_colabs, data = mod_colabs_dnr_aicc, family = binomial)
full_mod_extper_dnr <- glm(mod_form_extper, data = mod_extper_dnr_aicc, family = binomial)


# Model Selection # 
### Use MuMIn::dredge() to automate model selection process across all possible subsets
# w/o having to manually code/write all formulas of interest. Created model form code
# above b/c dredge() by default only tests main (additive) effects, ie. need to add manually.

options(na.action = "na.fail") # required for MuMIn::dredge()

# RLL
dredge_colabs_rll <- dredge(full_model_colabs_rll, rank = "AICc")
dredge_extper_rll <- dredge(full_model_extper_rll, rank = "AICc")

# DNR
dredge_colabs_dnr <- dredge(full_model_colabs_dnr, rank = "AICc")
dredge_extper_dnr <- dredge(full_model_extper_dnr, rank = "AICc")


### APPROACH B: MANUAL APPROACH

















################
### MODELING ###
################

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

