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

# Model and covariate fitting, evaluation
library(nnet)
library(stats)
library(broom)
library(car)
library(performance)
library(DescTools)
library(Metrics)
library(pROC)
library(rsample)
library(glmnet)
library(forcats)

# Spatial
library(sf)
library(units)

# Visualization
library(ggplot2)
library(viridis)
library(RColorBrewer)



###########################
### COVARIATE SELECTION ###
###########################

### Threshold filtering and penalized regression techniques to fix
# covariate sets for species/guild-specific modeling (ie. restrict 
# covariates to those relevant to species-specific block pool)


### --- FILTERING --- ##

Prepare_covars <- function(
    species_name,
    zf_summary = breeders_zf_summary,
    covars_df = wibba_modeling_covars,
    base_lc_threshold = 10,
    diff_lc_threshold = 2,
    diff_clim_threshold = 0.01,
    quantile_cut = 0.05
) {
  
  # Subset species detection data
  species_data <- zf_summary %>%
    filter(common_name == species_name)
  
  # Join covariates
  species_mod_data_raw <- species_data %>%
    left_join(covars_df, by = "atlas_block")
  
  
  # Identify landcover columns
  base_lc_cols <- grep("_2008$", names(species_mod_data_raw), value = TRUE)
  diff_lc_cols <- grep("_diff$", names(species_mod_data_raw), value = TRUE)
  
  # Filter landcover
  keep_base_lc <- base_lc_cols[
    sapply(species_mod_data_raw[base_lc_cols], function(x) quantile(x, quantile_cut, na.rm = TRUE) >= base_lc_threshold)
  ]
  
  keep_diff_lc <- diff_lc_cols[
    sapply(species_mod_data_raw[diff_lc_cols], function(x) quantile(abs(x), quantile_cut, na.rm = TRUE) >= diff_lc_threshold)
  ]
  
  
  # Identify climate columns
  base_clim_cols <- grep("_38yr$", names(species_mod_data_raw), value = TRUE)
  diff_clim_cols <- grep("_diff$", names(species_mod_data_raw), value = TRUE)
  
  keep_base_clim <- base_clim_cols
  keep_diff_clim <- diff_clim_cols[
    sapply(species_mod_data_raw[diff_clim_cols], function(x) quantile(abs(x), quantile_cut, na.rm = TRUE) >= diff_clim_threshold)
  ]
  
  # Non-filtered covars
  nonfilt_covars <- c(
    "pa_percent",
    "lat",
    "lon",
    grep("^sr", names(species_mod_data_raw), value = TRUE)
  )
  
  # Combine all covars to keep
  species_covars_raw <- c(
    keep_base_lc,
    keep_diff_lc,
    keep_base_clim,
    keep_diff_clim,
    nonfilt_covars
  )
  
  # Filter, z-scale covars
  species_mod_data_z <- species_mod_data_raw %>%
    dplyr::select(atlas_block, transition_state, all_of(species_covars_raw)) %>%
    mutate(
      across(
        .cols = -c(atlas_block, transition_state),
        .fns = ~ as.numeric(scale(.)),
        .names = "{.col}_z"
      )
    )
  
  # Return both raw and z-scaled data
  return(list(
    raw = species_mod_data_raw %>% dplyr::select(atlas_block, transition_state, all_of(species_covars_raw)),
    z_scaled = species_mod_data_z
  ))
}

### --- APPLY --- ###

bobo_covars <- Prepare_covars("Bobolink")
head(bobo_covars$raw)
head(bobo_covars$z_scaled)






################################################################################
# Set up flexible species filter
species_name <- "Bobolink"

species_data <- breeders_zf_summary %>%
  filter(common_name == species_name)

species_mod_data_raw <- species_data %>%
  left_join(wibba_modeling_coavrs, by = "atlas_block")

# Species-specific landcover filter
base_lc_cols <- grep("_2008$", names(species_mod_data_raw), value = TRUE) # single ref year
diff_lc_cols <- grep("_diff$", names(species_mod_data_raw), value = TRUE) # difference between ref years

base_lc_threshold <- 10 # min % cover
diff_lc_threshold <- 2  # min % change
quantile_cut <- 0.05  # quantile

keep_base_lc <- base_lc_cols[
  sapply(species_mod_data_raw[base_lc_cols], function(x) quantile(x, quantile_cut, na.rm = TRUE) >= base_lc_threshold)
]

keep_diff_lc <- diff_lc_cols[
  sapply(species_mod_data_raw[diff_lc_cols], function(x) quantile(abs(x), quantile_cut, na.rm = TRUE) >= diff_lc_threshold)
]

# Species-specific climate change filter
diff_clim_threshold <- 0.01

keep_diff_clim <- diff_clim_cols[
  sapply(species_mod_data_raw[diff_clim_cols], function(x) quantile(abs(x), quantile_cut, na.rm = TRUE) >= diff_clim_threshold)
]

# Non-filtered covars
keep_base_clim <- names(species_mod_data_raw) %>% # keep all base climate metric (3) values
  grep("_38yr$", ., value = TRUE)

nonfilt_covars <- c(
  "pa_percent", 
  "lat",
  "lon",
  grep("^sr", names(species_mod_data_raw), value = TRUE)
)


# Species-specific raw covariates, df for modeling 
species_covars_raw <- c(
  keep_base_lc,
  keep_diff_lc,
  keep_base_clim,
  keep_diff_clim,
  nonfilt_covars
)

species_mod_data_z <- species_mod_data_raw %>%
  dplyr::select(atlas_block, transition_state, all_of(species_covars_raw)) %>% 
  mutate(
    across(
      .cols = -c(atlas_block, transition_state),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  )

################################################################################



### --- DIAGNOSTICS, SELECTION --- ###

### Function to diagnose collinearity among, select best covariates per species
# to model (pairwise correlations, VIF, Elastic Net penalized regression).
# Established value |r| > 0.7 = highly correlated
# Established value VIF > 3-5 = moderate-strong collinearity

Diagnose_covars <- function(species_mod_data, alpha = 0.5, lambda_type = "lambda.min", verbose = TRUE) {
  
  # Prepare relevant covariates
  numeric_covars <- species_mod_data %>%
    dplyr::select(-atlas_block, -transition_state) %>%
    na.omit()
  
  
  # Pairwise correlations
  corr_matrix <- cor(numeric_covars)
  if(verbose) {
    corrplot::corrplot(corr_matrix, method = "color", type = "upper", tl.cex = 0.6)
  }
  

  # VIF
  vif_test <- lm(as.numeric(factor(species_mod_data$transition_state)) ~ ., data = numeric_covars)
  vif_vals <- car::vif(vif_test)
  if(verbose) {
    print("VIF values:")
    print(vif_vals)
  }
  

  # Elastic Net
  # Design matrices
  x_multi <- model.matrix(transition_state ~ . - 1, data = species_mod_data %>% dplyr::select(-atlas_block))
  y_multi <- species_mod_data$transition_state
  
  x_logi <- x_multi
  y_logi <- ifelse(species_mod_data$transition_state == "Presence", 1, 0)
  
  set.seed(123)
  
  # Multinomial elastic net
  enet_multi <- cv.glmnet(x = x_multi, y = y_multi, family = "multinomial", alpha = alpha, type.measure = "class")
  coef_multi <- coef(enet_multi, s = lambda_type)
  
  # Logistic elastic net
  enet_logi <- cv.glmnet(x = x_logi, y = y_logi, family = "binomial", alpha = alpha, type.measure = "class")
  coef_logi <- coef(enet_logi, s = lambda_type)
  
  
  # Extract non-zero predictors
  selected_multi <- lapply(coef_multi, function(x) rownames(x)[x[,1] != 0])
  selected_logi <- rownames(coef_logi)[coef_logi[,1] != 0]
  
  if(verbose) {
    cat("\nSelected covariates (Multinomial):\n")
    print(selected_multi)
    cat("\nSelected covariates (Logistic):\n")
    print(selected_logi)
  }
  
  
  # Return results
  return(list(
    correlation_matrix = corr_matrix,
    vif = vif_vals,
    enet_multinomial = enet_multi,
    enet_logistic = enet_logi,
    selected_multinomial = selected_multi,
    selected_logistic = selected_logi
  ))
}


### --- APPLY --- ###

bobo_diagnostics <- Diagnose_covars(bobo_covars$z_scaled)
















################################################################################


# PAIRWISE CORRELATIONS #
# Established value |r| > 0.7 = highly correlated
numeric_covars <- species_mod_data %>%
  dplyr::select(ends_with("_z")) %>%
  na.omit()

corr_matrix <- cor(numeric_covars)
corrplot::corrplot(corr_matrix, method = "color", type = "upper", tl.cex = 0.6)


# VIF #
# Established value VIF > 3-5 = moderate-strong collinearity
vif_test <- lm(transition_state ~ ., data = numeric_covars)
car::vif(vif_test)


# ELASTIC NET #
# Penalized regression to identify more informative predictors

# Design matrix
x <- model.matrix(transition_state ~ . - 1, data = species_mod_data %>% select(ends_with("_z")))
y <- species_mod_data$transition_state

# Fit net 
# (lasso, ridge penalty balance alpha = 0.5)