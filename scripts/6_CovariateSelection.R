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
library(MuMIn)

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

### --- DATA WRANGLING --- ###

# Prepare df with data for species-specific modeling of covariates
wibba_mod_covars_z <- wibba_modeling_covars %>% # z-scaled covars df
  mutate(
    across(
      .cols = -c(atlas_block, transition_state),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, ends_with("_z"))

## FUNCTION ##
# --- Candidate selection function ---
explore_covariates <- function(df, response_var, covars) {
  
  results <- list()
  
  # Fit model for each single covariate
  for(cov in covars) {
    fmla <- as.formula(paste(response_var, "~", cov))
    mod <- glm(fmla, data = df, family = binomial())
    results[[cov]] <- broom::tidy(mod) %>%
      mutate(term = cov, AIC = AIC(mod))
  }
  
  # Fit model for all covariates together
  if(length(covars) > 1) {
    fmla_all <- as.formula(paste(response_var, "~", paste(covars, collapse = " + ")))
    mod_all <- glm(fmla_all, data = df, family = binomial())
    results[["all_covars"]] <- broom::tidy(mod_all) %>%
      mutate(term = "all_covars", AIC = AIC(mod_all))
  }
  
  results_df <- bind_rows(results)
  return(results_df)
}


## APPLY ##
spp_name <- "Red-bellied Woodpecker"

spp_mod_data <- breeders_zf_summary %>%
  filter(common_name == spp_name) %>%
  left_join(wibba_mod_covars_z %>% 
    dplyr::select(-transition_state), by = "atlas_block") %>%
  mutate(
    col = ifelse(det_Atlas1 == 0 & det_Atlas2 == 1, 1, 0),
    per = ifelse(det_Atlas1 == 1 & det_Atlas2 == 1, 1, 0),
    ext = ifelse(det_Atlas1 == 1 & det_Atlas2 == 0, 1, 0)
  )

# Separate out land, climate covars
land_covars <- grep("_2008_z$|_diff_z$", names(spp_mod_data), value = TRUE)
clim_covars <- grep("_38yr_z$|_diffyrs_z$", names(spp_mod_data), value = TRUE)

responses <- c("col", "per", "ext")


# Explore covariates
candidate_results <- map_dfr(responses, function(resp) {
  map_dfr(list(land = land_covars, climate = clim_covars), function(cset) {
    explore_covariates(spp_mod_data, resp, cset) %>%
      mutate(response = resp)
  }, .id = "covar_set")  # .id automatically creates "covar_set" column with list names
})

# Summarize best covariates by AIC
candidate_summary <- candidate_results %>%
  group_by(response, covar_set) %>%
  filter(AIC == min(AIC)) %>%
  select(response, covar_set, term, AIC) %>%
  arrange(response, covar_set, AIC)

candidate_summary








### --- DIAGNOSTICS --- ###

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















































############################# THE WEEDS ########################################

### Threshold filtering and penalized regression techniques to fix
# covariate sets for species/guild-specific modeling (ie. restrict 
# covariates to those relevant to species-specific block pool)


### --- FILTERING --- ##

Prepare_covars <- function(
    species_name,
    zf_summary = breeders_zf_summary,
    covars_df = wibba_modeling_covars,
    base_lc_threshold = 1,
    diff_lc_threshold = 2,
    diff_clim_threshold = 0.01,
    quantile_cut = 0.10
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
  diff_clim_cols <- grep("_diffyrs$", names(species_mod_data_raw), value = TRUE)
  
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
    ) %>%
    # Keep only z-scaled columns + identifiers
    dplyr::select(atlas_block, transition_state, ends_with("_z"))
  
  return(species_mod_data_z)

}


### --- APPLY --- ###

bobo_covars <- Prepare_covars("Bobolink")
head(bobo_covars)

rbwo_covars <- Prepare_covars("Red-bellied Woodpecker")
head(rbwo_covars)

colo_covars <- Prepare_covars("Common Loon")
head(colo_covars)



################################################################################
################NEW

Diagnose_covars <- function(species_mod_data, alpha = 0.5, lambda_type = "lambda.min", verbose = TRUE) {
  
  # Keep only z-scaled covariates + identifiers
  numeric_covars <- species_mod_data %>%
    dplyr::select(ends_with("_z")) %>%
    na.omit()
  
  # Ensure transition_state is aligned
  transition_vec <- species_mod_data$transition_state
  if(length(transition_vec) != nrow(numeric_covars)) {
    stop("Number of rows in transition_state does not match numeric covariates (after removing NAs)")
  }
  
  # Combine for VIF
  vif_df <- cbind(
    transition_numeric = as.numeric(factor(transition_vec)),
    numeric_covars
  )
  
  # Pairwise correlations
  corr_matrix <- cor(numeric_covars)
  if(verbose) {
    corrplot::corrplot(corr_matrix, method = "color", type = "upper", tl.cex = 0.6)
  }
  
  # VIF
  vif_test <- lm(transition_numeric ~ ., data = vif_df)
  vif_vals <- car::vif(vif_test)
  if(verbose) {
    cat("VIF values:\n")
    print(vif_vals)
  }
  
  # Elastic Net design matrices
  x_multi <- model.matrix(transition_state ~ . - 1, data = cbind(transition_state = transition_vec, numeric_covars))
  y_multi <- transition_vec
  
  x_logi <- x_multi
  y_logi <- ifelse(transition_vec == "Presence", 1, 0)
  
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
  
  # Return everything
  return(list(
    correlation_matrix = corr_matrix,
    vif = vif_vals,
    enet_multinomial = enet_multi,
    enet_logistic = enet_logi,
    selected_multinomial = selected_multi,
    selected_logistic = selected_logi
  ))
}

# --- APPLY ---
bobo_diagnostics <- Diagnose_covars(bobo_covars)

rbwo_diagnostics <- Diagnose_covars(rbwo_covars)

colo_diagnostics <- Diagnose_covars(colo_covars)





############################################# OLD 
### --- DIAGNOSTICS, SELECTION --- ###

### Function to diagnose collinearity among, select best covariates per species
# to model (pairwise correlations, VIF, Elastic Net penalized regression).
# Established value |r| > 0.7 = highly correlated
# Established value VIF > 3-5 = moderate-strong collinearity

Diagnose_covars <- function(species_mod_data, alpha = 0.5, lambda_type = "lambda.min", verbose = TRUE) {
  
  # Prepare relevant covariates
  numeric_covars <- species_mod_data %>%
    dplyr::select(-atlas_block, -transition_state) %>%
    dplyr::select(ends_with("_z"))  %>%
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

bobo_diagnostics <- Diagnose_covars(bobo_covars)