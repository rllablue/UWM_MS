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
library(corrplot)
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
library(reshape2)
library(RColorBrewer)



###########################
### COVARIATE SELECTION ###
###########################

### --- DATA WRANGLING --- ###

# DATAFRAMES #
# Z-scaled nmodeling covars
wibba_mod_covars_z <- wibba_modeling_covars %>% # z-scaled covars df
  mutate(
    across(
      .cols = -c(atlas_block, transition_state),
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  dplyr::select(atlas_block, transition_state, ends_with("_z"))


# Species-specific covar, detection df
spp_name <- "Red-bellied Woodpecker"

spp_mod_data <- breeders_zf_summary %>%
  filter(common_name == spp_name) %>%
  left_join(
    wibba_mod_covars_z %>% dplyr::select(-transition_state),
    by = "atlas_block"
  ) %>%
  mutate(
    col = ifelse(det_Atlas1 == 0 & det_Atlas2 == 1, 1, 0),
    per = ifelse(det_Atlas1 == 1 & det_Atlas2 == 1, 1, 0),
    ext = ifelse(det_Atlas1 == 1 & det_Atlas2 == 0, 1, 0)
  )


# COVARIATES #
# Combo base, diffs
land_covars_all <- grep("_2008_z$|_diff_z$", names(spp_mod_data), value = TRUE)
clim_covars_all <- grep("_38yr_z$|_diffyrs_z$", names(spp_mod_data), value = TRUE)

# Base years
land_covars_base <- grep("_2008_z", names(spp_mod_data), value = TRUE)
clim_covars_base <- grep("_38yr_z", names(spp_mod_data), value = TRUE)

# Diff among ref years 
land_covars_diffs <- grep("_diff_z$", names(spp_mod_data), value = TRUE)
clim_covars_diffs <- grep("_diffyrs_z$", names(spp_mod_data), value = TRUE)



### --- COVAR SELECTION --- ###

# --- DIAGNOSTIC TESTS --- #

# PAIRWISE CORRS: LAND #
Explore_high_corr_grouped <- function(df, covars, threshold = 0.7) {
  
  assign_group <- function(var_name) {
    if (is.na(var_name)) return(NA_character_)
    
    # Define groups in priority order
    group_priority <- list(
      developed_lower  = c("developed_open", "developed_low"),
      developed_upper  = c("developed_med", "developed_high"),
      pasture_crop     = c("pasture", "cropland"),
      developed_total  = c("developed_open", "developed_low", "developed_med", "developed_high"),
      forest_total     = c("forest_deciduous", "forest_evergreen", "forest_mixed"),
      wetlands_total   = c("wetlands_woody", "wetlands_herb")
    )
    
    for(g in names(group_priority)) {
      patterns <- group_priority[[g]]
      if(any(sapply(patterns, function(x) grepl(x, var_name, fixed = TRUE)), na.rm = TRUE)) {
        return(g)
      }
    }
    
    return(NA_character_)
  }
  
  # Subset numeric covariates safely
  numeric_cov <- df %>% dplyr::select(any_of(covars)) %>% na.omit()
  
  if(ncol(numeric_cov) < 2) {
    message("Not enough variables to compute correlations.")
    return(list(pairs = tibble(), summary = tibble()))
  }
  
  # Correlation matrix
  cor_mat <- cor(numeric_cov)
  
  # Upper triangle
  cor_pairs <- which(abs(cor_mat) > threshold & abs(cor_mat) < 1, arr.ind = TRUE)
  
  cor_df <- tibble(
    var1 = colnames(cor_mat)[cor_pairs[,1]],
    var2 = colnames(cor_mat)[cor_pairs[,2]],
    r = cor_mat[cor_pairs]
  )
  
  # Assign groups safely
  cor_df <- cor_df %>%
    rowwise() %>%
    mutate(
      group = {
        g1 <- assign_group(var1)
        g2 <- assign_group(var2)
        if(!is.na(g1) && !is.na(g2) && g1 == g2) g1 else NA_character_
      }
    ) %>%
    ungroup()
  
  # Group summary
  group_summary <- cor_df %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    summarise(
      n_pairs = n(),
      vars_involved = paste(unique(c(var1, var2)), collapse = ", "),
      .groups = "drop"
    )
  
  message("High correlation pairs (>", threshold, "):")
  print(cor_df)
  
  message("\nGroup-level summary of correlated variables:")
  print(group_summary)
  
  return(list(pairs = cor_df, summary = group_summary))
}


# APPLY #
land_corr_info <- Explore_high_corr_grouped(spp_mod_data, land_covars_base, threshold = 0.7)
### From base keep: developed_upper, developed_lower, pasture_crop, wetlands_total, forest_otal OR substypes (and only maybe mixed)

land_corr_info <- Explore_high_corr_grouped(spp_mod_data, land_covars_diffs, threshold = 0.7)
### For diffs keep: developed_upper, developed_lower, forest subtypes OR forest_total, agg pasture_crop


### Ref base year against ref yr diffs 
threshold <- 0.7

# Base + diff covars
base_covars_uncorr <- c(
  "water_open_2008_z", "developed_lower_2008_z", "developed_upper_2008_z",
  "pasture_crop_2008_z", "forest_evergreen_2008_z", "forest_deciduous_2008_z", 
  "forest_mixed_2008_z", "grassland_2008_z", "shrub_scrub_2008_z", 
  "barren_land_2008_z", "wetlands_total_2008_z"
)

diff_covars_uncorr <- c(
  "water_open_diff_z", "pasture_crop_diff_z", 
  "forest_evergreen_diff_z", "forest_deciduous_diff_z",
  "forest_mixed_diff_z", "grassland_diff_z", "wetlands_total_diff_z"
)

all_covars_uncorr <- c(base_covars_uncorr, diff_covars_uncorr)

## --- CORRELATION CHECK --- ##
numeric_cov <- spp_mod_data %>% select(all_of(all_covars_uncorr)) %>% na.omit()
cor_mat <- cor(numeric_cov)

# Flatten correlation matrix
cor_pairs <- as.data.frame(as.table(cor_mat)) %>%
  filter(Var1 != Var2) %>%
  mutate(r_abs = abs(Freq)) %>%
  filter(r_abs > threshold) %>%
  arrange(desc(r_abs))

# Add a suggested "keep" based on ecological group
# Here you can manually or programmatically decide, e.g. always keep one base + one diff per type
cor_pairs <- cor_pairs %>%
  mutate(
    group = case_when(
      grepl("developed_lower", Var1) | grepl("developed_lower", Var2) ~ "developed_lower",
      grepl("developed_upper", Var1) | grepl("developed_upper", Var2) ~ "developed_upper",
      grepl("pasture_crop", Var1) | grepl("pasture_crop", Var2) ~ "pasture_crop",
      grepl("forest_evergreen", Var1) | grepl("forest_evergreen", Var2) ~ "forest_evergreen",
      grepl("grassland", Var1) | grepl("grassland", Var2) ~ "grassland",
      grepl("wetlands_total", Var1) | grepl("wetlands_total", Var2) ~ "wetlands_total",
      TRUE ~ NA_character_
    )
  )

cor_pairs %>%
  select(group, Var1, Var2, Freq) %>%
  arrange(group)


# Covariates to keep from correlation assesment 
# Land
land_base_final <- c(
  "water_open_2008_z",
  "barren_land_2008_z",
  "shrub_scrub_2008_z",
  "developed_lower_2008_z",
  "developed_upper_2008_z",
  "forest_deciduous_2008_z",
  "forest_evergreen_2008_z",
  "forest_mixed_2008_z",
  "wetlands_total_2008_z",
  "grassland_2008_z",
  "pasture_crop_2008_z"
)

land_diff_final <- c(
  "water_open_diff_z",
  "developed_total_diff_z",
  "forest_total_diff_z",     
  "wetlands_total_diff_z",   
  "pasture_crop_diff_z", 
  "grassland_diff_z"
)



# --- PAIRWISE CORRS: CLIMATE --- #

# Combine base + diff covariates
clim_covars <- c(
  "tmin_38yr_z", "tmax_38yr_z", "prcp_38yr_z",
  "tmin_diffyrs_z", "tmax_diffyrs_z", "prcp_diffyrs_z"
)

# Subset numeric data
clim_numeric <- spp_mod_data %>% dplyr::select(all_of(clim_covars)) %>% na.omit()

# Correlation matrix
clim_cor_mat <- cor(clim_numeric)

# Flatten and filter for high correlations ( > 0.7)
cor_pairs <- as.data.frame(as.table(clim_cor_mat)) %>%
  filter(Var1 != Var2) %>%
  mutate(r_abs = abs(Freq)) %>%
  filter(r_abs > 0.7) %>%
  arrange(desc(r_abs))

print(cor_pairs)

# Heatmap
cor_melt <- melt(clim_cor_mat, varnames = c("Var1", "Var2"), value.name = "cor")

ggplot(cor_melt, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limit = c(-1,1), name = "Pearson r") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Heatmap of Climate Covariates")

# All pairs 
cor_pairs_all <- as.data.frame(as.table(clim_cor_mat)) %>%
  rename(Var1 = Var1, Var2 = Var2, cor = Freq) %>%
  arrange(Var1, Var2)

print(cor_pairs_all)



# Base on correlations, keep:
# Base: tmax, prcp; COULD keep tmin (corr to tmax 0.8) if really wanted to preserve
# Diff: Low-mod corr among, could keep all if necessary


clim_base_final <- c("tmin_38yr_z", "tmax_38yr_z", "prcp_38yr_z")
clim_diff_final <- c("tmin_diffyrs_z", "tmax_diffyrs_z", "prcp_diffyrs_z")




#######

# VIF #

# Land - base
land_base_vif <- check_vif(spp_mod_data, land_base_final, threshold = 5)

# Land - diff
land_diff_vif <- check_vif(spp_mod_data, land_diff_final, threshold = 5)

# Climate - base
clim_base_vif <- check_vif(spp_mod_data, clim_base_final, threshold = 5)

# Climate - diff
clim_diff_vif <- check_vif(spp_mod_data, clim_diff_final, threshold = 5)

# Optionally, print summary
list(
  land_base = land_base_vif,
  land_diff = land_diff_vif,
  clim_base = clim_base_vif,
  clim_diff = clim_diff_vif
)


### KEEP ###

### Land ###
### BASE: "water_open_2008_z", "barren_land_2008_z", "shrub_scrub_2008_z",
# "developed_lower_2008_z"**,"developed_upper_2008_z"** (could group into develoepd_total),
# "forest_deciduous_2008_z"**, "forest_evergreen_2008_z"**, "forest_mixed_2008_z"** (could group forest_total),
# "wetlands_total_2008_z", "grassland_2008_z", "pasture_crop_2008_z"
### DIFF: "water_open_diff_z", "forest_total_diff_z", "wetlands_total_diff_z",
# "pasture_crop_diff_z"**, "grassland_diff_z"**, "developed_total_diff_z"**,

### Climate
# BASE: tmax* (OR tmin), prcp
# DIFF: could include all



# SPP-SPECIFIC CANDIDATE RANKING # 

### SPP-SPECIFIC COVARS ##

land_covars_final <- c("water_open_2008_z", "barren_land_2008_z", "shrub_scrub_2008_z",
"developed_total_2008_z", "forest_deciduous_2008_z", "forest_evergreen_2008_z", "forest_mixed_2008_z",
"wetlands_total_2008_z", "grassland_2008_z", "pasture_crop_2008_z", "water_open_diff_z",
"forest_total_diff_z", "wetlands_total_diff_z", "developed_total_diff_z")

clim_covars_final <- c("tmax_38yr_z", "prcp_38yr_z", "tmax_diffyrs_z", "prcp_diffyrs_z")


# --- Function to rank covariates by AIC ---
Explore_covars <- function(df, response_var, covars, max_combo_size = 5) {
  
  results <- list()
  
  # Single covariate models
  for (cov in covars) {
    fmla <- as.formula(paste(response_var, "~", cov))
    mod <- glm(fmla, data = df, family = binomial())
    results[[cov]] <- tibble(
      covariate = cov,
      model_type = "single",
      AIC = AIC(mod)
    )
  }
  
  # Covariate combinations
  if (length(covars) > 1) {
    combos <- unlist(
      lapply(2:max_combo_size, function(k) combn(covars, k, simplify = FALSE)),
      recursive = FALSE
    )
    
    for (combo in combos) {
      combo_name <- paste(combo, collapse = " + ")
      fmla_combo <- as.formula(paste(response_var, "~", combo_name))
      mod_combo <- glm(fmla_combo, data = df, family = binomial())
      results[[combo_name]] <- tibble(
        covariate = combo_name,
        model_type = paste0("combo_", length(combo)),
        AIC = AIC(mod_combo)
      )
    }
  }
  
  # Full model
  if (length(covars) > 2) {
    fmla_all <- as.formula(paste(response_var, "~", paste(covars, collapse = " + ")))
    mod_all <- glm(fmla_all, data = df, family = binomial())
    results[["all_covars"]] <- tibble(
      covariate = "all_covars",
      model_type = "all_covars",
      AIC = AIC(mod_all)
    )
  }
  
  bind_rows(results) %>% mutate(response = response_var)
}

# --- Rank covariates for each species response ---
Rank_covariates <- function(df, land_covars_final, clim_covars_final,
                            responses = c("col", "per", "ext"),
                            max_combo_size = 5) {
  
  rank_candidates <- map_dfr(responses, function(resp) {
    imap_dfr(
      list(land = land_covars_final, climate = clim_covars_final),
      function(cset, set_name) {
        Explore_covars(df, resp, cset, max_combo_size = max_combo_size) %>%
          mutate(covar_set = set_name)
      })
  })
  
  top5_covars <- rank_candidates %>%
    group_by(response, covar_set) %>%
    arrange(AIC, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  list(all = rank_candidates, top5 = top5_covars)
}

# --- Run ranking ---
rank_results <- Rank_covariates(
  df = spp_mod_data,
  land_covars_final = land_covars_final,
  clim_covars_final = clim_covars_final,
  responses = c("col", "per", "ext"),
  max_combo_size = 5
)

# --- View top-ranked covariates ---
rank_results$top5 %>% arrange(response, covar_set, AIC)







## APPLY ## 

# --- STEP 1: filter correlated variables ---
land_covars_clean <- Remove_high_corr(spp_mod_data, land_covars_base)
clim_covars_clean <- Remove_high_corr(spp_mod_data, clim_covars_base)

# --- STEP 2: check VIF only if applicable ---
land_covars_final <- check_vif(spp_mod_data, land_covars_clean$kept)
clim_covars_final <- check_vif(spp_mod_data, clim_covars_clean$kept)  # or use $kept if multiple vars

# --- STEP 3: run model ranking ---
covar_results <- Rank_covariates(
  df = spp_mod_data,
  land_covars_final = land_covars_final$kept,
  clim_covars_final = clim_covars_final$kept,
  responses = c("col", "per", "ext")
)

# --- STEP 4: view results ---
View(covar_results$top5)

















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