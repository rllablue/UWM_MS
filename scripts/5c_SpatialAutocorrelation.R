if (!require("pacman")) install.packages("pacman")
pacman::p_load(psych, usdm)


# Core
library(dplyr)
library(tidyverse)
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
library(lme4)
library(pscl)
library(AICcmodavg)
library(MuMIn)
options(na.action = "na.fail")
library(arm)
library(ncf)
library(DHARMa)
library(spdep)
library(gstat)
library(mgcv)

# Visualization
library(ggplot2)
library(viridis)
library(gt)
library(webshot2)
library(mgcViz)



################################################################################ SPATIAL KERNELS (5.26.26)

######################
### LOCALIZED SAC ###
#####################

# Species-specific occupancy structure based on Atlas 1 occupied source populations 
### within a species-defined dispersal neighborhood.


# 1) Join Data
### Species-specific occupancy data, Atlas block geometry (EPSG:3071)

spp_alpha <- "CAWA"
spp_name <- "Canada Warbler"

spp_blocks_sf <- blocks_rll_sf %>% # join data
  left_join(
    mod_data_all %>%
      filter(alpha_code == spp_alpha) %>%
      dplyr::select(
        atlas_block,
        alpha_code,
        common_name,
        det_Atlas1,
        det_Atlas2,
        transition_state
      ),
    by = "atlas_block"
  )


spp_blocks_sf <- spp_blocks_sf %>% # create centroids
  mutate(
    centroid = st_centroid(geom)
  ) %>%
  mutate(
    x = st_coordinates(centroid)[,1],
    y = st_coordinates(centroid)[,2],
    x_km = x / 1000,
    y_km = y / 1000
  )


# 2) Define Dispersal Neighborhood

# 2a) Neighborhood Size
max_disp_km <- 20 # species-specific max dispersal distance
max_disp_m  <- max_disp_km * 1000 # EPSG:3071 = meters

sigma_km <- 20 # for Gaussian kernels
sigma_m  <- sigma_km * 1000

# 2b) A1 occupancy status (dispersal source)
a1_occu <- spp_blocks_sf$det_Atlas1 # vector



# 3) Define Neighborhood Structure

# 3a) Pairwise Centroid Distances
dist_matrix <- st_distance(
  spp_blocks_sf$centroid
)

dist_matrix <- units::drop_units(dist_matrix)

diag(dist_matrix) <- NA_real_ # remove "self-neighbors"

dist_km <- dist_matrix / 1000


# 3b) Neighborhood
# A1-occupied neighbors within dispersal radius
neighbor_matrix <- (
  !is.na(dist_matrix) &
    dist_matrix <= max_disp_m
) & (
  a1_occu[col(dist_matrix)] == 1
)

# All possible neighbors within dispersal radius
all_neighbor_matrix <- (
  !is.na(dist_matrix) &
    dist_matrix <= max_disp_m
)


# 4) Inverse-Distance Weighting (IDW)
### # w_ij = 1 / d_ij
idw_matrix <- 1 / dist_km

idw_matrix[!neighbor_matrix] <- 0


# 5) Gaussian kernel weighting
kernel_matrix <- exp( # gaussian weights for A1-occupied neighbors
  -(dist_matrix^2) / (2 * sigma_m^2)
)

kernel_matrix[!neighbor_matrix] <- 0


kernel_max_matrix <- exp( # maximum gaussian neighborhood influence 
  -(dist_matrix^2) / (2 * sigma_m^2)
)

kernel_max_matrix[!all_neighbor_matrix] <- 0



# 6) Collect SAC Metrics

# M1) Occupied Neighbor Count
spp_blocks_sf$n_occ_neighbors <- rowSums(
  neighbor_matrix,
  na.rm = TRUE
)


# M2) IDW Connectivity (km^-1)
spp_blocks_sf$idw_connectivity <- rowSums(
  idw_matrix,
  na.rm = TRUE
)


# M3) Raw Gaussian Connectivity
spp_blocks_sf$kernel_connectivity_raw <- rowSums(
  kernel_matrix,
  na.rm = TRUE
)


# M4) Standardized Gaussian Connectivity
kernel_max <- rowSums(
  kernel_max_matrix,
  na.rm = TRUE
)

spp_blocks_sf$kernel_connectivity_std <- ifelse(
  kernel_max > 0,
  spp_blocks_sf$kernel_connectivity_raw / kernel_max,
  0
)


# M5) Distance to Nearest Occupied Source Block (km)
nearest_occ_dist <- apply(
  dist_matrix,
  1,
  function(x) {
    
    vals <- x[
      !is.na(x) &
        x <= max_disp_m &
        a1_occu == 1
    ]
    
    if(length(vals) == 0) {
      return(NA_real_)
    } else {
      return(min(vals))
    }
  }
)

spp_blocks_sf$nearest_occ_dist_km <- (
  nearest_occ_dist / 1000
)


# M6) Mean Distance to Occupied Source Blocks (km)
mean_occ_dist <- apply(
  dist_matrix,
  1,
  function(x) {
    
    vals <- x[
      !is.na(x) &
        x <= max_disp_m &
        a1_occu == 1
    ]
    
    if(length(vals) == 0) {
      return(NA_real_)
    } else {
      return(mean(vals))
    }
  }
)

spp_blocks_sf$mean_occ_dist_km <- (
  mean_occ_dist / 1000
)


# M7) Source Availability Indicator
### Binary: 1 = at least one A1-occupied source block within dispersal radius; 0 = none within radius
spp_blocks_sf$source_available <- ifelse(
  spp_blocks_sf$n_occ_neighbors > 0,
  1,
  0
)




# 7) Diagnostic Plots 

# D1) Occupied neighbor count
### M1: n_occ_neighbors

ggplot(spp_blocks_sf) +
  
  geom_sf(aes(fill = n_occ_neighbors)) +
  
  scale_fill_viridis_c() +
  
  labs(
    title = paste0(
      spp_name,
      ": Occupied Neighbor Count (",max_disp_km, " km)"
    ),
    fill = "Neighbor count"
  ) +
  
  theme_minimal()



# D2) IDW connectivity
### M2: idw_connectivity (km^-1)

ggplot(spp_blocks_sf) +
  
  geom_sf(aes(fill = idw_connectivity)) +
  
  scale_fill_viridis_c() +
  
  labs(
    title = paste0(
      spp_name,
      ": IDW Connectivity (",max_disp_km, " km)"
    ),
    fill = "IDW"
  ) +
  
  theme_minimal()



# D3) Raw Gaussian connectivity
### M3: kernel_connectivity_raw

ggplot(spp_blocks_sf) +
  
  geom_sf(aes(fill = kernel_connectivity_raw)) +
  
  scale_fill_viridis_c() +
  
  labs(
    title = paste0(
      spp_name,
      ": Raw Gaussian Connectivity (",max_disp_km, " km)"
    ),
    fill = "Raw kernel"
  ) +
  
  theme_minimal()



# D4) Standardized Gaussian connectivity (approx. 0-1)
### M4: kernel_connectivity_std

ggplot(spp_blocks_sf) +
  
  geom_sf(aes(fill = kernel_connectivity_std)) +
  
  scale_fill_viridis_c(
    limits = c(0, 1)
  ) +
  
  labs(
    title = paste0(
      spp_name,
      ": Standardized Gaussian Connectivity (",max_disp_km, " km)"
    ),
    fill = "Standardized kernel"
  ) +
  
  theme_minimal()



# D5) Distance to nearest occupied source block (km)
### M5: nearest_occ_dist_km

ggplot(spp_blocks_sf) +
  
  geom_sf(aes(fill = nearest_occ_dist_km)) +
  
  scale_fill_viridis_c(
    na.value = "grey90"
  ) +
  
  labs(
    title = paste0(
      spp_name,
      " Distance to Nearest Occupied Source (",max_disp_km, " km)"
    ),
    fill = "Distance (km)"
  ) +
  
  theme_minimal()



# D6) Mean distance to occupied source blocks (km)
### M6: mean_occ_dist_km

ggplot(spp_blocks_sf) +
  
  geom_sf(aes(fill = mean_occ_dist_km)) +
  
  scale_fill_viridis_c(
    na.value = "grey90"
  ) +
  
  labs(
    title = paste0(
      spp_name,
      ": Mean Distance to Occupied Sources (",max_disp_km, " km)"
    ),
    fill = "Mean distance (km)"
  ) +
  
  theme_minimal()



# D7) Source availability (binary)
### M7: source_available

ggplot(spp_blocks_sf) +
  
  geom_sf(aes(fill = factor(source_available))) +
  
  scale_fill_viridis_d(
    option = "viridis",
    na.value = "grey90"
  ) +
  
  labs(
    title = paste0(
      spp_name,
      ": Source Population Availability (",max_disp_km, " km)"
    ),
    fill = "Source available"
  ) +
  
  theme_minimal()



# D8) SAC histograms

ggplot(
  spp_blocks_sf,
  aes(x = idw_connectivity)
) +
  
  geom_histogram(
    bins = 30
  ) +
  
  labs(
    title = paste0(
      spp_name,
      ": IDW Connectivity Distribution (",max_disp_km, " km)"
    ),
    x = "IDW connectivity"
  ) +
  
  theme_minimal()



# D9) SAC correlations

sac_metrics <- spp_blocks_sf %>%
  
  st_drop_geometry() %>%
  
  dplyr::select(
    n_occ_neighbors,
    idw_connectivity,
    kernel_connectivity_raw,
    kernel_connectivity_std,
    nearest_occ_dist_km,
    mean_occ_dist_km,
    source_available
  )

cor(
  sac_metrics,
  use = "pairwise.complete.obs"
)



# D10) SAC summaries

summary(spp_blocks_sf$n_occ_neighbors)

summary(spp_blocks_sf$idw_connectivity)

summary(spp_blocks_sf$kernel_connectivity_std)

summary(spp_blocks_sf$nearest_occ_dist_km)

summary(spp_blocks_sf$mean_occ_dist_km)






### MODEL SPATIAL TERM ###

# Use std Gaussian kernel (smoothed distance decay function tuned w sigma) 
### as model covariate term to capture influence of plausible source populations
###  on occupancy in a given block.
# Not "just" a nuisance term, but biologically informed population connectivity metric; 
### i.e. not just correcting residuals, modeling dispersal limitation, rescue effects,
### colonization pressure, neighborhood density.


# Define SAC term
sac_covar <- c(
  "kernel_connectivity_std"
)

sac_vars_df <- spp_blocks_sf %>%
  
  st_drop_geometry() %>%
  
  dplyr::select(
    atlas_block,
    all_of(sac_covar)
  )


data_dir_sac <- lapply(
  data_dir,
  function(df) {
    
    left_join(
      df,
      sac_vars_df,
      by = "atlas_block"
    )
  }
)



# Fit models
pa_sac_glm_models <- lapply(
  names(pa_glm_models),
  function(nm) {

    glm(
      formula = pa_sac_formulas[[nm]],
      data    = data_dir_sac[[nm]],
      family  = binomial
    )
  }
)

names(pa_sac_glm_models) <- names(pa_glm_models) 
pa_sac_glm_summaries <- lapply(pa_sac_glm_models, summary) # model output
pa_sac_glm_summaries


pa_sac_predictions <- lapply( # predicted probabilities
  pa_sac_glm_models,
  function(model) {
    
    predict(
      model,
      type = "response"
    )
  }
)


pa_sac_residuals <- lapply( # Pearson's residuals
  pa_sac_glm_models,
  function(model) {
    
    residuals(
      model,
      type = "pearson"
    )
  }
)



# Model comparison
pa_sac_model_comparison <- lapply(
  names(pa_glm_models),
  function(nm) {
    
    base_mod <- pa_glm_models[[nm]]
    sac_mod  <- pa_sac_glm_models[[nm]]
    
    data.frame(
      model = c(
        "Baseline_PA",
        "Localized_SAC"
      ),
      
      AICc = c(
        MuMIn::AICc(base_mod),
        MuMIn::AICc(sac_mod)
      ),
      
      delta_AICc = c(
        0,
        MuMIn::AICc(sac_mod) -
          MuMIn::AICc(base_mod)
      ),
      
      deviance_explained = c(
        1 - (base_mod$deviance / base_mod$null.deviance),
        1 - (sac_mod$deviance  / sac_mod$null.deviance)
      )
    )
  }
)

names(pa_sac_model_comparison) <- names(pa_glm_models)

pa_sac_model_comparison




### SEMIVARIOGRAM RESIDUAL ANALYSIS ###

### DHARMa standardized quantile residuals for semivariogram construction; compare
# spatial structure across true null, pa null, full pa, pa with sac term

# Construct candidate models for comparison
GetResponseName <- function(model_name) {
  
  strsplit(model_name, "_")[[1]][2]
}


# True Null
pa_true_null_models <- lapply(
  names(pa_glm_models),
  function(nm) {
    
    response <- GetResponseName(nm)
    
    glm(
      as.formula(
        paste(response, "~ 1")
      ),
      data   = data_dir[[nm]],
      family = binomial
    )
  }
)

names(pa_true_null_models) <- names(pa_glm_models)


# Null PA
pa_paonly_models <- lapply(
  names(pa_glm_models),
  function(nm) {
    
    response <- GetResponseName(nm)
    
    glm(
      as.formula(
        paste(response, "~ pa_prop")
      ),
      data   = data_dir[[nm]],
      family = binomial
    )
  }
)

names(pa_paonly_models) <- names(pa_glm_models)


# Full PA
pa_full_models <- pa_glm_models


# SAC
pa_localized_sac_models <- pa_sac_glm_models



# Model Sets
model_sets <- list(
  
  true_null    = pa_true_null_models,
  pa_null      = pa_paonly_models,
  full_pa      = pa_full_models,
  localized_sac = pa_localized_sac_models
)


coords_df <- spp_blocks_sf %>%
  st_drop_geometry() %>%
  dplyr::select(atlas_block, x_km, y_km)

# sanity check
stopifnot(!any(is.na(coords_df$x_km)))
stopifnot(!any(is.na(coords_df$y_km)))





# Semivariogram (km space)
set.seed(123)

extract_dharma <- function(models, model_type) {
  
  bind_rows(lapply(names(models), function(nm) {
    
    mod <- models[[nm]]
    
    sim <- DHARMa::simulateResiduals(mod, plot = FALSE)
    
    data.frame(
      atlas_block = data_dir[[nm]]$atlas_block,
      residual    = sim$scaledResiduals,
      model_name  = nm,
      model_type  = model_type
    )
  }))
}


dharma_residuals <- bind_rows(
  extract_dharma(model_sets$true_null,     "true_null"),
  extract_dharma(model_sets$pa_null,       "pa_null"),
  extract_dharma(model_sets$full_pa,       "full_pa"),
  extract_dharma(model_sets$localized_sac, "localized_sac")
)


# attach coords
dharma_residuals <- dharma_residuals %>%
  left_join(coords_df, by = "atlas_block")

# hard safety check
stopifnot(!any(is.na(dharma_residuals$x_km)))
stopifnot(!any(is.na(dharma_residuals$y_km)))


# construct semivariogram
semivar_df <- dharma_residuals %>%
  
  group_by(model_type, model_name) %>%
  
  group_modify(~ {
    
    dat <- .x
    
    # convert to spatial object
    gdf <- sf::st_as_sf(
      dat,
      coords = c("x_km", "y_km"),
      crs = NA_real_
    )
    
    # empirical variogram of residual spatial structure
    v <- gstat::variogram(
      residual ~ 1,
      locations = gdf,
      cutoff = 300,  # km
      width  = 25    # km bins
    )
    
    as.data.frame(v)
  }) %>%
  
  ungroup()






# Visualize
semivar_df <- semivar_df %>%
  
  mutate(
    
    process = case_when(
      grepl("_col$", model_name) ~ "Colonization",
      grepl("_ext$", model_name) ~ "Extirpation",
      TRUE ~ "All transitions"
    ),
    
    model_type = factor(
      model_type,
      levels = c(
        "true_null",
        "pa_null",
        "full_pa",
        "localized_sac"
      ),
      labels = c(
        "True null",
        "PA-only null",
        "Full PA model",
        "Localized SAC model"
      )
    )
  ) %>%
  
  arrange(model_type, model_name, dist) %>%
  
  filter(!is.na(gamma))


# 0-300 km
ggplot(
  semivar_df,
  aes(
    x = dist,
    y = gamma,
    color = model_type,
    group = interaction(model_type, model_name, process)
  )
) +
  
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  
  facet_wrap(~process) +
  
  scale_color_manual(
    values = c(
      "True null" = "grey40",
      "PA-only null" = "orange",
      "Full PA model" = "turquoise4",
      "Localized SAC model" = "purple4"
    )
  ) +
  
  labs(
    title = "Residual Spatial Structure Across Model Types",
    x = "Distance between blocks (km)",
    y = "Semivariance of DHARMa residuals",
    color = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )


# 0-100 km
ggplot(
  semivar_df %>% filter(dist <= 100),
  aes(
    x = dist,
    y = gamma,
    color = model_type,
    group = interaction(model_type, model_name, process)
  )
) +
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  geom_point(size = 1.6, na.rm = TRUE) +
  
  facet_wrap(~process) +
  
  scale_color_manual(
    values = c(
      "True null" = "grey40",
      "PA-only null" = "orange",
      "Full PA model" = "turquoise4",
      "Localized SAC model" = "purple4"
    )
  ) +
  
  labs(
    title = "Residual Spatial Structure (0–100 km)",
    x = "Distance between blocks (km)",
    y = "Semivariance of DHARMa residuals",
    color = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )





# Summarize semivariance
semivar_summary <- semivar_df %>%
  
  group_by(process, model_type) %>%
  
  summarise(
    
    mean_semivariance =
      mean(gamma, na.rm = TRUE),
    
    max_semivariance =
      max(gamma, na.rm = TRUE),
    
    min_semivariance =
      min(gamma, na.rm = TRUE),
    
    semivariance_range =
      max(gamma, na.rm = TRUE) -
      min(gamma, na.rm = TRUE),
    
    .groups = "drop"
  )


semivar_summary




























############################################################## NEW W/ DHAARMa, SEMIVARIOGRAMS

### DATA ###
### Join raw projected lat, lon to data

data_dir_coords <- lapply(names(data_dir), function(nm) {
  
  dat <- data_dir[[nm]]
  
  dat <- dat %>%
    dplyr::left_join(
      covars_raw_rll %>%
        dplyr::select(atlas_block, lon_raw = lon, lat_raw = lat),
      by = "atlas_block"
    ) %>%
    dplyr::mutate(
      x_km = as.numeric(lon_raw) / 1000,
      y_km = as.numeric(lat_raw) / 1000
    ) %>%
    dplyr::filter(!is.na(x_km), !is.na(y_km))
  
  dat
})

names(data_dir_coords) <- names(data_dir)


### Coords df
block_centroids_df <- covars_raw_rll %>%
  dplyr::select(atlas_block, lon, lat) %>%
  distinct() %>%
  mutate(
    x_km = lon / 1000,
    y_km = lat / 1000
  )


block_centroids <- blocks_all_sf %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2],
    x_km = lon / 1000,
    y_km = lat / 1000
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(atlas_block, lon, lat, x_km, y_km)





### --- SPATIAL AUTOCORRELATION --- ###

# Goals:
# 1) Fit null + full PA models
# 2) Extract DHARMa standardized quantile residuals
# 3) Build semivariograms from DHARMa residuals


### TRUE NULL, PA NULL, FULL MODEL COMPARISON ###

# Helper: extract response variable name
GetResponseName <- function(model_name) {
  strsplit(model_name, "_")[[1]][2]
}



### TRUE NULL (response ~ 1)
pa_true_null_models <- lapply(names(pa_glm_models), function(nm) {
  
  response <- GetResponseName(nm)
  
  glm(
    as.formula(paste(response, "~ 1")),
    data = data_dir[[nm]],
    family = binomial
  )
})

names(pa_true_null_models) <- names(pa_glm_models)



### PA-ONLY NULL (response ~ pa_prop)
pa_paonly_models <- lapply(names(pa_glm_models), function(nm) {
  
  response <- GetResponseName(nm)
  
  glm(
    as.formula(paste(response, "~ pa_prop")),
    data = data_dir[[nm]],
    family = binomial
  )
})

names(pa_paonly_models) <- names(pa_glm_models)



### FULL MODEL (response ~ pa_prop + env, etc. covariates)
pa_full_models <- pa_glm_models





### EXTRACT, MAP DHARMa RESIDUALS ###
set.seed(123)

model_sets <- list(
  true_null = pa_true_null_models,
  pa_null   = pa_paonly_models,
  full      = pa_full_models
)


dharma_residuals <- lapply(names(model_sets), function(set_name) {
  
  models <- model_sets[[set_name]]
  
  out <- lapply(names(models), function(nm) {
    
    sim <- simulateResiduals(
      fittedModel = models[[nm]],
      plot = FALSE
    )
    
    data.frame(
      atlas_block = data_dir_coords[[nm]]$atlas_block,
      x_km = data_dir_coords[[nm]]$x_km,
      y_km = data_dir_coords[[nm]]$y_km,
      residual = sim$scaledResiduals,
      species_model = nm,
      model_type = set_name
    )
  })
  
  bind_rows(out)
})

dharma_residuals <- bind_rows(dharma_residuals)





### SEMIVARIOGRTAMS ###

semivar_df <- dharma_residuals %>%
  
  group_by(species_model, model_type) %>%
  
  group_modify(~ {
    
    dat <- .x
    
    coordinates <- dat[, c("x_km", "y_km")]
    
    gdf <- sf::st_as_sf(
      dat,
      coords = c("x_km", "y_km"),
      crs = 3071
    )
    
    vgm_obj <- gstat::variogram(
      residual ~ 1,
      locations = gdf,
      cutoff = 300,
      width = 25
    )
    
    as.data.frame(vgm_obj)
    
  }) %>%
  
  ungroup()



### PLOT ###

# Tidy
semivar_df <- semivar_df %>%
  
  mutate(
    
    process = case_when(
      grepl("_col$", species_model) ~ "Colonization",
      grepl("_ext$", species_model) ~ "Extirpation"
    ),
    
    model_type = factor(
      model_type,
      levels = c("true_null", "pa_null", "full"),
      labels = c(
        "True null",
        "PA-only null",
        "Full model"
      )
    )
  )



# Plot
ggplot(
  semivar_df,
  aes(x = dist, y = gamma, color = model_type)
) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 2) +
  
  facet_wrap(~process) +
  
  scale_color_manual(
    values = c(
      "True null" = "grey40",
      "PA-only null" = "orange",
      "Full model" = "turquoise4"
    )
  ) +
  
  labs(
    title = "Residual Spatial Structure Across Model Types",
    x = "Distance (km)",
    y = "Semivariance",
    color = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )



# Summarize

semivar_summary <- semivar_df %>%
  
  group_by(process, model_type) %>%
  
  summarise(
    mean_semivariance = mean(gamma, na.rm = TRUE),
    max_semivariance = max(gamma, na.rm = TRUE),
    .groups = "drop"
  )

semivar_summary


################################################################################
################# SMOOTH LATENT SPATIAL PROCESS ################################


### Fit SAC-smoothing GAM and compare fit, residuals with GLM

pa_full_models <- pa_glm_models


# Fit GAMs
spatial_gam_models <- lapply(names(pa_full_models), function(nm) {
  
  # Original GLM
  glm_mod <- pa_full_models[[nm]]
  
  # Original formula text
  base_formula <- formula(glm_mod)
  
  # Add spatial smooth
  gam_formula <- update(
    base_formula,
    . ~ . + s(x_km, y_km, k = 30) # k is flexible, controls "wiggle"
  )
  
  gam(
    formula = gam_formula,
    data    = data_dir_coords[[nm]],
    family  = binomial(link = "logit"),
    method  = "REML"
  )
})

names(spatial_gam_models) <- names(pa_full_models)


spatial_gam_summaries <- lapply(
  spatial_gam_models,
  summary
)

spatial_gam_summaries



### Model Comparison
model_comparison <- lapply(names(pa_full_models), function(nm) {
  
  glm_mod <- pa_full_models[[nm]]
  gam_mod <- spatial_gam_models[[nm]]
  
  data.frame(
    model = c("GLM", "Spatial_GAM"),
    
    AIC = c(
      MuMIn::AICc(glm_mod),
      MuMIn::AICc(gam_mod)
    ),
    
    deviance_explained = c(
      1 - glm_mod$deviance / glm_mod$null.deviance,
      
      summary(gam_mod)$dev.expl
    )
  )
})

names(model_comparison) <- names(pa_full_models)

model_comparison



### Residual Comparison
set.seed(123)

gam_dharma_residuals <- lapply(names(spatial_gam_models), function(nm) {
  
  sim <- simulateResiduals(
    fittedModel = spatial_gam_models[[nm]],
    plot = FALSE
  )
  
  data.frame(
    atlas_block = data_dir_coords[[nm]]$atlas_block,
    x_km = data_dir_coords[[nm]]$x_km,
    y_km = data_dir_coords[[nm]]$y_km,
    residual = sim$scaledResiduals,
    species_model = nm,
    model_type = "spatial_gam"
  )
})

gam_dharma_residuals <- bind_rows(gam_dharma_residuals)

# compare with GLM residuals
dharma_residuals_compare <- bind_rows(
  dharma_residuals,
  gam_dharma_residuals
)


### Semivariograms

semivar_comp_df <- dharma_residuals_compare %>%
  
  group_by(species_model, model_type) %>%
  
  group_modify(~ {
    
    dat <- .x
    
    gdf <- sf::st_as_sf(
      dat,
      coords = c("x_km", "y_km"),
      crs = 3071
    )
    
    vgm_obj <- gstat::variogram(
      residual ~ 1,
      locations = gdf,
      cutoff = 300,
      width = 25
    )
    
    as.data.frame(vgm_obj)
    
  }) %>%
  
  ungroup()



### Visualize

### VISUALIZE ###

# Clean labels
semivar_comp_df <- semivar_comp_df %>%
  
  mutate(
    
    process = case_when(
      grepl("_col$", species_model) ~ "Colonization",
      grepl("_ext$", species_model) ~ "Extirpation"
    ),
    
    model_type_label = case_when(
      model_type == "true_null"   ~ "True null",
      model_type == "pa_null"     ~ "PA-only null",
      model_type == "full"        ~ "Full GLM",
      model_type == "spatial_gam" ~ "Spatial GAM"
    )
  )



### PLOT ###

ggplot(
  semivar_comp_df,
  aes(
    x = dist,
    y = gamma,
    color = model_type_label,
    group = interaction(model_type_label, species_model)
  )
) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 2) +
  
  facet_wrap(~process) +
  
  scale_color_manual(
    values = c(
      "True null"    = "grey40",
      "PA-only null" = "orange",
      "Full GLM"     = "turquoise4",
      "Spatial GAM"  = "darkorchid"
    )
  ) +
  
  labs(
    title = "Residual Spatial Structure Across Model Types",
    x = "Distance (km)",
    y = "Semivariance",
    color = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )



# Summarize

semivar_comp_summary <- semivar_comp_df %>%
  
  group_by(process, model_type_label) %>%
  
  summarise(
    mean_semivariance = mean(gamma, na.rm = TRUE),
    max_semivariance = max(gamma, na.rm = TRUE),
    .groups = "drop"
  )

semivar_comp_summary












### SPATIAL AUTOCOVARIATE TERM #################################################
### Using distance-beased weights

### V1

coords_mat <- as.matrix(
  coords_df[, c("x_km", "y_km")]
)

# neighbors within 100 km; set ditance
nb <- dnearneigh(
  coords_mat,
  d1 = 0,
  d2 = 100
)


distances <- nbdists(nb, coords_mat)

inv_dist <- lapply(
  distances,
  function(x) 1 / x
)


inv_dist <- lapply(
  distances,
  function(x) 1 / (x + 0.001)
)


dat$auto_cov <- lag.listw(
  lw,
  dat[[response]],
  zero.policy = TRUE
)







## DATA ####################################################### OLD WIP W/ PEARSON'S RES


### --- GLOBAL AUTOCORRELATION --- ###
# --> Is there spatial structure in residuals?

### MORAN'S I ###
### + kNN sensitivity analysis

k_seq <- c(3, 5, 8, 10, 15, 20, 30)

pa_moran_k <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir_coords[[nm]]
  res <- pa_residuals[[nm]]
  coords <- cbind(dat$x_km, dat$y_km)
  
  out <- lapply(k_seq, function(k_val) {
    
    nb <- knn2nb(knearneigh(coords, k = k_val))
    lw <- nb2listw(nb, style = "W")
    
    test <- moran.test(res, lw)
    
    data.frame(
      model = nm,
      k = k_val,
      Moran_I = unname(test$estimate[1]),
      p_value = test$p.value
    )
  })
  
  bind_rows(out)
})

pa_moran_k <- bind_rows(pa_moran_k)




# visualize neighbor-definition sensitivity
ggplot(pa_moran_k, aes(x = k, y = Moran_I, color = model)) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  scale_x_continuous(breaks = k_seq) +
  
  scale_color_manual(
    values = c(
      "RLL_col" = "turquoise",
      "RLL_ext" = "tomato"
    )
  ) +
  
  labs(
    title = "Moran's I Sensitivity to Neighborhood Size (kNN)",
    x = "k (number of neighbors)",
    y = "Moran's I"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )







### CORRELOGRAM: BINNED ###
pa_correlog <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir_coords[[nm]]
  res <- pa_residuals[[nm]]
  
  ncf::correlog(
    x = dat$x_km,
    y = dat$y_km,
    z = res,
    increment = 50,
    resamp = 999
  )
})

names(pa_correlog) <- names(pa_glm_models)



correlog_df <- bind_rows(lapply(names(pa_correlog), function(nm) {
  
  obj <- pa_correlog[[nm]]
  
  data.frame(
    model    = nm,
    distance = obj$mean.of.class,
    moran_I  = obj$correlation,
    p_value  = obj$p,
    n_pairs  = obj$n
  )
}))



ggplot(correlog_df, 
       aes(x = distance, y = moran_I, color = model)) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  scale_color_manual(
    values = c(
      "RLL_col" = "turquoise",
      "RLL_ext" = "tomato"
    ),
    labels = c(
      "Colonization",
      "Extirpation"
    )
  ) +
  
  labs(
    title = "Binned Spatial Correlogram",
    x = "Distance (km)",
    y = "Moran's I",
    color = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )







### CORRELOG: SPLINE ###
pa_spline <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir_coords[[nm]]
  res <- pa_residuals[[nm]]
  
  ncf::spline.correlog(
    x = dat$x_km,
    y = dat$y_km,
    z = res,
    resamp = 999
  )
})

names(pa_spline) <- names(pa_glm_models)



spline_df <- bind_rows(lapply(names(pa_spline), function(nm) {
  
  obj <- pa_spline[[nm]]
  
  # extract safely as vectors
  dist_vals  <- as.numeric(obj$real$predicted$x)
  corr_vals  <- as.numeric(obj$real$predicted$y)
  lower_vals <- as.numeric(obj$boot$boot.summary$predicted$y["0.025", ])
  upper_vals <- as.numeric(obj$boot$boot.summary$predicted$y["0.975", ])
  
  data.frame(
    model = nm,
    distance = dist_vals,
    corr = corr_vals,
    lower = lower_vals,
    upper = upper_vals
  )
}))


# OPT 1: half-decay distance
### what is distance where corr is mostly gone?
### half-decay: distance where corr drops to half its max
half_decay_range <- spline_df %>%
  group_by(model) %>%
  summarise(
    peak_I = max(corr, na.rm = TRUE),
    target = peak_I / 2,
    
    range_km = distance[which.min(abs(corr - target))]
  )


  
# OPT 2: Effective spatial independence range
### first distance where CI overlaps 0 for the final time
zero_range <- spline_df %>%
  group_by(model) %>%
  arrange(distance) %>%
  mutate(
    non_sig = (lower <= 0 & upper >= 0),
    was_sig = !non_sig,
    transition = was_sig & lead(non_sig)
  ) %>%
  summarise(
    range_km = ifelse(any(transition),
                      distance[which(transition)[1]],
                      NA)
  )


# Visualization
ggplot(spline_df, aes(x = distance, y = corr, color = model)) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = model),
    alpha = 0.2,
    color = NA
  ) +
  
  geom_line(linewidth = 1.2) +
  
  geom_vline(
    data = half_decay_range,
    aes(xintercept = range_km, color = "black"),
    linetype = "dotted",
    linewidth = 1
  ) +
  
  geom_vline(
    data = zero_range,
    aes(xintercept = range_km, color = "green"),
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "RLL_col" = "turquoise",
    "RLL_ext" = "tomato"
  )) +
  
  scale_fill_manual(values = c(
    "RLL_col" = "turquoise",
    "RLL_ext" = "tomato"
  )) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  
  labs(
    title = "Spline Spatial Correlogram",
    x = "Distance (km)",
    y = "Correlation"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )



# Zoomed ggplot
ggplot(spline_df, aes(x = distance, y = corr, color = model)) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = model),
    alpha = 0.2,
    color = NA
  ) +
  
  geom_line(linewidth = 1.2) +
  
  geom_vline(
    data = half_decay_range,
    aes(xintercept = range_km),
    linetype = "dotted",
    linewidth = 1,
    color = "orchid"
  ) +
  
  geom_vline(
    data = zero_range,
    aes(xintercept = range_km),
    linetype = "dashed",
    linewidth = 1,
    color = "springgreen"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "RLL_col" = "turquoise",
    "RLL_ext" = "tomato"
  )) +
  
  scale_fill_manual(values = c(
    "RLL_col" = "turquoise",
    "RLL_ext" = "tomato"
  )) +
  
  coord_cartesian(xlim = c(0, 200)) +
  
  labs(
    title = "Spline Spatial Correlogram (0–200 km)",
    x = "Distance (km)",
    y = "Correlation"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )












### LISA ###

k <- 6

lisa_list <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir[[nm]]
  res <- pa_residuals[[nm]]
  
  coords <- cbind(dat$lon, dat$lat)
  
  knn <- knearneigh(coords, k = k)
  nb  <- knn2nb(knn)
  lw  <- nb2listw(nb, style = "W")
  
  lisa <- spdep::localmoran(res, lw)
  
  data.frame(
    model = nm,
    lon   = dat$lon,
    lat   = dat$lat,
    residual = res,
    local_I  = lisa[, 1],
    p_value  = lisa[, 5]
  )
})


lisa_df <- do.call(rbind, lisa_list)


lagged_res <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir[[nm]]
  res <- pa_residuals[[nm]]
  
  coords <- cbind(dat$lon, dat$lat)
  
  knn <- knearneigh(coords, k = k)
  nb  <- knn2nb(knn)
  lw  <- nb2listw(nb, style = "W")
  
  lag <- lag.listw(lw, res)
  
  data.frame(
    model = nm,
    lag_residual = lag
  )
})

lag_df <- do.call(rbind, lagged_res)

lisa_df$lag_residual <- lag_df$lag_residual


lisa_df$cluster <- "Not significant"

sig <- lisa_df$p_value < 0.05

lisa_df$cluster[sig & lisa_df$residual > 0 & lisa_df$lag_residual > 0] <- "High-High"
lisa_df$cluster[sig & lisa_df$residual < 0 & lisa_df$lag_residual < 0] <- "Low-Low"
lisa_df$cluster[sig & lisa_df$residual > 0 & lisa_df$lag_residual < 0] <- "High-Low"
lisa_df$cluster[sig & lisa_df$residual < 0 & lisa_df$lag_residual > 0] <- "Low-High"




ggplot(lisa_df, aes(x = lon, y = lat, fill = cluster)) +
  
  geom_point(shape = 21, size = 1.6, color = "transparent") +
  
  facet_wrap(~model) +
  
  scale_fill_manual(
    values = c(
      "High-High" = "#d7191c",
      "Low-Low"   = "#2c7bb6",
      "High-Low"  = "#fdae61",
      "Low-High"  = "#abd9e9",
      "Not significant" = "grey85"
    )
  ) +
  
  coord_equal() +
  theme_void()











### AUTOCOVARIATE MODEL TERM ###

# use the same k / neighborhood logic as your Moran + LISA work
# here we use inverse-distance style autocovariance via autocov_dist()

pa_autocov <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir[[nm]]
  
  response <- strsplit(nm, "_")[[1]][2]
  y <- dat[[response]]
  
  coords <- cbind(
    dat$lon,
    dat$lat
  )
  
  autocov <- spdep::autocov_dist(
    y,
    coords,
    nbs = 1,
    type = "inverse",
    style = "W",
    zero.policy = TRUE
  )
  
  dat$autocov_pa <- autocov
  
  dat
})

names(pa_autocov) <- names(pa_glm_models)



pa_glm_models_auto <- lapply(names(pa_autocov), function(nm) {
  
  response <- strsplit(nm, "_")[[1]][2]
  
  glm(
    formula = update(
      ExtractPAFormula(reference_global_models[[nm]], response),
      . ~ . + autocov_pa
    ),
    data = pa_autocov[[nm]],
    family = binomial
  )
})

names(pa_glm_models_auto) <- names(pa_autocov)
summary(pa_glm_models_auto)





### STATISTICAL FILTERS — PA MODEL + SPATIAL AUTOCOVARIATE ###

### predictor set used in final PA + SAC model
### (base predictors from selected PA model + pa_prop + autocov_pa)

GetPAAutoCovars <- function(model_obj) {
  
  # pull predictor names directly from fitted glm
  vars <- names(coef(model_obj))
  
  # remove intercept
  vars <- vars[vars != "(Intercept)"]
  
  vars
}


### VISUALIZE CORRELATIONS INCLUDING AUTOCOVARIATE

pa_auto_corr_results <- lapply(names(pa_glm_models_auto), function(nm) {
  
  dat <- pa_autocov[[nm]]
  
  covs_use <- GetPAAutoCovars(pa_glm_models_auto[[nm]])
  
  VisualizeCorrelations(
    data  = dat,
    covs  = covs_use,
    label = paste(nm, "+ autocov")
  )
})

names(pa_auto_corr_results) <- names(pa_glm_models_auto)


### NUMERICALLY IDENTIFY HIGH CORRELATIONS
### especially important for checking autocov_pa inflation

pa_auto_highcorrs <- lapply(names(pa_glm_models_auto), function(nm) {
  
  dat <- pa_autocov[[nm]]
  
  covs_use <- GetPAAutoCovars(pa_glm_models_auto[[nm]])
  
  cat("\n============================\n")
  cat("High correlations for:", nm, "\n")
  cat("============================\n")
  
  GetHighCorrs(
    data = dat,
    covs = covs_use,
    threshold = 0.7
  )
})

names(pa_auto_highcorrs) <- names(pa_glm_models_auto)


### VIF CHECK INCLUDING AUTOCOVARIATE

vif_pa_models_auto <- lapply(pa_glm_models_auto, car::vif)

vif_pa_models_auto


### Example:
vif_pa_models_auto[["RLL_col"]]
vif_pa_models_auto[["RLL_ext"]]












# MODEL COMPARISON #
### Check performance of spatial vs non-spatial term model
AIC(pa_glm_models[["RLL_col"]])
AIC(pa_glm_models_auto[["RLL_col"]])

summary(pa_glm_models_auto[["RLL_col"]])



res_auto <- residuals(
  pa_glm_models_auto[["RLL_col"]],
  type = "pearson"
)

dat <- pa_autocov[["RLL_col"]]

coords <- cbind(dat$lon, dat$lat)

knn <- knearneigh(coords, k = 6)
nb  <- knn2nb(knn)
lw  <- nb2listw(nb, style = "W")

moran.test(
  res_auto,
  lw,
  zero.policy = TRUE
)



# AUC
pa_auc_auto <- lapply(names(pa_glm_models_auto), function(nm) {
  
  model <- pa_glm_models_auto[[nm]]
  
  # extract observed response from model frame
  observed <- model$model[[1]]
  
  # predicted probabilities
  predicted <- predict(model, type = "response")
  
  # calculate AUC
  auc_value <- pROC::auc(observed, predicted)
  
  data.frame(
    model = nm,
    AUC = as.numeric(auc_value)
  )
})

pa_auc_auto <- do.call(rbind, pa_auc_auto)
pa_auc_auto



# Pseudo-R2
pa_r2_auto <- lapply(names(pa_glm_models_auto), function(nm) {
  
  model <- pa_glm_models_auto[[nm]]
  
  r2_vals <- pscl::pR2(model)
  
  data.frame(
    model = nm,
    McFadden_R2 = r2_vals["McFadden"]
  )
})

pa_r2_auto <- do.call(rbind, pa_r2_auto)
pa_r2_auto


# Combine
pa_model_metrics_auto <- merge(
  pa_auc_auto,
  pa_r2_auto,
  by = "model"
)

pa_model_metrics_auto
