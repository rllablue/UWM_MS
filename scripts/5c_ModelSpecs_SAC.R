### WIP NOTES
# 5.28: Move to Spatial Kernels workflow
# 6.3: Complete SAC workflow; select covariate term, run models and diagnostics
# 6.8: Build out max dispersal distance step, i.e. run for multiple km buffers and 
# provide stats justification for selecting scale of effect (include as supplementary material)


### Script: 5a_ModelSpecs_SAC.R
### Purpose:


######################
### --- SET UP --- ###
######################


### --- LOAD PACKAGES --- ###

if (!require("pacman")) install.packages("pacman")
pacman::p_load(psych, 
               usdm,
               )

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



### --- DATA WRANGLING --- ###

# A) Join species occupancy, Atlas block geometry 
### Native geometry: Wisconsin Transverse Mercator (EPSG:3071, units = meters)
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


# B) Represent blocks via their centroids 
### Necessary for pairwise distance comparison
spp_blocks_sf <- spp_blocks_sf %>% 
  mutate(
    centroid = st_centroid(geom)
  ) %>%
  mutate(
    x = st_coordinates(centroid)[,1],
    y = st_coordinates(centroid)[,2],
    x_km = x / 1000,
    y_km = y / 1000
  )




#############################
### --- SPATIAL SPECS --- ###
#############################

### Block relationships and spatial weights are calculated from different subsets
# among colonization, extirpation responses. These weights and the subsequently
# derived SAC term must therefore be created for each species with this script.


### --- A) Define biological neighborhood structure

# Identify block occupancy status in A1
### i.e. identify potential source populations for species
a1_occu <- spp_blocks_sf$det_Atlas1 # vector


# Species-specific maximum dispersal radius
### constrain biologically plausible neighborhood extent 
max_disp_km <- 20
max_disp_m  <- max_disp_km * 1000 # EPSG:3071



### --- B) Define pairwise, statistical neighborhood structure 

# Pairwise centroid distances
dist_matrix <- st_distance(
  spp_blocks_sf$centroid
)

dist_matrix <- units::drop_units(dist_matrix)


# Remove "self-as-neighbor" relationship
diag(dist_matrix) <- NA_real_


# Transform distance matrix (m to km) for plot legibility
### Atlas blocks = 5 x 5 x km (scale of model inference)
dist_km <- dist_matrix / 1000


# Define criteria for "source-neighbor" status
### A1-occupied blocks within max dispersal radius
neighbor_matrix <- (
  !is.na(dist_matrix) &
    dist_matrix <= max_disp_m
) & (
  a1_occu[col(dist_matrix)] == 1
)


# All blocks within max dispersal radius
### i.e. source-neighbor and non-source-neighbors
all_neighbor_matrix <- (
  !is.na(dist_matrix) &
    dist_matrix <= max_disp_m
)




#########################
### --- SAC MODEL --- ###
#########################

### --- DEFINE MODEL TERM --- ###

### --- A) Raw Gaussian Kernel Weighting
### Smoothed distance-decay function tuned by sigma; influence declines continuously
### w increased distance; decay process less abrupt than IDW.
### sigma controls rate of decay, i.e. neighborhood smoothness, connectivity
# small sigma = adjacent populations influence greatly; larger = 
sigma_km <- 15
sigma_m  <- sigma_km * 1000

kernel_matrix <- exp(-(dist_matrix^2) / (2 * sigma_m^2))
kernel_matrix[!neighbor_matrix] <- 0 # only source-neighbors

spp_blocks_sf$kernel_connectivity_raw <- rowSums(kernel_matrix, na.rm = TRUE) 


### --- B) Standardized Gaussian Kernel Weighting
# Maximum Gaussian neighborhood influence, i.e. If all neighboring blocks (within
# dispersal radius) were A1-occupied, then what would be the max possible spatial 
# influence on a given block?
kernel_max_matrix <- exp(-(dist_matrix^2) / (2 * sigma_m^2))
kernel_max_matrix[!all_neighbor_matrix] <- 0

kernel_max <- rowSums(kernel_max_matrix, na.rm = TRUE)

spp_blocks_sf$kernel_connectivity_std <- ifelse(
  kernel_max > 0,
  spp_blocks_sf$kernel_connectivity_raw / kernel_max,
  0
)




### --- FIT MODEL --- ###

### --- A) Construct model

### A1) Append z-scaled SAC covariate to data directory/each model covariate pool
data_dir$RLL_col <- data_dir$RLL_col %>%
  left_join(
    spp_blocks_sf %>% st_drop_geometry() %>% 
      dplyr::select(atlas_block, kernel_connectivity_std),
    by = "atlas_block"
  ) %>%
  mutate(kernel_std = as.numeric(scale()))

data_dir$RLL_ext <- data_dir$RLL_ext %>%
  left_join(
    spp_blocks_sf %>% st_drop_geometry() %>% 
      dplyr::select(atlas_block, kernel_connectivity_std),
    by = "atlas_block"
  )




### A2) Update model formula
## Helper: parse PA model string, add SAC covariate
ExtractSACFormula <- function(ref_df, response) {
  
  if (is.null(ref_df) || nrow(ref_df) == 0) {
    stop("Empty reference model for response: ", response)
  }
  
  as.formula(
    paste(response, "~", ref_df$Modnames[1], "+ kernel_std")
  )
}


sac_glm_models <- lapply(names(pa_glm_models), function(nm) {
  
  response <- strsplit(nm, "_")[[1]][2]
  
  glm(
    formula = ExtractSACFormula(pa_glm_models[[nm]], response),
    data    = data_dir[[nm]],
    family  = binomial
  )
})

names(sac_glm_models) <- names(pa_glm_models)

sac_glm_summaries <- lapply(sac_glm_models, summary)
sac_glm_summaries








### --- EVALUATE MODEL --- ###









### SCRIPT NOTE: SAC covariate term was derived, assessed, and selected from below workflow.

#############################
### --- ASSESSING SAC --- ###
#############################

### Neighborhood structure 1) conditioned on proximity to Atlas 1-occupied blocks 
# (i.e. source populations w dispersal potential) and 2) within a species-specific 
# maximum dispersal radius.

### Steps
# 1) Data wrangling, spatial specs
# 2) Pairwise neighborhood structure
# 3) Pairwise spatial weighting schemes
# 4) Global SAC testing
# 5) Block-level SAC metrics / covariates



### --- SPATIAL SPECS --- ###

### --- A) Quantify source-neighbor metrics

# A1) Count-based source-neighbor metrics
spp_blocks_sf$n_occ_neighbors <- rowSums(
  neighbor_matrix, 
  na.rm = TRUE
)


# A2) Binary source-neighbor metrics
### 1 = at least one A1-occupied source block within dispersal radius; 
# 0 = none within radius; NOT a pairwise metric
spp_blocks_sf$source_available <- ifelse(
  spp_blocks_sf$n_occ_neighbors > 0,
  1,
  0
)


# A3) Raw distance-based source-neighbor metrics (km)
nearest_occ_dist <- apply(dist_matrix, 1, function(x) { # dist to nearest source-neighbor
  vals <- x[!is.na(x) & x <= max_disp_m & a1_occu == 1]
  if(length(vals) == 0) NA_real_ else min(vals)
})

mean_occ_dist <- apply(dist_matrix, 1, function(x) { # mean dist to nearest source-neighbor
  vals <- x[!is.na(x) & x <= max_disp_m & a1_occu == 1]
  if(length(vals) == 0) NA_real_ else mean(vals)
})

spp_blocks_sf$nearest_occ_dist_km <- nearest_occ_dist / 1000
spp_blocks_sf$mean_occ_dist_km    <- mean_occ_dist / 1000



### A4) Weighted distance-based source-neighbor metrics
# matrices define how neighboring blocks influence each other spatially, i.e. 
# pairwise relationships (NOT block-level statistics)

# A4a) Inverse-distance weighting (IDW)
### weight function: w_ij = 1 / d_ij
### Nearby, occupied source populations more strongly influence block-level response
### than more distance source populations (un-smoothed)
idw_matrix <- 1 / dist_km
idw_matrix[!neighbor_matrix] <- 0

spp_blocks_sf$idw_connectivity <- rowSums(idw_matrix, na.rm = TRUE)


# A4b) Raw Gaussian kernel weighting
### Smoothed distance-decay function tuned by sigma; influence declines continuously
### w increased distance; decay process less abrupt than IDW.
### sigma controls rate of decay, i.e. neighborhood smoothness, connectivity
# small sigma = adjacent populations influence greatly; larger = 
sigma_km <- 15
sigma_m  <- sigma_km * 1000

kernel_matrix <- exp(-(dist_matrix^2) / (2 * sigma_m^2))
kernel_matrix[!neighbor_matrix] <- 0 # only source-neighbors

spp_blocks_sf$kernel_connectivity_raw <- rowSums(kernel_matrix, na.rm = TRUE) 

# Standardized Gaussian kernel weighting
# Maximum Gaussian neighborhood influence, i.e. If all neighboring blocks (within
# dispersal radius) were A1-occupied, then what would be the max possible spatial 
# influence on a given block?
kernel_max_matrix <- exp(-(dist_matrix^2) / (2 * sigma_m^2))
kernel_max_matrix[!all_neighbor_matrix] <- 0

kernel_max <- rowSums(kernel_max_matrix, na.rm = TRUE)

spp_blocks_sf$kernel_connectivity_std <- ifelse(
  kernel_max > 0,
  spp_blocks_sf$kernel_connectivity_raw / kernel_max,
  0
)




### --- GLOBAL SAC DIAGNOSTICS --- ###
### NOT models; spatial structure checks

### --- A) Quantitative diagnostics

# A1) Correlations among metrics
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

sac_cor_matrix <- cor(sac_metrics, use = "pairwise.complete.obs")
round(sac_cor_matrix, 2)

sac_cor_pairs <- as.data.frame(as.table(sac_cor_matrix)) %>% # pairwise cor table
  
  rename(
    metric_1 = Var1,
    metric_2 = Var2,
    correlation = Freq
  ) %>%
  
  filter(metric_1 != metric_2) %>% # remove self-correlations
  
  rowwise() %>% # remove duplicate pairs
  
  mutate(
    pair_id = paste(sort(c(metric_1, metric_2)), collapse = "__")
  ) %>%
  
  ungroup() %>%
  
  distinct(pair_id, .keep_all = TRUE) %>%
  
  dplyr::select(
    -pair_id
  ) %>%
  
  filter(
    abs(correlation) >= 0.7 # cor threshold
  ) %>%
  
  arrange(desc(correlation))

sac_cor_pairs


# A2) Univariate distributions
summary(spp_blocks_sf$n_occ_neighbors)
summary(spp_blocks_sf$idw_connectivity)
summary(spp_blocks_sf$kernel_connectivity_std)
summary(spp_blocks_sf$nearest_occ_dist_km)
summary(spp_blocks_sf$mean_occ_dist_km)



# B) Diagnostic Plots 

# DP1) Occupied neighbor count
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



# DP2) IDW connectivity
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



# DP3) Raw Gaussian connectivity
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



# DP4) Standardized Gaussian connectivity (approx. 0-1)
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



# DP5) Distance to nearest occupied source block (km)
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



# DP6) Mean distance to occupied source blocks (km)
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



# DP7) Source availability (binary)
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



# DP8) SAC histograms

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





###################################
### --- MODEL SPATIAL TERM --- ###
##################################

### Select representative term/s to include in model as biologically informed covariate
# representing/modeling population connectivity (i.e. block susceptibility to dispersal, propogule
# pressure, rescue effects, etc.)

### Used here: kernel_connectivity_std (std gaussian kernel weighting)

# A) Data wrangling
### create df of all SAC terms for ease of substitution
data_dir_sac <- lapply(
  data_dir,
  function(df) {
    
    left_join(
      df,
      spp_blocks_sf %>%
        st_drop_geometry() %>%
        dplyr::select(
          atlas_block,
          x,
          y,
          source_available,
          n_occ_neighbors,
          nearest_occ_dist_km,
          mean_occ_dist_km,
          idw_connectivity,
          kernel_connectivity_raw,
          kernel_connectivity_std
        ),
      by = "atlas_block"
    )
  }
)


# B) Fit, index comparison models
# Helper: Extract response name
GetResponseName <- function(model_name) {
  
  strsplit(model_name, "_")[[1]][2]
}


# B1) Null models
null_models <- lapply(
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

names(null_models) <- names(pa_glm_models)


# B2) PA-only models
paonly_models <- lapply(
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

names(paonly_models) <- names(pa_glm_models)


# B3) Full environmental-PA models: pa_glm_models
envpa_models <- pa_glm_models


# C) Fit model w SAC term
### append term to best full models from previous script (pa_glm_models)
sac_envpa_models <- lapply(
  names(pa_glm_models),
  function(nm) {

    glm(
      formula = update(
        formula(pa_glm_models[[nm]]),
        . ~ . + kernel_connectivity_std # flexible covariate inclusion
      ),
      data    = data_dir_sac[[nm]],
      family  = binomial
    )
  }
)

names(sac_envpa_models) <- names(pa_glm_models) 


# A5) Index models
model_sets <- list(
  null_models,
  paonly_models,
  envpa_models,
  sac_envpa_models
)


# C) Model comparison
### fit and outputs (coefficients, etc.) among null, pa-only, env-pa, and sac-env-pa models

mcfadden_r2 <- function(model, null_model) {
  
  ll_model <- as.numeric(logLik(model))
  ll_null  <- as.numeric(logLik(null_model))
  
  1 - (ll_model / ll_null)
}


model_set_comparison <- lapply(
  names(pa_glm_models),
  function(nm) {
    
    # model object definitions
    null_mod  <- null_models[[nm]]
    paonly_mod    <- paonly_models[[nm]]
    envpa_mod  <- envpa_models[[nm]]
    sac_envpa_mod   <- sac_envpa_models[[nm]]
    
    models = list( # labels
      Null_Model = null_mod,
      PAOnly_Model = paonly_mod,
      EnvPA_Model = envpa_mod,
      SAC_EnvPA_Model = sac_envpa_mod
    )
    
   # fit-summaries per model object 
    model_details <- lapply(models, function(mod) {
      list(
        fit     = mod,
        summary = summary(mod),
        glance  = broom::glance(mod),
        tidy    = broom::tidy(mod)
      )
    })
    
    # helper: mcfadden r2
    mcfadden_wrapper <- function(mod, null) {
      tryCatch(
        mcfadden_r2(mod, null),
        error = function(e) NA_real_
      )
    }

    
   comparison <- data.frame(
      model = names(models),
      
      AICc = sapply(models, MuMIn::AICc),
      
      deviance_explained = c(
        1 - (null_mod$deviance / null_mod$null.deviance),
        1 - (paonly_mod$deviance  / paonly_mod$null.deviance),
        1 - (envpa_mod$deviance / envpa_mod$null.deviance),
        1 - (sac_envpa_mod$deviance / sac_envpa_mod$null.deviance)
      ),
      
      mcfadden_R2v_null = c(
        0,
        mcfadden_wrapper(paonly_mod, null_mod),
        mcfadden_wrapper(envpa_mod, null_mod),
        mcfadden_wrapper(sac_envpa_mod, null_mod)
      )
    )
    
    # deltas AICc, R2 relative to envpa_mod
    comparison$delta_AICc <- comparison$AICc - comparison$AICc[comparison$model == "EnvPA_Model"]
    comparison$delta_mcfadden_R2 <- comparison$mcfadden_R2 - comparison$mcfadden_R2[comparison$model == "EnvPA_Model"]
    
    list(
      models = model_details,
      comparison_table = comparison
    )
  }
)
    

names(model_set_comparison) <- names(pa_glm_models)
model_set_comparison



# D) Visualize 
# pivot to long df

# Figure 1: Comparison metrics

comparison_df <- bind_rows(
  lapply(model_set_comparison, function(x) {
    x$comparison_table
  }),
  .id = "response"
)


comparison_long <- comparison_df %>%
  tidyr::pivot_longer(
    cols = c(delta_AICc, deviance_explained, mcfadden_R2v_null),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = dplyr::recode(metric,
                           delta_AICc = "Δ AICc",
                           deviance_explained = "Deviance Explained",
                           mcfadden_R2v_null = "McFadden R²"
    )
  )


label_size <- 4
label_fontface <- "bold"

ggplot(comparison_long, aes(x = model, y = value)) +
  
  geom_col(aes(fill = model), show.legend = FALSE) +
  
  geom_text(
    aes(label = round(value, 3)),
    position = position_stack(vjust = 0.5),  # <- centers text in bar
    size = label_size,
    fontface = label_fontface,
    color = "black"  # improves contrast inside bars
  ) +
  
  facet_grid(metric ~ response, scales = "free_y") +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 10)
  ) +
  
  labs(
    x = NULL,
    y = "Value"
  )




# Figure 2:  PA response curves
MakeCurve <- function(model,
                      focal_var,
                      data,
                      model_name = NULL,
                      response_name = NULL,
                      n = 100) {
  
  mf <- model.frame(model)
  vars <- all.vars(formula(model))[-1]
  
  base_vals <- setNames(
    lapply(vars, function(v) {
      
      if (v == focal_var) return(NULL)
      
      if (is.numeric(mf[[v]])) {
        median(mf[[v]], na.rm = TRUE)
      } else {
        names(which.max(table(mf[[v]])))
      }
    }),
    vars
  )
  
  base_vals <- base_vals[!sapply(base_vals, is.null)]

  
  nd <- as.data.frame(base_vals)
  nd <- nd[rep(1, n), ]
  
  nd[[focal_var]] <- seq(
    min(mf[[focal_var]], na.rm = TRUE),
    max(mf[[focal_var]], na.rm = TRUE),
    length.out = n
  )

  pred <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
  
  tibble::tibble(
    x = nd[[focal_var]],
    fit = plogis(pred$fit),
    lwr = plogis(pred$fit - 1.96 * pred$se.fit),
    upr = plogis(pred$fit + 1.96 * pred$se.fit),
    model = model_name,
    response = response_name
  )
}



pa_curve_df <- bind_rows(
  
  # Colonization
  MakeCurve(
    envpa_models$RLL_col,
    focal_var = "pa_prop",
    data = data_dir$RLL_col,
    model_name = "Env + PA",
    response_name = "Colonization"
  ),
  
  MakeCurve(
    sac_envpa_models$RLL_col,
    focal_var = "pa_prop",
    data = data_dir$RLL_col,
    model_name = "SAC + Env + PA",
    response_name = "Colonization"
  ),
  
  # Extinction
  MakeCurve(
    envpa_models$RLL_ext,
    focal_var = "pa_prop",
    data = data_dir$RLL_ext,
    model_name = "Env + PA",
    response_name = "Extinction"
  ),
  
  MakeCurve(
    sac_envpa_models$RLL_ext,
    focal_var = "pa_prop",
    data = data_dir$RLL_ext,
    model_name = "SAC + Env + PA",
    response_name = "Extinction"
  )
)



ggplot(pa_curve_df, aes(x = x, y = fit, color = model, fill = model)) +
  
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  
  geom_line(linewidth = 1) +
  
  facet_wrap(~response) +
  
  scale_x_continuous(labels = scales::percent) +
  
  labs(
    x = "Protected Area",
    y = "Predicted Probability",
    title = "Marginal PA Response (Comparable Model Structures Only)"
  ) +
  
  theme_minimal()




# Figure 3: SAC response curves
sac_term_col <- "kernel_connectivity_std"
sac_term_ext <- "kernel_connectivity_std"



sac_curve_df <- bind_rows(
  
  MakeCurve(
    sac_envpa_models$RLL_col,
    sac_term_col,
    data_dir$RLL_col
  ) %>% mutate(
    response="Colonization"
  ),
  
  MakeCurve(
    sac_envpa_models$RLL_ext,
    sac_term_ext,
    data_dir$RLL_ext
  ) %>% mutate(
    response="Extinction"
  )
)


ggplot(
  sac_curve_df,
  aes(
    x = x,
    y = fit
  )
) +
  
  geom_ribbon(
    aes(
      ymin=lwr,
      ymax=upr
    ),
    alpha=.2
  ) +
  
  geom_line(
    linewidth=1
  ) +
  
  facet_wrap(
    ~response,
    scales="free_x"
  ) +
  
  theme_minimal() +
  
  labs(
    title="Spatial Autocorrelation Effect",
    x="SAC Covariate",
    y="Predicted Probability"
  )



# Figure 4: PA, SAC coefficient comparison
coef_df <- bind_rows(
  
  broom::tidy(paonly_models$RLL_col,
              conf.int=TRUE) %>%
    mutate(
      model="PAOnly",
      response="Colonization"
    ),
  
  broom::tidy(envpa_models$RLL_col,
              conf.int=TRUE) %>%
    mutate(
      model="EnvPA",
      response="Colonization"
    ),
  
  broom::tidy(sac_envpa_models$RLL_col,
              conf.int=TRUE) %>%
    mutate(
      model="SAC_EnvPA",
      response="Colonization"
    ),
  
  broom::tidy(paonly_models$RLL_ext,
              conf.int=TRUE) %>%
    mutate(
      model="PAOnly",
      response="Extinction"
    ),
  
  broom::tidy(envpa_models$RLL_ext,
              conf.int=TRUE) %>%
    mutate(
      model="EnvPA",
      response="Extinction"
    ),
  
  broom::tidy(sac_envpa_models$RLL_ext,
              conf.int=TRUE) %>%
    mutate(
      model="SAC_EnvPA",
      response="Extinction"
    )
)


coef_df <- coef_df %>%
  filter(
    term %in% c(
      "pa_prop",
      sac_term_col
    )
  )



ggplot(
  coef_df,
  aes(
    x = estimate,
    y = model,
    color = model
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = .2
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_grid(
    response ~ term
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Coefficient",
    y = NULL,
    title = "PA and SAC Effect Sizes"
  )





#########################################
### --- RESIDUAL SAC DIAGNOSTICS --- ###
########################################

### Evaluate whether adding biologically informed SAC covariates reduces
# spatial structure in model residuals.
### SAC Metrics: 
# Moran's I: global SAC diagnostic; assess similarity of spatially adjacent residuals
# Semivariance (-variograms): quantifies distance-decay (dissimilarity)structure of SAC across various scales
### Residual Types:
# DHARMa Residuals: standardized [U(0,1]), rank-based/quantile
# Pearson Residuals: sd units, approx. N(mean, sd)


### Workflow:
# 1) Fit nested model hierarchy
### Null, PA-only, Full-PA, SAC models
# 2) Extract residuals
### Pearson, DHARMa
# 3) Quantify residual spatial structure
### Define spatial weighting schemes (Moran's I)
### Morans' I, Semivariance/Semivariograms
# 4) Assess whether SAC term modulates residual structure
# 5) If needed, autoregressive/spatial model approach


# Indexed models
#model_sets <- list(
#  null_models,
#  paonly_models,
#  envpa_models,
#  sac_envpa_models
#)


### --- SPATIAL STRUCTURES --- ###

### A) Define spatial coordinate structures

# A1) Full block coords set (n = 3337)
coords_lookup <- spp_blocks_sf %>%
  st_drop_geometry() %>%
  dplyr::select(
    atlas_block,
    x_km,
    y_km
  )

# A2) Response-specific coords sets
# Must partition block coord pool, as each model runs on block subset

coords_by_model <- lapply( # partition blocks
  names(pa_glm_models),
  function(nm) {
    
    blocks <- data_dir[[nm]]$atlas_block
    
    coords_lookup %>%
      filter(atlas_block %in% blocks) %>%
      right_join(
        tibble(atlas_block = blocks),
        by = "atlas_block"
      ) %>%
      arrange(match(atlas_block, blocks))
  })

names(coords_by_model) <- names(pa_glm_models)


distance_matrices <- lapply( # response-specific centroid matrices
  coords_by_model,
  function(coords) {
    
    d <- as.matrix(
      dist(
        coords[, c("x_km", "y_km")]
      )
    )
    
    diag(d) <- NA_real_
    d
  }
)


### B) Response-specific neighborhood structures
# within max dispersal distance
neighbor_matrices <- lapply(
  distance_matrices,
  function(dmat) {
    
    !is.na(dmat) & dmat <= max_disp_km
  }
)


### C) Response-specific spatial weighting schemes
# Necessary for Moran's I; After accounting for environment/climate/etc. in model, 
# do residuals remain spatially autocorrelated AMONG biologically relevant neighborhood, 
# i.e. DO use max dispersal radius, do NOT limit neighbors those occupied in A1

### C1) Binary adjacency weights
# all neighbors influence equally
binary_weights <- lapply(names(neighbor_matrices), function(nm) {
  
  mat <- neighbor_matrices[[nm]] * 1  # convert TRUE/FALSE → 1/0
  
  spdep::mat2listw(
    mat,
    style = "W",
    zero.policy = TRUE
  )
})

names(binary_weights) <- names(neighbor_matrices)


### C2) IDW spatial weights 
# nearby neighbors influence more than faraway neighbors
idw_weights <- lapply(names(distance_matrices), function(nm) {
  
  dmat <- distance_matrices[[nm]]
  
  idw <- ifelse(
    dmat > 0 & !is.na(dmat),
    1 / dmat,
    0
  )
  
  idw[!neighbor_matrices[[nm]]] <- 0
  
  spdep::mat2listw(
    idw,
    style = "W",
    zero.policy = TRUE
  )
})

names(idw_weights) <- names(distance_matrices)


### C3) Gaussian kernel weights
# assessed residual SAC for for source- and non-source neighbors
# above GKW uses A1-occupancy to create static map of neighbor-summed influence on a given 
# block, then creating covariate to model
# from residual assessment, need to understand all existing spatial clustering 
# (though still confined within max dispersal distance), regardless of previous occupancy
gaussian_weights <- lapply(names(distance_matrices), function(nm) {
  
  dmat <- distance_matrices[[nm]]
  
  k <- exp(-(dmat^2) / (2 * sigma_km^2))
  k[!neighbor_matrices[[nm]]] <- 0
  
  row_sums <- rowSums(k, na.rm = TRUE)
  
  # preserve matrix structure explicitly
  for (i in seq_len(nrow(k))) {
    if (row_sums[i] > 0) {
      k[i, ] <- k[i, ] / row_sums[i]
    } else {
      k[i, ] <- 0
    }
  }
  
  k[is.na(k)] <- 0
  
  spdep::mat2listw(
    k,
    style = "W",
    zero.policy = TRUE
  )
})

names(gaussian_weights) <- names(distance_matrices)



### C4) Index weight sets
weight_sets <- lapply(names(pa_glm_models), function(nm) {
  
  list(
    binary   = binary_weights[[nm]],
    idw      = idw_weights[[nm]],
    gaussian = gaussian_weights[[nm]]
  )
})

names(weight_sets) <- names(pa_glm_models)




### --- EXTRACT MODEL RESIDUALS --- ###

### A) DHARMa residuals
set.seed(123)

extract_dharma <- function(models, model_type) {
  
  stopifnot(length(models) > 0)
  
  bind_rows(lapply(names(models), function(nm) {
    
    mod <- models[[nm]]
    sim <- DHARMa::simulateResiduals(mod, plot = FALSE)
    
    tibble(
      atlas_block   = data_dir[[nm]]$atlas_block,
      residual      = sim$scaledResiduals,
      model_name    = nm,
      model_type    = model_type,
      residual_type = "DHARMa"
    )
  }))
}


### B) Pearson residuals
extract_pearson <- function(models, model_type) {
  
  bind_rows(lapply(names(models), function(nm) {
    
    mod <- models[[nm]]
    
    tibble(
      atlas_block   = data_dir[[nm]]$atlas_block,  # ← same fix here
      residual      = residuals(mod, type = "pearson"),
      model_name    = nm,
      model_type    = model_type,
      residual_type = "Pearson"
    )
  }))
}


### C) Combine residuals
model_sets <- list(
  null_mod      = null_models,
  paonly_mod    = paonly_models,
  envpa_mod     = envpa_models,
  sac_envpa_mod = sac_envpa_models
)


dharma_residuals <- bind_rows(
  extract_dharma(model_sets$null_mod,     "true_null"),
  extract_dharma(model_sets$paonly_mod,       "pa_only"),
  extract_dharma(model_sets$envpa_mod,       "env_pa"),
  extract_dharma(model_sets$sac_envpa_mod, "local_sac")
)

pearson_residuals <- bind_rows(
  extract_pearson(model_sets$null_mod, "true_null"),
  extract_pearson(model_sets$paonly_mod, "pa_only"),
  extract_pearson(model_sets$envpa_mod, "env_pa"),
  extract_pearson(model_sets$sac_envpa_mod, "local_sac")
)


all_residuals <- bind_rows(
  dharma_residuals,
  pearson_residuals
) %>%
  left_join(
    coords_lookup,
    by = "atlas_block"
  )


all_residuals <- all_residuals %>%
  mutate(
    model_type = case_when(
      model_type %in% c("true_null") ~ "true_null",
      model_type %in% c("pa_only") ~ "pa_only",
      model_type %in% c("env_pa") ~ "env_pa",
      model_type %in% c("local_sac", "SAC_EnvPA_Model", "sac_envpa") ~ "local_sac",
      TRUE ~ model_type
    )
  )



### --- TEST MODEL RESIDUALS --- ###

### --- A) Moran's I
# MI+: observed nearby residuals are more similar than expected;
# MI~0: little reamining SAC;
# significant p-values: residual SAC present

model_names <- names(weight_sets)


# A1) Analyze (3 weighting schemes)
morans_results <- bind_rows(
  
  lapply(model_names, function(mn) {
    
    bind_rows(
      
      lapply(unique(all_residuals$model_type), function(mt) {
        
        bind_rows(
          
          lapply(unique(all_residuals$residual_type), function(rt) {
            
            dat <- all_residuals %>%
              filter(
                model_name == mn,
                model_type == mt,
                residual_type == rt
              )
            
            coords <- coords_by_model[[mn]]
            
            # --- strict alignment ---
            dat <- dat %>%
              right_join(
                tibble(atlas_block = coords$atlas_block),
                by = "atlas_block"
              ) %>%
              arrange(match(atlas_block, coords$atlas_block)) %>%
              filter(!is.na(residual))
            
            stopifnot(length(dat$residual) == length(coords$atlas_block))
            
            bind_rows(
              lapply(names(weight_sets[[mn]]), function(wname) {
                
                w <- weight_sets[[mn]][[wname]]
                
                mi <- spdep::moran.test(dat$residual, w, zero.policy = TRUE)
                
                tibble(
                  model_name    = mn,
                  model_type    = mt,
                  residual_type = rt,
                  weight_type   = wname,
                  morans_I      = unname(mi$estimate[1]),
                  expected_I    = unname(mi$estimate[2]),
                  variance_I    = unname(mi$estimate[3]),
                  p_value       = mi$p.value
                )
              })
            )
          })
        )
      })
    )
  })
)


# A2) Quantitative results

morans_results <- morans_results %>%
  
  mutate(
    process = case_when(
      grepl("_col$", model_name) ~ "Colonization",
      grepl("_ext$", model_name) ~ "Extirpation",
      TRUE ~ "All transitions"
    ),
    
    model_type = factor(
      model_type,
      levels = c("true_null", "pa_only", "env_pa", "local_sac"),
      labels = c(
        "Null Model",
        "PA-only Model",
        "Environment-PA Model",
        "Environment-PA w SAC Model"
      )
    )
  )


# A3) Delta Moran's I (correct baseline)

null_baseline <- morans_results %>%
  filter(model_type == "Null Model") %>%
  dplyr::select(weight_type, residual_type, process, null_I = morans_I)


delta_morans_results <- morans_results %>%
  left_join(
    null_baseline,
    by = c("weight_type", "residual_type", "process")
  ) %>%
  mutate(delta_morans_I = morans_I - null_I)


delta_summary <- delta_morans_results %>%
  
  group_by(weight_type, residual_type, model_type) %>%
  
  summarise(
    mean_delta_morans_I = mean(delta_morans_I, na.rm = TRUE),
    sd_delta_morans_I   = sd(delta_morans_I, na.rm = TRUE),
    min_delta           = min(delta_morans_I, na.rm = TRUE),
    max_delta           = max(delta_morans_I, na.rm = TRUE),
    .groups = "drop"
  )


# A4) Visualize

# Raw Moran's I
ggplot(
  morans_results,
  aes(
    x = model_type,
    y = morans_I,
    fill = residual_type
  )
) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_grid(process ~ weight_type) +
  labs(
    title = "Residual Moran's I Across Models and Weight Structures",
    x = NULL,
    y = "Moran's I",
    fill = "Residual Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )


# Delta Moran's I
ggplot(
  delta_morans_results,
  aes(
    x = model_type,
    y = delta_morans_I,
    fill = residual_type
  )
) +
  geom_col(position = position_dodge()) +
  facet_grid(process ~ weight_type) +
  labs(
    title = "Reduction in Spatial Autocorrelation (Δ Moran's I vs Null)",
    x = NULL,
    y = "Δ Moran's I",
    fill = "Residual type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )





### --- B) Semivariance, Semivariograms
# Quantify how residual similarity changes with distance
# Plot interp: low semivariance at short distances = nearby blocks w similar residuals
# = evidence of positive SAC;
# increasing semivariance w distance = spatial structure decays across space;
# Flat semivariogram = little remaining SAC

# B1) Analysis
### --- B) Semivariance, Semivariograms --- ###
# Quantify spatial structure in residuals across distance classes

semivar_df <- all_residuals %>%
  
  group_by(residual_type, model_type, model_name) %>%
  
  group_modify(~ {
    
    dat <- .x
    
    # safety check
    if (nrow(dat) < 5) return(tibble())
    
    gdf <- sf::st_as_sf(
      dat,
      coords = c("x_km", "y_km"),
      remove = FALSE
    )
    
    cutoff_i <- 25  # km ecological scale
    
    v <- gstat::variogram(
      residual ~ 1,
      locations = gdf,
      cutoff = cutoff_i,
      width  = 5
    )
    
    as.data.frame(v)
  }) %>%
  
  ungroup() %>%
  
  filter(!is.na(gamma)) %>%
  
  mutate(
    
    process = case_when(
      grepl("_col$", model_name) ~ "Colonization",
      grepl("_ext$", model_name) ~ "Extirpation",
      TRUE ~ "All transitions"
    ),
    
    # IMPORTANT: assumes you already standardized model_type upstream
    model_type = factor(
      model_type,
      levels = c("true_null", "pa_only", "env_pa", "local_sac"),
      labels = c(
        "Null Model",
        "PA-only Model",
        "Env-PA Model",
        "SAC-Env-PA Model"
      )
    )
  )


# B2b) Summarize
semivar_summary <- semivar_df %>%
  
  group_by(residual_type, process, model_type) %>%
  
  summarise(
    mean_semivariance = mean(gamma, na.rm = TRUE),
    max_semivariance  = max(gamma, na.rm = TRUE),
    min_semivariance  = min(gamma, na.rm = TRUE),
    semivariance_range = max_semivariance - min_semivariance,
    .groups = "drop"
  )

semivar_summary


# B2c) Visualize
# Full distance (see above, cutoff = )
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
  geom_point(size = 1.5, na.rm = TRUE) +
  
  facet_grid(residual_type ~ process) +
  
  scale_color_manual(values = c(
    "Null Model" = "grey40",
    "PA-only Model" = "orange",
    "Env-PA Model" = "turquoise4",
    "SAC-Env-PA Model" = "darkorchid"
  )) +
  
  labs(
    title = "Residual Spatial Structure Across Model Types",
    x = "Distance Between Blocks (km)",
    y = "Semivariance",
    color = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")



# Zoom/Clipped distance
ggplot(
  semivar_df %>% filter(dist <= 25),
  aes(
    x = dist,
    y = gamma,
    color = model_type,
    group = interaction(model_type, model_name, process)
  )
) +
  
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  
  facet_wrap(~process) +
  
  scale_color_manual(values = c(
    "Null Model" = "grey40",
    "PA-only Model" = "orange",
    "Env-PA Model" = "turquoise4",
    "SAC-Env-PA Model" = "darkorchid"
  )) +
  
  labs(
    title = "Residual Spatial Structure (0–25 km)",
    x = "Distance (km)",
    y = "Semivariance",
    color = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")




# B3) Summarize pairwise binning

### np = number of pairwise comparisons contributing to each
# semivariance estimate.
# Important because unstable bins with few comparisons can produce
# noisy semivariance estimates.


# B3a) Analyze
pairs_summary <- semivar_df %>%
  
  filter(model_type == "Env-PA Model") %>%
  
  group_by(process) %>%
  
  summarise(
    mean_pairs = mean(np, na.rm = TRUE),
    min_pairs  = min(np, na.rm = TRUE),
    max_pairs  = max(np, na.rm = TRUE),
    sd_pairs   = sd(np, na.rm = TRUE),
    .groups = "drop"
  )

pairs_summary <- semivar_df %>%
  
  filter(model_type == "Env-PA Model") %>%
  
  group_by(process) %>%
  
  summarise(
    mean_pairs = mean(np, na.rm = TRUE),
    min_pairs  = min(np, na.rm = TRUE),
    max_pairs  = max(np, na.rm = TRUE),
    sd_pairs   = sd(np, na.rm = TRUE),
    .groups = "drop"
  )

pairs_summary


pair_bins_summary <- semivar_df %>%
  
  filter(model_type == "Env-PA Model") %>%
  
  dplyr::select(process, dist, np) %>%
  distinct() %>%
  arrange(process, dist)


ggplot(pair_bins_summary,
       aes(x = dist, y = np)) +
  
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  
  facet_wrap(~process, scales = "free_y") +
  
  labs(
    title = "Pairwise Comparisons per Distance Bin",
    x = "Distance Bin (km)",
    y = "Number of Pairs"
  ) +
  
  theme_minimal(base_size = 13)




#### Overall SAC Reduction ###
# SACRI: SAC reduction index; SACRI = MoransI from fitted model / MoransI from null model
# ---> how much residual SAC did a covariate term remove?

baseline_I <- morans_results %>%
  filter(model_type == "Null Model") %>%
  group_by(residual_type, weight_type) %>%
  summarise(I0 = mean(morans_I, na.rm = TRUE), .groups = "drop")

srei_results <- morans_results %>%
  
  left_join(baseline_I,
            by = c("residual_type", "weight_type")) %>%
  
  mutate(
    SREI = (I0 - morans_I) / (abs(I0) + 1e-6)
  )

srei_summary <- srei_results %>%
  group_by(model_type, weight_type) %>%
  summarise(mean_SREI = mean(SREI, na.rm = TRUE), .groups = "drop")

srei_summary



# Plot
ggplot(srei_summary,
       aes(x = model_type, y = mean_SREI, fill = weight_type)) +
  
  geom_col(position = position_dodge(width = 0.8)) +
  
  labs(
    title = "Spatial Autocorrelation Reduction Efficiency Index (SREI)",
    y = "SREI (1 - I_model / I_null)",
    x = NULL
  ) +
  
  theme_minimal(base_size = 13)
