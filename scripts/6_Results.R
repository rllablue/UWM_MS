### Store, compare, visualize results 

library(caret)


### Assess all fitted models per response

cand_model_sets <- list(
  
  ENV                   = global_glm_models,
  
  ENV_PA                = pa_glm_models,
  ENV_PA_SAC            = sac_glm_models,
  
  ENV_PA_GAP            = pasub_glm_models,
  ENV_PA_GAP_SAC        = sac2_glm_models,
  
  INT_ENV_PA_SAC        = int_glm_models,
  INT_ENV_PA_GAP_SAC    = int2_glm_models
  
)

# Model Fit Comparison
CompareModels <- function(response_name, model_sets){
  
  mods <- lapply(model_sets, `[[`, response_name)
  
  tab <- data.frame(
    Model  = names(mods),
    K      = sapply(mods, function(x) attr(logLik(x), "df")),
    LogLik = sapply(mods, logLik),
    AICc   = sapply(mods, AICc),
    stringsAsFactors = FALSE
  )
  
  # AICc
  tab <- tab[order(tab$AICc), ]
  tab$Delta_AICc <- tab$AICc - min(tab$AICc)
  relLik <- exp(-0.5 * tab$Delta_AICc)
  tab$Weight <- relLik / sum(relLik)
  tab$Deviance_Explained <- sapply(
    mods[tab$Model],
    function(m){
      1 - (m$deviance / m$null.deviance)
    }
  )
  
  # McFadden pseudo-R2
  tab$McFadden_R2 <- sapply(
    mods[tab$Model],
    function(m){
      pscl::pR2(m)[["McFadden"]]
    }
  )
  
  # AUC
  tab$AUC <- sapply(
    mods[tab$Model],
    function(m){
      
      y <- model.response(model.frame(m))
      p <- fitted(m)
      
      as.numeric(
        pROC::auc(
          pROC::roc(
            response = y,
            predictor = p,
            quiet = TRUE
          )
        )
      )
    }
  )
  
  # Brier score
  tab$Brier <- sapply(
    mods[tab$Model],
    function(m){
      
      y <- model.response(model.frame(m))
      p <- fitted(m)
      
      mean((y - p)^2)
    }
  )
  
  
  rownames(tab) <- NULL
  
  tab
}


col_comp <- CompareModels("RLL_col", cand_model_sets)
ext_comp <- CompareModels("RLL_ext", cand_model_sets)

col_comp
ext_comp

# top candidates, i.e. deltaAICc < 2
top_col_models <- col_comp %>%
  filter(Delta_AICc <= 2) %>%
  pull(Model)

top_ext_models <- ext_comp %>%
  filter(Delta_AICc <= 2) %>%
  pull(Model)

top_col_models
top_ext_models


# Coefficient Comparison
AssessTopModels <- function(response_name,
                            supported_models,
                            model_sets){
  
  mods <- lapply(
    supported_models,
    function(x) model_sets[[x]][[response_name]]
  )
  
  names(mods) <- supported_models
  
  mods
}


top_col <- AssessTopModels(
  "RLL_col",
  top_col_models,
  cand_model_sets
)

top_ext <- AssessTopModels(
  "RLL_ext",
  top_ext_models,
  cand_model_sets
)


#############################
### Top Model Diagnostic Comparisons ###
#############################

top_models <- c(
  top_col,
  top_ext
)

models <- top_models


### --- 1A: Coefficient Summary

coefs_summary <- bind_rows(
  lapply(names(models), function(nm){
    
    broom::tidy(
      models[[nm]],
      conf.int = TRUE
    ) %>%
      mutate(model = nm)
    
  })
)

coefs_summary


### --- 1B: Potentially Uninformative Parameters

uninf_pars <- coefs_summary %>%
  
  filter(term != "(Intercept)") %>%
  
  mutate(
    
    ci_overlaps_zero =
      conf.low <= 0 &
      conf.high >= 0,
    
    non_significant =
      p.value > 0.05,
    
    potentially_uninformative =
      ci_overlaps_zero &
      non_significant
    
  )

# summarize
uninf_summary <- uninf_pars %>%
  
  group_by(model) %>%
  
  summarise(
    
    n_predictors = n(),
    
    n_uninformative =
      sum(
        potentially_uninformative,
        na.rm = TRUE
      ),
    
    pct_uninformative =
      100 *
      n_uninformative /
      n_predictors,
    
    .groups = "drop"
    
  )

uninf_summary


### --- 1C: Likelihood Ratio Tests

lrt_results <- bind_rows(
  lapply(names(models), function(nm){
    
    drop1(
      models[[nm]],
      test = "Chisq"
    ) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term") %>%
      mutate(model = nm)
    
  })
)

lrt_results


### --- 1D: Calibration
calibration_df <- bind_rows(
  
  lapply(names(models), function(nm){
    
    y <- model.response(
      model.frame(models[[nm]])
    )
    
    p <- fitted(models[[nm]])
    
    tibble(
      observed = y,
      predicted = p,
      model = nm
    )
    
  })
  
)

# bin predictions
calibration_summary <- calibration_df %>%
  
  mutate(
    decile =
      ntile(predicted, 10)
  ) %>%
  
  group_by(
    model,
    decile
  ) %>%
  
  summarise(
    
    mean_pred =
      mean(predicted),
    
    obs_rate =
      mean(observed),
    
    n = n(),
    
    .groups = "drop"
    
  )

ggplot(
  calibration_summary,
  aes(
    mean_pred,
    obs_rate
  )
) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  
  geom_point(size = 3) +
  
  geom_line() +
  
  facet_wrap(~model) +
  
  labs(
    x = "Predicted Probability",
    y = "Observed Frequency",
    title = "Calibration Curves"
  ) +
  
  theme_minimal()

calibration_summary



### --- 1D: DHARMa Residual Diagnostics

set.seed(123)

resids_dharma <- lapply(
  models,
  DHARMa::simulateResiduals,
  plot = FALSE
)

# Uniformity
resids_uniformity <- bind_rows(
  
  lapply(
    names(resids_dharma),
    function(nm){
      
      tst <- DHARMa::testUniformity(
        resids_dharma[[nm]],
        plot = FALSE
      )
      
      tibble(
        model = nm,
        statistic = tst$statistic,
        p_value = tst$p.value
      )
    })
)

# Dispersion
resids_dispersion <- bind_rows(
  
  lapply(
    names(resids_dharma),
    function(nm){
      
      tst <- DHARMa::testDispersion(
        resids_dharma[[nm]],
        plot = FALSE
      )
      
      tibble(
        model = nm,
        statistic = tst$statistic,
        p_value = tst$p.value
      )
    })
)

# Outliers
resids_outlier <- bind_rows(
  
  lapply(
    names(resids_dharma),
    function(nm){
      
      tst <- DHARMa::testOutliers(
        resids_dharma[[nm]],
        plot = FALSE
      )
      
      tibble(
        model = nm,
        p_value = tst$p.value
      )

    })
)

resids_uniformity
resids_dispersion
resids_outlier




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





































### --- 6: STORE EFFECTS --- ###

# STORING #
pa_results_df <- data.frame()


# Extract, Collect Base PA effects 
ExtractPACoefficients <- function(model, model_name, spp_name) {
  
  parts <- strsplit(model_name, "_")[[1]]
  block    <- parts[1]
  response <- parts[2]
  
  df <- broom::tidy(model)
  
  if (!"pa_prop" %in% df$term) {
    return(NULL)
  }
  
  df %>%
    filter(term == "pa_prop") %>%
    mutate(
      species   = spp_name,
      response  = response,
      block     = block,
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      sig       = ifelse(conf.low > 0 | conf.high < 0, "yes", "no")
    ) %>%
    dplyr::select(species, block, response, estimate, std.error, conf.low, conf.high, p.value, sig)
}


species_results <- dplyr::bind_rows(
  lapply(names(pa_glm_models), function(nm) {
    ExtractPACoefficients(pa_glm_models[[nm]], nm, spp_name)
  })
)


if (!"species" %in% names(rll_results_df)) {
  
  rll_results_df <- species_results
  
} else {
  
  rll_results_df <- rll_results_df %>%
    dplyr::filter(species != spp_name) %>%
    dplyr::bind_rows(species_results)
  
}

