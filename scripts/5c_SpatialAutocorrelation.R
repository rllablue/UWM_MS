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

# Visualization
library(ggplot2)
library(viridis)
library(gt)
library(webshot2)



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






### SPATIAL AUTOCORRELATION ###

# Goals:
# 1) Fit null + full PA models
# 2) Extract DHARMa standardized quantile residuals
# 3) Build semivariograms from DHARMa residuals


### --- TRUE NULL, PA NULL, FULL MODEL COMPARISON --- ###

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






