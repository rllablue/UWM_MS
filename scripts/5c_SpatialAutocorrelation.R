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
library(spdep)

# Visualization
library(ggplot2)
library(viridis)
library(gt)
library(webshot2)




########### WIP 4/26/26


### MORAN'S I ###
pa_moran <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir[[nm]]
  res <- pa_residuals[[nm]]
  
  coords <- cbind(dat$lon, dat$lat)
  
  # k-nearest neighbors (k = 8 is a common default for grids)
  knn <- spdep::knearneigh(coords, k = 8)
  nb  <- spdep::knn2nb(knn)
  lw  <- spdep::nb2listw(nb, style = "W")
  
  test <- spdep::moran.test(res, lw)
  
  data.frame(
    model = nm,
    Moran_I = as.numeric(test$estimate[1]),
    p_value = test$p.value
  )
})

pa_moran <- do.call(rbind, pa_moran)
pa_moran



### CORRELOGRAM
pa_correlog <- lapply(names(pa_glm_models), function(nm) {
  
  dat <- data_dir[[nm]]
  res <- pa_residuals[[nm]]
  
  correlog(
    x = dat$lon,
    y = dat$lat,
    z = res,
    increment = 0.5,  # adjust if needed for your grid spacing
    resamp = 0
  )
})

names(pa_correlog) <- names(pa_glm_models)


correlog_df <- do.call(rbind, lapply(names(pa_correlog), function(nm) {
  
  obj <- pa_correlog[[nm]]
  
  data.frame(
    model    = nm,
    distance = obj$mean.of.class,
    moran_I  = obj$correlation
  )
}))




ggplot(correlog_df, aes(x = distance, y = moran_I, color = model)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  scale_color_manual(
    values = c(
      "RLL_col" = "#2C7BB6",
      "RLL_ext" = "#D7191C"
    ),
    labels = c(
      "Colonization",
      "Extinction"
    )
  ) +
  
  labs(
    title = "Spatial Correlogram of PA Model Residuals",
    x = "Distance",
    y = "Moran's I (by distance class)",
    color = NULLe
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "bottom",
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





# Compare to non-spatial model
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

library(pROC)

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
library(pscl)

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












# Models per response x block set
global_glm_models
# $RLL_col, $RLL_ext, $DNR_col, $DNR_ext



# Atlas blocks geometry
blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # load full block map
  rename(atlas_block = BLOCK_ID)
crs(blocks_all_sf)

blocks_background <- blocks_all_sf %>%
  filter(atlas_block %in% union(blocks_rll, blocks_dnr))




### MORAN'S I ###

RunMoranResiduals <- function(model,
                              model_df,
                              blocks_sf,
                              blocks_rll,
                              blocks_dnr,
                              model_name,
                              k = 4) {
  
  # Determine pool
  pool_blocks <- if (grepl("RLL", model_name)) {
    blocks_rll
  } else {
    blocks_dnr
  }
  
  background_sf <- blocks_sf %>%
    filter(atlas_block %in% pool_blocks)
  
  sf_obj <- model_df %>%
    left_join(
      blocks_sf %>% dplyr::select(atlas_block, geometry),
      by = "atlas_block"
    ) %>%
    st_as_sf() %>%
    arrange(atlas_block)
  
  sf_obj$resid <- residuals(model, type = "pearson")
  
  coords <- st_coordinates(st_centroid(sf_obj))
  nb <- knn2nb(knearneigh(coords, k = k))
  lw <- nb2listw(nb, style = "W")
  
  moran_res <- moran.test(sf_obj$resid, lw)
  
  # Labels
  type_label <- ifelse(grepl("col", model_name),
                       "Colonization", "Extinction")
  set_label  <- ifelse(grepl("RLL", model_name),
                       "RLL", "DNR")
  
  cat("\n============================\n")
  cat(set_label, "-", type_label, "\n")
  print(moran_res)
  
  # Moran plot
  moran.plot(sf_obj$resid, lw,
             main = paste(set_label, "-", type_label))
  
  # Residual map
  print(
    ggplot() +
      geom_sf(data = background_sf,
              fill = NA,
              color = "grey70") +
      geom_sf(data = sf_obj,
              aes(fill = resid)) +
      scale_fill_gradient2(midpoint = 0) +
      ggtitle(paste(set_label, "-", type_label, "Residuals")) +
      theme_minimal()
  )
  
  return(moran_res)
}


moran_results <- lapply(names(global_glm_models), function(nm) {
  
  RunMoranResiduals(
    model        = global_glm_models[[nm]],
    model_df     = data_dir[[nm]],
    blocks_sf    = blocks_all_sf,
    blocks_rll   = blocks_rll,
    blocks_dnr   = blocks_dnr,
    model_name   = nm
  )
})

names(moran_results) <- names(global_glm_models)
moran_results