if (!require("pacman")) install.packages("pacman")
pacman::p_load(psych, usdm)


# Core
library(dplyr)
library(tidyverse)
library(magrittr)
library(purrr)
library(stringr)
library(readxl)
library(ncf)

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

# Visualization
library(ggplot2)
library(viridis)
library(gt)
library(webshot2)




# Models per response x block set
global_glm_models
# $RLL_col, $RLL_ext, $DNR_col, $DNR_ext



# Atlas blocks geometry
blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # load full block map
  rename(atlas_block = BLOCK_ID)
# filter by: blocks_rll, blocks_dnr (chr list)


RunMoranResiduals <- function(model, model_df, blocks_sf) {
  
  # 1. Attach geometry
  sf_obj <- model_df %>%
    left_join(
      blocks_sf %>% dplyr::select(atlas_block, geometry),
      by = "atlas_block"
    ) %>%
    st_as_sf()
  
  # 2. Sort to ensure alignment
  sf_obj <- sf_obj %>%
    arrange(atlas_block)
  
  # 3. Extract Pearson residuals
  sf_obj$resid <- residuals(model, type = "pearson")
  
  # 4. Build neighbors (queen contiguity)
  nb <- poly2nb(sf_obj, queen = TRUE)
  
  # Check isolated blocks
  if (any(card(nb) == 0)) {
    message("Some blocks have zero neighbors.")
  }
  
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # 5. Moran's I
  moran.test(sf_obj$resid, lw, zero.policy = TRUE)
}




moran_results <- lapply(names(global_glm_models), function(nm) {
  
  RunMoranResiduals(
    model    = global_glm_models[[nm]],
    model_df = data_dir[[nm]],
    blocks_sf = blocks_all_sf
  )
})

names(moran_results) <- names(global_glm_models)







RunMoranResidual <- function(model, model_df, blocks_sf, k = 4) {
  
  sf_obj <- model_df %>%
    left_join(
      blocks_sf %>% select(atlas_block, geometry),
      by = "atlas_block"
    ) %>%
    st_as_sf() %>%
    arrange(atlas_block)
  
  sf_obj$resid <- residuals(model, type = "pearson")
  
  # Use centroids for kNN
  coords <- st_coordinates(st_centroid(sf_obj))
  
  nb <- knn2nb(knearneigh(coords, k = k))
  lw <- nb2listw(nb, style = "W")
  
  moran.test(sf_obj$resid, lw)
}
