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

