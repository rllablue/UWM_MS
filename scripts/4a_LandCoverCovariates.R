#############
### SETUP ###
#############

## --- PACKAGES --- ##

library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(stringr)
library(ggplot2)
library(ggridges)
library(ggsignif)
library(sf)
library(raster)
library(terra)
library(exactextractr)
library(GGally)
library(corrplot)


## -- DATAFRAMES --- ###

# Create new df with only covariate values to be used in modeling 
wibba_covars_raw <- wibba_summary_comp %>%
  dplyr::select(atlas_block)


### --- FUNCTIONS --- ###

# Helper to clip buffered atlas outline from full NLCD 
Crop_mask_buffer <- function(raster_obj, polygon_sf, buffer_dist = 1000) {
  bbox_sf <- st_as_sfc(st_bbox(polygon_sf)) # Create bbox polygon from input sf polygon
  buffered_bbox_sf <- st_buffer(bbox_sf, dist = buffer_dist) # Buffer the bbox (units same as polygon crs)
  buffered_bbox_sf_proj <- st_transform(buffered_bbox_sf, crs(raster_obj)) # Reproject buffered bbox to raster CRS
  buffered_bbox_vect <- vect(buffered_bbox_sf_proj) # Convert to terra::vector
  raster_cropped <- terra::crop(raster_obj, buffered_bbox_vect) # Crop, mask raster
  raster_masked <- terra::mask(raster_cropped, buffered_bbox_vect)
  return(raster_masked)
}



#########################
### LAND COVER (NLCD) ###
#########################

# Results in raw % land type coverage per block in landcover and modeling covariate df


### --- SINGLE REFERENCE YEAR --- ###
# Land cover for intermediate year (2008) between Atlas periods

# --- FILES, SET UP --- #

# Load NLCD raster (crs custom AEA)
nlcd_2008_raw <- rast("data/maps/nlcd/NLCD_LndCov_2008.tif")

# Clip, mask NLCD raster to buffered Atlas outline
wi_nlcd_2008_buffered <- Crop_mask_buffer(nlcd_2008_raw, atlas_outline_crs3071, buffer_dist = 1000)

# Reproject cropped wi rasters back to block crs (EPSG:3071)
wi_nlcd_2008_crs3071 <- terra::project(wi_nlcd_2008_buffered, crs(blocks_comp_shp), method = "near")

# Convert terra raster to raster object (exactextractr compatability)
wi_nlcd_2008_rast <- raster::raster(wi_nlcd_2008_crs3071)
crs(wi_nlcd_2008_rast) <- crs(wi_nlcd_2008_crs3071) # force crs


# --- FUNCTION --- #

# Function to extract land cover 
Extract_landcover <- function(nlcd_raster, blocks_sf, year_suffix) {
  
  # Extract raster values per polygon
  nlcd_extract <- exactextractr::exact_extract(nlcd_raster, blocks_sf)
  valid_classes <- c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95)
  
  # Summarize class coverage by block
  landcover_summary <- purrr::map2_dfr(
    nlcd_extract,
    blocks_sf$atlas_block,
    ~ {
      tibble(class = .x$value, coverage = .x$coverage_fraction) %>%
        filter(class %in% valid_classes) %>%
        group_by(class) %>%
        summarize(class_cover = sum(coverage, na.rm = TRUE), .groups = "drop") %>%
        mutate(atlas_block = .y,
          percent_cover = 100 * class_cover / sum(class_cover))
    }
  )
  
  # Pivot wide, rename, group
  landcover_named <- landcover_summary %>%
    dplyr::select(atlas_block, class, percent_cover) %>%
    tidyr::pivot_wider(
      names_from = class,
      values_from = percent_cover,
      names_prefix = "nlcd_"
    ) %>%
    mutate(
      across(
        starts_with("nlcd_"), ~replace_na(.x, 0)
      )
    ) %>%
    rename(
      water_open        = nlcd_11,
      developed_open    = nlcd_21,
      developed_low     = nlcd_22,
      developed_med     = nlcd_23,
      developed_high    = nlcd_24,
      barren_land       = nlcd_31,
      forest_deciduous  = nlcd_41,
      forest_evergreen  = nlcd_42,
      forest_mixed      = nlcd_43,
      shrub_scrub       = nlcd_52,
      grassland         = nlcd_71,
      pasture           = nlcd_81,
      cropland          = nlcd_82,
      wetlands_woody    = nlcd_90,
      wetlands_herb     = nlcd_95
    ) %>%
    mutate(
      developed_lower   = developed_open + developed_low,
      developed_upper   = developed_med + developed_high,
      pasture_crop      = pasture + cropland,
      developed_total   = developed_open + developed_low + developed_med + developed_high,
      forest_total      = forest_deciduous + forest_evergreen + forest_mixed,
      wetlands_total    = wetlands_woody + wetlands_herb
    )
  
  # Combine raw + z, append year suffix
  final_df <- landcover_named %>%
    rename_with(~ paste0(., "_", year_suffix), -atlas_block)
  
  return(final_df)
}


### --- APPLY --- ###

landcover_2008 <- Extract_landcover(
  nlcd_raster = wi_nlcd_2008_rast,
  blocks_sf = blocks_comp_shp,
  year_suffix = "2008"
)


# OPTIONAL: # PARALLELIZE PROCESS #
num_cores <- detectCores() # examine cores
print(num_cores)

future::plan(multisession, workers = 10) # set # cores for processing

Parallel_land2008 <- future({ # Assign process to run in parallel
  Extract_landcover(
    nlcd_raster = wi_nlcd_2008_rast,
    blocks_sf = blocks_comp_shp,
    year_suffix = "2008"
  )
})

resolved(Parallel_land2008) # Check completion
landcover_summary_2008 <- value(Parallel_land2008) # store results





### --- COVER DIFFERENCE --- ###

# Land cover change between two years (2015, 1995) overlapping with Atlas periods

### SET-UP ###

# Load NLCD rasters (crs custom AEA)
nlcd_1995_raw <- rast("data/maps/nlcd/NLCD_LndCov_1995.tif")
nlcd_2015_raw <- rast("data/maps/nlcd/NLCD_LndCov_2015.tif")

# Clip, mask NLCD rasters to buffered Atlas outline
wi_nlcd_1995_buffered <- Crop_mask_buffer(nlcd_1995_raw, atlas_outline_crs3071, buffer_dist = 1000)
wi_nlcd_2015_buffered <- Crop_mask_buffer(nlcd_2015_raw, atlas_outline_crs3071, buffer_dist = 1000)

# Reproject cropped wi rasters back to block crs (EPSG:3071)
wi_nlcd_1995_crs3071 <- terra::project(wi_nlcd_1995_buffered, crs(blocks_comp_shp), method = "near")
wi_nlcd_2015_crs3071 <- terra::project(wi_nlcd_2015_buffered, crs(blocks_comp_shp), method = "near")

# Save clipped, buffered rasters
writeRaster(wi_nlcd_1995_crs3071, "outputs/maps/wi_nlcd_1995_buff_crs3071.tif", overwrite = TRUE)
writeRaster(wi_nlcd_2015_crs3071, "outputs/maps/wi_nlcd_2015_buff_crs3071.tif", overwrite = TRUE)

# Convert terra raster to raster object (exactextractr compatability)
wi_nlcd_1995_rast <- raster::raster(wi_nlcd_1995_crs3071)
crs(wi_nlcd_1995_rast) <- crs(wi_nlcd_1995_crs3071) # force crs

wi_nlcd_2015_rast <- raster::raster(wi_nlcd_2015_crs3071)
crs(wi_nlcd_2015_rast) <- crs(wi_nlcd_2015_crs3071)


### EXTRACTION ###

# Function to extract land cover per block
Extract_landcover_diff <- function(raster_early, raster_late, blocks_sf, years_suffix = c("1995", "2015")) {
  
  # Extract single raster
  Extract_single <- function(nlcd_raster, blocks_sf) {
    nlcd_extract <- exactextractr::exact_extract(nlcd_raster, blocks_sf)
    valid_classes <- c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95)
    
    landcover_summary <- purrr::map2_dfr(
      nlcd_extract, blocks_sf$atlas_block,
      ~ {
        tibble(class = .x$value, coverage = .x$coverage_fraction) %>%
          filter(class %in% valid_classes) %>%
          group_by(class) %>%
          summarize(class_cover = sum(coverage, na.rm = TRUE), .groups = "drop") %>%
          mutate(
            atlas_block = .y,
            percent_cover = 100 * class_cover / sum(class_cover)
          )
      }
    )
    
    landcover_named <- landcover_summary %>%
      dplyr::select(atlas_block, class, percent_cover) %>%
      pivot_wider(
        names_from = class,
        values_from = percent_cover,
        names_prefix = "nlcd_"
      ) %>%
      mutate(
        across(
          starts_with("nlcd_"), ~replace_na(.x, 0)
        )
      ) %>%
      rename(
        water_open        = nlcd_11,
        developed_open    = nlcd_21,
        developed_low     = nlcd_22,
        developed_med     = nlcd_23,
        developed_high    = nlcd_24,
        barren_land       = nlcd_31,
        forest_deciduous  = nlcd_41,
        forest_evergreen  = nlcd_42,
        forest_mixed      = nlcd_43,
        shrub_scrub       = nlcd_52,
        grassland         = nlcd_71,
        pasture           = nlcd_81,
        cropland          = nlcd_82,
        wetlands_woody    = nlcd_90,
        wetlands_herb     = nlcd_95
      ) %>%
      mutate(
        developed_lower   = developed_open + developed_low,
        developed_upper   = developed_med + developed_high,
        pasture_crop      = pasture + cropland,
        developed_total   = developed_open + developed_low + developed_med + developed_high,
        forest_total      = forest_deciduous + forest_evergreen + forest_mixed,
        wetlands_total    = wetlands_woody + wetlands_herb
      )
    
    landcover_named
  }
  
  ys_early <- years_suffix[1]
  ys_late  <- years_suffix[2]
  
  land_early <- Extract_single(raster_early, blocks_sf)
  land_late  <- Extract_single(raster_late,  blocks_sf)
  
  land_early <- rename_with(land_early, ~ paste0(., "_", ys_early), -atlas_block)
  land_late  <- rename_with(land_late,  ~ paste0(., "_", ys_late), -atlas_block)
  
  base_vars <- c(
    "water_open", "developed_open", "developed_low", "developed_med", "developed_high",
    "barren_land", "forest_deciduous", "forest_evergreen", "forest_mixed",
    "shrub_scrub", "grassland", "pasture", "cropland",
    "wetlands_woody", "wetlands_herb",
    "developed_lower","developed_upper","pasture_crop",
    "developed_total","forest_total","wetlands_total"
  )
  
  # Calculate differences
  land_diff <- land_early %>%
    inner_join(land_late, by = "atlas_block")
  
  for(var in base_vars){
    land_diff[[paste0(var, "_diff")]] <- land_diff[[paste0(var, "_", ys_late)]] - land_diff[[paste0(var, "_", ys_early)]]
  }
  
  # Return all three dfs
  list(
    landcover_early = land_early,
    landcover_late  = land_late,
    landcover_diff  = land_diff
  )
}



### --- APPLY --- ################################# Note: need to fix above function to not include 2015, 1995 values in _diff df
# List with all data
landcover_9515 <- Extract_landcover_diff(
    raster_early = wi_nlcd_1995_rast,
    raster_late  = wi_nlcd_2015_rast,
    blocks_sf    = blocks_comp_shp
)
  



# OPTIONAL: PARALLELIZE PROCESS #

future::plan(multisession, workers = 10) 

# Assign process to run in parallel
Parallel_land9515 <- future({
  Extract_landcover_diff(
    raster_early = wi_nlcd_1995_rast,
    raster_late  = wi_nlcd_2015_rast,
    blocks_sf    = blocks_comp_shp
  )
})

resolved(Parallel_land9515)
landcover_summary_9515 <- value(Parallel_land9515)






### --- SUMMARIZE --- ###

# List to dfs
landcover_1995 <- landcover_9515[[1]]
landcover_2015 <- landcover_9515[[2]]
landcover_diff <- landcover_9515[[3]]
  
landcover_diff <- landcover_diff %>%
  dplyr::select(-ends_with("_1995"), -ends_with("_2015"))


# Join to modeling covariate df
wibba_covars_raw <- wibba_covars_raw %>%
  left_join(
    landcover_2008 %>%
      rename_with(~ str_replace(.x, "_2008$", "_base")),
    by = "atlas_block"
  ) %>%
  left_join(
    landcover_diff,
    by = "atlas_block"
  )





###################
### STATS TESTS ###
###################

### --- VIF --- ###

# Multi-collinearity testing

covars_all <- wibba_modeling_comp %>%
  dplyr::select(-atlas_block)

cor_mat <- cor(covars_all, use = "pairwise.complete.obs")

corrplot::corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  tl.cex = 0.7,
  tl.col = "black",
  order = "hclust",
  addCoef.col = "black",
  number.cex = 0.5
)



# Simple linear model including all covariates
vif_mod <- lm(
  # dummy dependent variable (not used, just for VIF computation)
  as.numeric(rep(1, nrow(wibba_modeling_comp))) ~ .,
  data = covars_all
)

# Compute VIFs
vif_values <- car::vif(vif_mod)
vif_values <- sort(vif_values, decreasing = TRUE)
print(vif_values)





#################
### VISUALIZE ###
#################

# -----------------------------
# 1. Define landcover types to exclude and grouped classes
# -----------------------------
exclude_types <- c("developed_open", "developed_low", "developed_med",
                   "developed_high", "pasture", "cropland")
grouped_classes <- c("developed_lower", "developed_upper", "pasture_crop")

# -----------------------------
# 2. Combine and tidy individual year data
# -----------------------------
landcover_z_long <- bind_rows(
  wi_landcover1995 %>% dplyr::mutate(year = 1995),
  wi_landcover2008 %>% dplyr::mutate(year = 2008),
  wi_landcover2015 %>% dplyr::mutate(year = 2015)
) %>%
  dplyr::select(atlas_block, year, dplyr::ends_with("_z")) %>%
  tidyr::pivot_longer(
    cols = dplyr::ends_with("_z"),
    names_to = "landcover_type",
    values_to = "z_score"
  ) %>%
  dplyr::mutate(
    landcover_type = stringr::str_remove(landcover_type, "_z"),
    landcover_type_clean = stringr::str_replace_all(landcover_type, "_", " "),
    category = dplyr::if_else(landcover_type %in% grouped_classes, "Grouped", "Individual")
  ) %>%
  dplyr::filter(!landcover_type %in% exclude_types)

# Order landcover types: grouped first, then alphabetically
landcover_order <- landcover_z_long %>%
  dplyr::distinct(landcover_type, category) %>%
  dplyr::arrange(category, landcover_type) %>%
  dplyr::pull(landcover_type)

landcover_z_long <- landcover_z_long %>%
  dplyr::mutate(landcover_type = factor(landcover_type, levels = landcover_order))

# -----------------------------
# 3. Ridge plot for individual years
# -----------------------------
ggplot2::ggplot(landcover_z_long, aes(x = z_score, y = as.factor(year), fill = as.factor(year))) +
  ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2) +
  ggplot2::facet_wrap(~ landcover_type_clean, scales = "free", ncol = 4) +
  ggplot2::scale_fill_brewer(palette = "Dark2") +
  ggplot2::labs(
    title = "Distribution of Land Cover Z-Scores Across Years",
    x = "Z-score",
    y = "Year"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold"),
    panel.spacing = grid::unit(0.6, "lines")
  )

# -----------------------------
# 4. Prepare 1995 → 2015 difference data
# -----------------------------
landcover_diff_9515 <- landcover_z_long %>%
  dplyr::filter(year %in% c(1995, 2015)) %>%
  dplyr::select(atlas_block, year, landcover_type, z_score) %>%
  tidyr::pivot_wider(names_from = year, values_from = z_score) %>%
  dplyr::mutate(diff_9515 = `2015` - `1995`) %>%
  dplyr::select(atlas_block, landcover_type, diff_9515)

# Wilcoxon signed-rank test for differences
wilcox_diff <- landcover_diff_9515 %>%
  dplyr::group_by(landcover_type) %>%
  dplyr::summarise(
    p_value = stats::wilcox.test(diff_9515, mu = 0, paired = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# -----------------------------
# 5. Ridge plot of 1995 → 2015 differences
# -----------------------------
max_diff <- max(landcover_diff_9515$diff_9515, na.rm = TRUE)

ggplot2::ggplot(landcover_diff_9515, aes(x = diff_9515, y = landcover_type, fill = landcover_type)) +
  ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_text(
    data = wilcox_diff,
    aes(y = landcover_type, label = sig_label),
    x = max_diff,
    inherit.aes = FALSE,
    hjust = 1.1,
    size = 4
  ) +
  ggplot2::labs(
    title = "Change in Landcover Z-Scores (1995 & 2015) Across Blocks",
    x = "Z-score difference (2015 - 1995)",
    y = "Landcover Type"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

# -----------------------------
# 6. Friedman and pairwise Wilcoxon tests
# -----------------------------
landcover_types <- landcover_z_long %>%
  dplyr::distinct(landcover_type) %>%
  dplyr::pull()

# Friedman test
friedman_results <- purrr::map_dfr(landcover_types, function(lc) {
  mat <- landcover_z_long %>%
    dplyr::filter(landcover_type == lc) %>%
    dplyr::select(atlas_block, year, z_score) %>%
    tidyr::pivot_wider(names_from = year, values_from = z_score) %>%
    dplyr::arrange(atlas_block) %>%
    dplyr::select(-atlas_block) %>%
    as.matrix()
  res <- stats::friedman.test(mat)
  tibble::tibble(
    landcover_type = lc,
    chi_sq = unname(res$statistic),
    df = unname(res$parameter),
    p_value = res$p.value
  )
})

friedman_results


# Pairwise Wilcoxon
pairwise_wilcoxon_results <- purrr::map_dfr(landcover_types, function(lc) {
  lc_data <- landcover_z_long %>%
    dplyr::filter(landcover_type == lc) %>%
    dplyr::select(atlas_block, year, z_score) %>%
    tidyr::pivot_wider(names_from = year, values_from = z_score) %>%
    dplyr::arrange(atlas_block)
  mat <- as.matrix(dplyr::select(lc_data, -atlas_block))
  yrs <- colnames(mat)
  combs <- combn(yrs, 2, simplify = FALSE)
  pw_tests <- purrr::map_dfr(combs, function(pair) {
    x <- mat[, pair[1]]
    y <- mat[, pair[2]]
    test <- stats::wilcox.test(x, y, paired = TRUE, exact = FALSE)
    tibble::tibble(
      landcover_type = lc,
      year1 = pair[1],
      year2 = pair[2],
      W_stat = unname(test$statistic),
      p_value_raw = test$p.value
    )
  })
  pw_tests %>% dplyr::mutate(p_value_adj = p.adjust(p_value_raw, method = "bonferroni"))
})

pairwise_wilcox_results





### REAL CHANGE v SIGNIFICANCE ###

# -----------------------------
# 1. Prepare difference data (2015 - 1995)
# -----------------------------
exclude_types <- c("developed_open", "developed_low", "developed_med",
                   "developed_high", "pasture", "cropland")
grouped_classes <- c("developed_lower", "developed_upper", "pasture_crop")

land_diff_long <- wi_land9515_diff %>%
  dplyr::select(atlas_block, dplyr::ends_with("_diff")) %>%
  tidyr::pivot_longer(
    cols = dplyr::ends_with("_diff"),
    names_to = "landcover_type",
    values_to = "diff_percent"
  ) %>%
  dplyr::mutate(
    landcover_type = stringr::str_remove(landcover_type, "_diff"),
    landcover_type_clean = stringr::str_replace_all(landcover_type, "_", " "),
    category = dplyr::if_else(landcover_type %in% grouped_classes, "Grouped", "Individual"),
    ecologically_significant = abs(diff_percent) >= 1,
    change_direction = dplyr::case_when(
      diff_percent > 0 ~ "Increase",
      diff_percent < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  ) %>%
  dplyr::filter(!landcover_type %in% exclude_types)

# Order landcover types
landcover_order <- land_diff_long %>%
  dplyr::distinct(landcover_type, category) %>%
  dplyr::arrange(category, landcover_type) %>%
  dplyr::pull(landcover_type)

land_diff_long <- land_diff_long %>%
  dplyr::mutate(landcover_type = factor(landcover_type, levels = landcover_order))

# -----------------------------
# 2. Wilcoxon significance
# -----------------------------
wilcox_9515_diff <- land_diff_long %>%
  dplyr::group_by(landcover_type) %>%
  dplyr::summarise(
    p_value = stats::wilcox.test(diff_percent, mu = 0, paired = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# -----------------------------
# 3. Median, min, max, and % blocks ≥1%
# -----------------------------
summary_annotations <- land_diff_long %>%
  dplyr::group_by(landcover_type, landcover_type_clean) %>%
  dplyr::summarise(
    median_diff = median(diff_percent, na.rm = TRUE),
    min_diff    = min(diff_percent, na.rm = TRUE),
    max_diff    = max(diff_percent, na.rm = TRUE),
    pct_ecol_sig = mean(ecologically_significant, na.rm = TRUE) * 100,
    .groups = "drop"
  )
# -----------------------------
# 4. Ridge plot colored by direction of change
# -----------------------------
ggplot(land_diff_long, aes(x = diff_percent, y = landcover_type_clean, fill = change_direction)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(
    data = summary_annotations,
    aes(x = min_diff, y = landcover_type_clean),
    color = "#e41a1c", shape = 21, fill = "#e41a1c", size = 2
  ) +
  geom_point(
    data = summary_annotations,
    aes(x = max_diff, y = landcover_type_clean),
    color = "#4daf4a", shape = 21, fill = "#4daf4a", size = 2
  ) +
  # Median + min/max + % blocks ≥1% annotation
  geom_text(
    data = summary_annotations,
    aes(
      x = max(land_diff_long$diff_percent, na.rm = TRUE),
      y = landcover_type_clean,
      label = paste0(
        "Median: ", round(median_diff, 2), "%\n",
        "Min: ", round(min_diff, 2), "%  Max: ", round(max_diff, 2), "%\n",
        round(pct_ecol_sig, 1), "% ≥1%"
      )
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    size = 3.5
  ) +
  # Wilcoxon significance on the left
  geom_text(
    data = wilcox_9515_diff %>%
      dplyr::mutate(landcover_type_clean = stringr::str_replace_all(landcover_type, "_", " ")),
    aes(
      x = min(land_diff_long$diff_percent, na.rm = TRUE),
      y = landcover_type_clean,
      label = sig_label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Increase" = "#4daf4a", "Decrease" = "#e41a1c", "No change" = "grey70"),
    name = "Direction of Change"
  ) +
  labs(
    title = "Distribution of Landcover % Change Across Blocks (Between 1995 and 2015)",
    x = "% Cover Change (2015-1995)",
    y = "Landcover Type",
    subtitle = "Black dashed @ 0 = no change; Blue dashed = ±1% ecological threshold\nAsterisks = Wilcoxon significance"
  ) +
  theme_ridges() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )