#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

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


## -- MODELING COVARIATE DF --- ###

# Create new summary df with only covariate values to be used in modeling 

wibba_modeling_comp <- wibba_summary_comp %>%
  dplyr::select(atlas_block)



#########################
### LAND COVER (NLCD) ###
#########################

### --- CLIP RASTERS --- ###

# Create function to clip buffered atlas outline from full NLCD 
crop_mask_buffered_bbox <- function(raster_obj, polygon_sf, buffer_dist = 1000) {
  bbox_sf <- st_as_sfc(st_bbox(polygon_sf)) # Create bbox polygon from input sf polygon
  buffered_bbox_sf <- st_buffer(bbox_sf, dist = buffer_dist) # Buffer the bbox (units same as polygon crs)
  buffered_bbox_sf_proj <- st_transform(buffered_bbox_sf, crs(raster_obj)) # Reproject buffered bbox to raster CRS
  buffered_bbox_vect <- vect(buffered_bbox_sf_proj) # Convert to terra::vector
  raster_cropped <- terra::crop(raster_obj, buffered_bbox_vect) # Crop, mask raster
  raster_masked <- terra::mask(raster_cropped, buffered_bbox_vect)
  return(raster_masked)
}


### --- COVER DIFFERENCE --- ###

# Land cover change between two years (2015, 1995) overlapping with Atlase periods


### SET-UP ###
# Load NLCD rasters (crs custom AEA)
nlcd_1995_raw <- rast("data/maps/nlcd/NLCD_LndCov_1995.tif")
nlcd_2015_raw <- rast("data/maps/nlcd/NLCD_LndCov_2015.tif")

# Clip, mask NLCD rasters to buffered Atlas outline
wi_nlcd_1995_buffered <- crop_mask_buffered_bbox(nlcd_1995_raw, atlas_outline_crs3071, buffer_dist = 1000)
wi_nlcd_2015_buffered <- crop_mask_buffered_bbox(nlcd_2015_raw, atlas_outline_crs3071, buffer_dist = 1000)

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
Extract_landcover_diff <- function(raster_early, raster_late, blocks_sf) {
  
  # Extract single raster
  extract_single <- function(nlcd_raster, blocks_sf) {
    
    nlcd_extract <- exactextractr::exact_extract(nlcd_raster, blocks_sf)
    valid_classes <- c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95)
    
    landcover_summary <- purrr::map2_dfr(
      nlcd_extract,
      blocks_sf$atlas_block,
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
      mutate(across(starts_with("nlcd_"), ~replace_na(.x, 0))) %>%
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
      )
    
    # Grouped categories
    landcover_grouped <- landcover_named %>%
      mutate(
        developed_lower     = developed_open + developed_low,
        developed_upper     = developed_med + developed_high,
        pasture_crop        = pasture + cropland,
      )
    
    # Z-scores for both individual and grouped columns
    landcover_z <- landcover_grouped %>%
      mutate(
        across(
          .cols = -atlas_block,             # all numeric columns except the ID
          .fns = ~ {
            z <- as.numeric(scale(.))
            z[is.na(z)] <- 0               # replace NAs (from SD=0) with 0
            z
          },
          .names = "{.col}_z"
        )
      )
    
    return(landcover_z)
  }
  
  # Extract early and late years
  land_early <- extract_single(raster_early, blocks_sf)
  land_late  <- extract_single(raster_late,  blocks_sf)
  
  # Compute differences
  land_diff <- land_early %>%
    inner_join(land_late, by = "atlas_block", suffix = c("_early", "_late")) %>%
    # Individual diffs
    mutate(
      water_open_diff       = water_open_late       - water_open_early,
      developed_open_diff   = developed_open_late   - developed_open_early,
      developed_low_diff    = developed_low_late    - developed_low_early,
      developed_med_diff    = developed_med_late    - developed_med_early,
      developed_high_diff   = developed_high_late   - developed_high_early,
      barren_land_diff      = barren_land_late      - barren_land_early,
      forest_deciduous_diff = forest_deciduous_late - forest_deciduous_early,
      forest_evergreen_diff = forest_evergreen_late - forest_evergreen_early,
      forest_mixed_diff     = forest_mixed_late     - forest_mixed_early,
      shrub_scrub_diff      = shrub_scrub_late      - shrub_scrub_early,
      grassland_diff        = grassland_late        - grassland_early,
      pasture_diff          = pasture_late          - pasture_early,
      cropland_diff         = cropland_late         - cropland_early,
      wetlands_woody_diff   = wetlands_woody_late   - wetlands_woody_early,
      wetlands_herb_diff    = wetlands_herb_late    - wetlands_herb_early
    ) %>%
    # Grouped diffs
     mutate(
      developed_lower_diff    = developed_open_diff + developed_low_diff,
      developed_upper_diff    = developed_med_diff + developed_high_diff,
      pasture_crop_diff       = pasture_diff + cropland_diff,
    ) %>%
    # Z-scores for all diffs
    mutate(
      across(
        .cols = ends_with("_diff"),
        .fns = ~ as.numeric(scale(.)),
        .names = "{.col}_z"
      )
    )
  
  # Return all three dataframes
  return(list(
    landcover_early = land_early,
    landcover_late  = land_late,
    landcover_diff  = land_diff %>%
      dplyr::select(atlas_block, ends_with("_diff"), ends_with("_diff_z"))
  ))
}

# Apply function
land9515_all <- Extract_landcover_diff(
  raster_early = wi_nlcd_1995_rast,
  raster_late  = wi_nlcd_2015_rast,
  blocks_sf    = blocks_comp_shp
)

# Create separate dfs (all w/ individual, grouped raw, z-standardized cover % per block)
wi_landcover1995 <- land9515_all$landcover_early    
wi_landcover2015 <- land9515_all$landcover_late      
wi_land9515_diff <- land9515_all$landcover_diff     





### --- SINGLE REFERENCE YEAR --- ###
# Land cover for intermediate year (2008) between Atlas periods

# Load NLCD raster (crs custom AEA)
nlcd_2008_raw <- rast("data/maps/nlcd/NLCD_LndCov_2008.tif")

# Clip, mask NLCD raster to buffered Atlas outline
wi_nlcd_2008_buffered <- crop_mask_buffered_bbox(nlcd_2008_raw, atlas_outline_crs3071, buffer_dist = 1000)

# Reproject cropped wi rasters back to block crs (EPSG:3071)
wi_nlcd_2008_crs3071 <- terra::project(wi_nlcd_2008_buffered, crs(blocks_comp_shp), method = "near")

# Convert terra raster to raster object (exactextractr compatability)
wi_nlcd_2008_rast <- raster::raster(wi_nlcd_2008_crs3071)
crs(wi_nlcd_2008_rast) <- crs(wi_nlcd_2008_crs3071) # force crs

# Function to extract land cover 
Extract_landcover <- function(nlcd_raster, blocks_sf) {
  
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
  
  # Pivot to wide format
  landcover_wide <- landcover_summary %>%
    dplyr::select(atlas_block, class, percent_cover) %>%
    tidyr::pivot_wider(
      names_from = class,
      values_from = percent_cover,
      names_prefix = "nlcd_"
    ) %>%
    mutate(across(starts_with("nlcd_"), ~replace_na(.x, 0)))
  
  # Rename NLCD classes
  landcover_named <- landcover_wide %>%
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
    )
  
  # Add grouped totals
  landcover_grouped <- landcover_named %>%
    mutate(
      developed_lower     = developed_open + developed_low,
      developed_upper     = developed_med + developed_high,
      pasture_crop        = pasture + cropland,
    )
  
  # Z-score both individual and grouped columns
  landcover_z <- landcover_grouped %>%
    mutate(
      across(
        .cols = -atlas_block,             # all numeric columns except the ID
        .fns = ~ {
          z <- as.numeric(scale(.))
          z[is.na(z)] <- 0               # replace NAs (from SD=0) with 0
          z
        },
        .names = "{.col}_z"
      )
    )
  
  # Create results df (individual, grouped and raw, z-standardized cover % per block)
  final_df <- landcover_grouped %>%
    left_join(
      landcover_z %>% dplyr::select(atlas_block, ends_with("_z")),
      by = "atlas_block"
    )
  
  return(final_df)
}

# Extract land cover for reference year (2008), create summary
wi_landcover2008 <- Extract_landcover(wi_nlcd_2008_rast, blocks_comp_shp)



########################
### GROUP, SUMMARIZE ###
########################

# Select relevant covariates for modeling
landcover_cols_2008z <- c(
  "atlas_block",
  "water_open_z",
  "forest_deciduous_z",
  "forest_evergreen_z",
  "forest_mixed_z",
  "shrub_scrub_z",
  "grassland_z",
  "wetlands_woody_z",
  "wetlands_herb_z",
  "developed_lower_z",
  "developed_upper_z",
  "pasture_crop_z"
)

landcover_cols_diffz <- c(
  "atlas_block",
  "water_open_diff_z",
  "forest_deciduous_diff_z",
  "forest_evergreen_diff_z",
  "forest_mixed_diff_z",
  "shrub_scrub_diff_z",
  "grassland_diff_z",
  "wetlands_woody_diff_z",
  "wetlands_herb_diff_z",
  "developed_lower_diff_z",
  "developed_upper_diff_z",
  "pasture_crop_diff_z"
)

# New dfs for filtered, renamed covariates
wiland_2008z_filt <- wi_landcover2008 %>%
  dplyr::select(dplyr::all_of(landcover_cols_2008z)) %>%
  dplyr::rename_with(
    ~ paste0(., "_08"),
    -atlas_block
  )

wiland_9515z_filt <- wi_land9515_diff %>%
  dplyr::select(dplyr::all_of(landcover_cols_diffz)) %>%
  dplyr::rename_with(
    ~ paste0(., "_9515"),
    -atlas_block
  )

# Merge with modeling df
wibba_modeling_comp <- wibba_modeling_comp %>%
  dplyr::left_join(wiland_2008z_filt, by = "atlas_block") %>%
  dplyr::left_join(wiland_9515z_filt, by = "atlas_block")











###################
### STATS TESTS ###
###################

### --- VIF --- ###

# Multicollinearity testing

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




# --- Gather data for easier plotting ---
cor_long <- wibba_modeling_comp %>%
  # Select only z-scored landcover columns with year suffixes
  dplyr::select(matches("_(08|9515)$")) %>%
  # Reshape long for comparison
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("type", "year"),
    names_pattern = "(.*)_(08|9515)$"
  ) %>%
  tidyr::pivot_wider(names_from = year, values_from = value) %>%
  # Convert to numeric just in case (fix for 'x must be numeric')
  dplyr::mutate(
    `08` = as.numeric(`08`),
    `9515` = as.numeric(`9515`)
  ) %>%
  # Drop NA pairs and variables with no variation
  dplyr::filter(!is.na(`08`), !is.na(`9515`)) %>%
  dplyr::group_by(type) %>%
  dplyr::filter(
    dplyr::n_distinct(`08`) > 1,
    dplyr::n_distinct(`9515`) > 1
  ) %>%
  dplyr::ungroup()

# --- Compute correlation coefficients for annotation ---
cor_labels <- cor_long %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(
    r = cor(`08`, `9515`, use = "pairwise.complete.obs"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = paste0("r = ", round(r, 2)))

# --- Plot ---
ggplot(cor_long, aes(x = `08`, y = `9515`)) +
  geom_point(alpha = 0.4, size = 1.2, color = "darkgray") +
  geom_smooth(method = "lm", se = FALSE, color = "darkblue", linewidth = 0.6) +
  facet_wrap(~type, scales = "free") +
  geom_text(
    data = cor_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.2,
    size = 3.5, color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Z-score (2008 land cover)",
    y = "Z-score (1995–2015 change)",
    title = "Relationships between static (2008) and change (1995–2015) land cover covariates"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )



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