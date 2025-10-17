#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(sf)
library(raster)
library(terra)
library(exactextractr)


#########################
### SET UP COVARIATES ###
#########################

## --- WIBBA BLOCKS --- ##

# Atlas blocks shp (crs EPSG:3071)
blocks_shp <- st_read("data/maps/wibba blocks/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # load full block map
  rename(atlas_block = BLOCK_ID)

blocks_comp_shp <- blocks_shp %>% # create sf object of comparable blocks
  filter(atlas_block %in% blocks_comp)

# Create Atlas outline polygon
atlas_outline_crs3071 <- blocks_shp %>% 
  st_union() %>%      # dissolve all blocks into one geometry
  st_as_sf()          # convert back to sf object
st_write(atlas_outline_crs3071, "outputs/maps/atlas_outline_crs3071.shp", append = FALSE)

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


## --- LAND COVER (NLCD) --- ##

## LAND COVER REFERENCE ##
# land cover in single intermediate year (2008) between Atlas periods

# Load NLCD raster (crs custom AEA)
nlcd_2008_raw <- rast("data/maps/nlcd/NLCD_LndCov_2008.tif")

# Clip, mask NLCD raster to buffered Atlas outline
wi_nlcd_2008_buffered <- crop_mask_buffered_bbox(nlcd_2008_raw, atlas_outline_crs3071, buffer_dist = 1000)

# Reproject cropped wi rasters back to block crs (EPSG:3071)
wi_nlcd_2008_crs3071 <- terra::project(wi_nlcd_2008_buffered, crs(blocks_comp_shp), method = "near")

# Convert terra raster to raster object (exactextractr compatability)
wi_nlcd_2008_rast <- raster::raster(wi_nlcd_2008_crs3071)
crs(wi_nlcd_2008_rast) <- crs(wi_nlcd_2008_crs3071) # force crs

# Function to extract landcover
extract_landcover <- function(nlcd_raster, blocks_sf) {
  nlcd_extract <- exactextractr::exact_extract(nlcd_raster, blocks_sf) # Extract raster object values, coverage fractions per sf polygon
  valid_classes <- c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95) # Define valid NLCD classes (N/A = class 12: perennial ice/snow)
  landcover_summary <- purrr::map2_dfr( # Summarize, repeat per block
    nlcd_extract,
    blocks_sf$atlas_block,
    ~ {
      df <- tibble(class = .x$value, coverage = .x$coverage_fraction) %>% # Filter valid classes
        filter(class %in% valid_classes) %>%
        group_by(class) %>%
        summarize(class_cover = sum(coverage, na.rm = TRUE), .groups = "drop") %>%
        mutate(atlas_block = .y) %>%
        mutate(percent_cover = 100 * class_cover / sum(class_cover))
      return(df)
    }
  )
  
  landcover_wide <- landcover_summary %>% # Pivot wide
    dplyr::select(atlas_block, class, percent_cover) %>%
    tidyr::pivot_wider(
      names_from = class,
      values_from = percent_cover,
      names_prefix = "nlcd_"
    ) %>%
    mutate(across(starts_with("nlcd_"), ~replace_na(.x, 0)))
  
  landcover_named <- landcover_wide %>% # NLCD codes to names
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
  
  return(landcover_named)
}

# Extract land cover for reference year (2008)
wi_landcover_2008 <- extract_landcover(wi_nlcd_2008_rast, blocks_comp_shp)

wi_landcover_2008 <- wi_landcover_2008 %>%
  rename_with(.fn = ~ paste0(.x, "_2008"), .cols = 2:16) %>%
  mutate(
    across(
      .cols = 2:16,
      .fns = ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  )

# Create land cover summary table
land_summary_comp <- wi_landcover_2008 

# Join standardized value to WIBBA summary
wibba_summary_comp <- wibba_summary_comp %>%
  left_join(
    wi_landcover_2008 %>% 
      dplyr::select(atlas_block, ends_with("_z")),
    by = "atlas_block"
  )


### --- LAND COVER DIFFERENCE --- ##

### NLCD ###
# change in cover between reference years for each Atlas period

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

# Extract land cover for reference years (1995, 2015)
wi_landcover_1995 <- extract_landcover(wi_nlcd_1995_rast, blocks_comp_shp)
wi_landcover_2015 <- extract_landcover(wi_nlcd_2015_rast, blocks_comp_shp)

# Add to land cover summary table
landcover_summary_comp <- wi_landcover_1995 %>%
  inner_join(wi_landcover_2015, by = "atlas_block", suffix = c("_1995", "_2015")) %>%
  mutate(
    water_open_diff       = water_open_2015       - water_open_1995,
    developed_open_diff   = developed_open_2015   - developed_open_1995,
    developed_low_diff    = developed_low_2015    - developed_low_1995,
    developed_med_diff    = developed_med_2015    - developed_med_1995,
    developed_high_diff   = developed_high_2015   - developed_high_1995,
    barren_land_diff      = barren_land_2015      - barren_land_1995,
    forest_deciduous_diff = forest_deciduous_2015 - forest_deciduous_1995,
    forest_evergreen_diff = forest_evergreen_2015 - forest_evergreen_1995,
    forest_mixed_diff     = forest_mixed_2015     - forest_mixed_1995,
    shrub_scrub_diff      = shrub_scrub_2015      - shrub_scrub_1995,
    grassland_diff        = grassland_2015        - grassland_1995,
    pasture_diff          = pasture_2015          - pasture_1995,
    cropland_diff         = cropland_2015         - cropland_1995,
    wetlands_woody_diff   = wetlands_woody_2015   - wetlands_woody_1995,
    wetlands_herb_diff    = wetlands_herb_2015    - wetlands_herb_1995
  ) %>%
  mutate(across(ends_with("_diff"), ~ scale(.)[,1], .names = "{.col}_z")) # z-standardize

# Join to WIBBA summary
wibba_summary_comp <- wibba_summary_comp %>%
  left_join(
    land_summary_comp %>% 
      dplyr::select(atlas_block, ends_with("_diff_z")),
    by = "atlas_block"
  )


### WISCLAND ###







## --- PROTECTED AREA (PAD) --- ##

# UPLOAD FILES # 
blocks_comp_shp_5070 <- blocks_comp_shp %>% # reproject atlas blocks
  st_transform(5070)  # NAD83, Conus Albers EPSG:5070 (ideal for area calcs)

pad_wi_5070 <- st_read("data/maps/pad wi/pad_wi_polygons.shp") %>% # align pad crs
  st_transform(5070)

# FULL PAD #
# Calculate PA/block w/ dissolving
pad_wi_diss <- pad_wi_5070 %>% 
  st_union() %>% 
  st_sf(geometry = .) # attempt to remedy inaccurate pa %s over 100

blocks_pad_diss <- st_intersection(blocks_comp_shp_5070, pad_wi_diss) %>%
  mutate(overlap_area = st_area(geometry))

# Summarize PA by blocks
blocks_pad_summary_diss <- blocks_pad_diss %>%
  group_by(atlas_block) %>%
  summarise(pa_area = sum(overlap_area), .groups = "drop") %>% 
  left_join(
    blocks_comp_shp_5070 %>%
      mutate(block_area = st_area(.)) %>%
      st_drop_geometry() %>%
      dplyr::select(atlas_block, block_area),
    by = "atlas_block"
  ) %>%
  mutate(
    pa_area = if_else(is.na(pa_area), set_units(0, m^2), pa_area),
    block_area_km2 = set_units(block_area, km^2),
    pa_area_km2 = set_units(pa_area, km^2),
    pa_percent = as.numeric(pa_area_km2 / block_area_km2 * 100),
    pa_z = as.numeric(scale(pa_percent))
  ) %>%
  dplyr::select(atlas_block, pa_area_km2, block_area_km2, pa_percent, pa_z)

# Join to wibba summary df
wibba_summary_comp <- wibba_summary_comp %>% 
  left_join(
    blocks_pad_summary_diss %>% dplyr::select(atlas_block, pa_z),
    by = "atlas_block"
  )


####################
# calculate total area of each block
block_areas <- blocks_wibba %>%
  mutate(block_area = st_area(geometry)) %>%
  dplyr::select(BLOCK_ID, block_area)

# join and summarize 
blocks_pad_summary <- left_join(
  block_areas, 
  st_drop_geometry(blocks_pad_summary), 
  by = "BLOCK_ID"
) %>%
  mutate(
    pa_area = if_else(is.na(pa_area), set_units(0, m^2), pa_area),
    block_area_km2 = set_units(block_area, km^2),
    pa_area_km2 = set_units(pa_area, km^2),
    pa_percent = as.numeric(pa_area_km2 / block_area_km2 * 100)
  ) %>%
  dplyr::select(BLOCK_ID, pa_area_km2, block_area_km2, pa_percent)

# QAQC for pa% errors
problem_blocks <- blocks_pad_summary %>%
  filter(pa_percent > 100) %>%
  arrange(desc(pa_percent))

print(problem_blocks)

# pa sumamry df
pa_summary <- blocks_pad_summary %>% # extract pa % from overall area summary
  st_drop_geometry() %>%
  dplyr::select(BLOCK_ID, pa_percent) %>%
  rename(atlas_block = BLOCK_ID) %>%
  mutate(
    pa_z = as.numeric(scale(pa_percent))
  )

# join to other summaries
wibba_summary <- wibba_summary %>%
  left_join(pa_summary %>% dplyr::select(atlas_block, pa_z), by = "atlas_block")

rcki_wibba_summary <- rcki_wibba_summary %>%
  left_join(blocks_pa_percent, by = "atlas_block")

rbwo_wibba_summary <- rbwo_wibba_summary %>%
  left_join(blocks_pa_percent, by = "atlas_block")

# FILTERED PAD #
# Calculate PA/block w/ no dissolving (allows for land owner filtering)

# pad filtering [WIP]
pad_wi_filtered <- pad_wi_5070 %>%
  filter(
    Own_Type %in% c("FED", "STAT", "NGO", "LOC"),
    Own_Name %in% c("NGO", "CNTY", "CITY", "SDNR", "NPS", "USFS", "FWS", "NRCS", "BLM", "USBR"),
    Loc_Own %in% c("The Nature Conservancy", "National Audubon Society")
  )

blocks_pad_ndiss <- st_intersection(blocks_comp_shp_5070, pad_wi_5070) %>%
  mutate(overlap_area = st_area(.))

## WIP WIP WIP ###
# summarize non-dissolved by block
blocks_pad_summary_nd <- blocks_pad_ndiss %>%
  group_by(atlas_block) %>%
  summarise(pa_area = sum(overlap_area), .groups = "drop")

# calculate total block area
block_areas <- blocks_comp_shp_5070 %>%
  mutate(block_area = st_area(geometry)) %>%
  dplyr::select(atlas_block, block_area)

# join and calculate PA percent
blocks_pad_summary_nd <- left_join(
  block_areas, 
  st_drop_geometry(blocks_pad_summary_nd), 
  by = "atlas_block"
) %>%
  mutate(
    pa_area = if_else(is.na(pa_area), set_units(0, m^2), pa_area),
    block_area_km2 = set_units(block_area, km^2),
    pa_area_km2 = set_units(pa_area, km^2),
    pa_percent = as.numeric(pa_area_km2 / block_area_km2 * 100),
    pa_z = as.numeric(scale(pa_percent))
  ) %>%
  dplyr::select(atlas_block, pa_area_km2, block_area_km2, pa_percent, pa_z)



# --- ELEVATION DATA (NED) --- #
# mean elevation per block done in ArcPro (1/3 arc-sec DEMs,~10m)

elevation_summary <- read_csv("maps/ned/ned_wibba_table.csv") # read in data calc'd in ArcPro

elevation_summary <- elevation_summary %>%
  rename(
    atlas_block = BLOCK_ID, #rename variables to match
    elev_mean = MEAN
  ) %>%
  mutate(atlas_block = as.character(atlas_block)) %>%
  mutate(
    elev_mean_z = as.numeric(scale(elev_mean)) # z-standard
  )

wibba_summary <- wibba_summary %>%
  left_join(elevation_summary %>% dplyr::select(atlas_block, elev_mean_z), by = "atlas_block") # join to summary df



### --- CLIMATE DATA --- ###

# crs integration of daymet and atlas maps
daymet_crs <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" # Daymet custom crs (Lambert Conformal Conic)
blocks_wibba_daymet <- st_transform(blocks_wibba, crs = daymet_crs)

# --- TEMPERATURE MAX --- #

# load daymet geoTifs
tmax_dir <- "maps/climate/tmax"
tmax_files <- list.files(path = tmax_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()

# QQ geotifs
tmax_test <- stack("maps/climate/tmax/daymet_v4_tmax_monavg_na_1995.tif")
plot(tmax_test[[6]])
crs(tmax_test)
nlayers(tmax_test)
res(tmax_test)
extent(tmax_test)

# extract summer (jun-aug) avg tmax per year per atlas block
tmax_extract_summer <- function(file_path) { # set up fx
  year <- as.numeric(gsub("\\D", "", basename(file_path)))
  tmax_stack <- raster::stack(file_path)
  tmax_summer_avg <- raster::calc(tmax_stack[[6:8]], fun = mean, na.rm = TRUE)
  tmax_vals <- exact_extract(tmax_summer_avg, blocks_wibba_daymet, 'mean')
  
  data.frame(
    atlas_block = blocks_wibba_daymet$BLOCK_ID,
    year = year,
    tmax_summer = tmax_vals
  )
}

tmax_annual_summer <- map_dfr(tmax_files, tmax_extract_summer) %>% # loop fx to all geotifs
  mutate(
    atlas_period = case_when( # define atlas periods
      year >= 1995 & year <= 2000 ~ "Atlas1",
      year >= 2015 & year <= 2019 ~ "Atlas2",
    )
  ) %>%
  mutate(
    year = as.numeric(stringr::str_extract(as.character(year), "\\d{4}$")) #fix year names
  )
  
# summary tables of results
tmax_block_period_avg <- tmax_annual_summer %>% # tmax by block and atlas period
  group_by(atlas_block, atlas_period) %>%
  summarise(
    tmax_avg_summer = mean(tmax_summer, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

tmax_block_wide <- tmax_block_period_avg %>% # pivot wide, calc differences
  pivot_wider(
    names_from = atlas_period,
    values_from = tmax_avg_summer,
    names_glue = "tmax_summer_{atlas_period}"
  ) %>%
  mutate(
    tmax_summer_diff = tmax_summer_Atlas2 - tmax_summer_Atlas1,
    tmax_summer_diff_z = scale(tmax_summer_diff)[, 1] # z-standardize
  )

# join to summary dfs
wibba_summary <- wibba_summary %>%
  left_join(
    dplyr::select(tmax_block_wide, atlas_block, tmax_summer_diff_z),
    by = "atlas_block"    
  )

rcki_wibba_summary <- rcki_wibba_summary %>% # join to species model dfs
  left_join(
    dplyr::select(tmax_block_wide, atlas_block, tmax_summer_diff_z),
    by = "atlas_block"
  )  


# --- PRECIPITATION --- #

# load daymet geoTifs
prcp_dir <- "maps/climate/prcp"
prcp_files <- list.files(path = prcp_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()

# QQ geotifs
prcp_test <- stack("maps/climate/prcp/daymet_v4_prcp_monttl_na_1995.tif")
plot(prcp_test[[6]])
crs(prcp_test)
nlayers(prcp_test)
res(prcp_test)
extent(prcp_test)

# extract summer (jun-aug) aprcp total per atlas block
prcp_extract_summer <- function(file_path) { # set up fx
  year <- as.numeric(gsub("\\D", "", basename(file_path)))
  prcp_stack <- raster::stack(file_path)
  prcp_summer_avg <- raster::calc(prcp_stack[[6:8]], fun = mean, na.rm = TRUE)
  prcp_vals <- exact_extract(prcp_summer_avg, blocks_wibba_daymet, 'mean')
  
  data.frame(
    atlas_block = blocks_wibba_daymet$BLOCK_ID,
    year = year,
    prcp_summer = prcp_vals
  )
}

prcp_annual_summer <- map_dfr(prcp_files, prcp_extract_summer) %>% # loop fx to all geotifs
  mutate(
    atlas_period = case_when( # define atlas periods
      year >= 1995 & year <= 2000 ~ "Atlas1",
      year >= 2015 & year <= 2019 ~ "Atlas2",
    )
  ) %>%
  mutate(
    year = as.numeric(stringr::str_extract(as.character(year), "\\d{4}$")) #fix year names
  )

# create tables of results
prcp_block_period_avg <- prcp_annual_summer %>% # prcp by block and atlas period
  group_by(atlas_block, atlas_period) %>%
  summarise(
    prcp_avg_summer = mean(prcp_summer, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

prcp_block_wide <- prcp_block_period_avg %>% # pivot wider, calc diffs
  pivot_wider(
    names_from = atlas_period,
    values_from = prcp_avg_summer,
    names_glue = "prcp_summer_{atlas_period}"
  ) %>%
  mutate(
    prcp_summer_diff = prcp_summer_Atlas2 - prcp_summer_Atlas1,
    prcp_summer_diff_z = scale(prcp_summer_diff)[, 1] # z-standardize
  )

#join to df
wibba_summary <- wibba_summary %>%
  left_join(
    dplyr::select(prcp_block_wide, atlas_block, prcp_summer_diff_z),
    by = "atlas_block"    
  )

rcki_wibba_summary <- rcki_wibba_summary %>% # join to species model dfs
  left_join(
    dplyr::select(prcp_block_wide, atlas_block, prcp_summer_diff_z),
    by = "atlas_block"    
  )