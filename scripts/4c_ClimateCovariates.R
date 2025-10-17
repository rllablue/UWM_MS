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



### --- CLIMATE DATA --- ###

### CLIMATE PROXY ###
# Longitude, Latitude centroids from atlas blocks as climate proxies

# Real and z-standard lat and lon values 
# z-standard. for mock-up model/climate effects proxy / will use real for spatial autocorrection
wibba_summary_comp <- wibba_summary_comp %>%
  left_join(
    blocks_comp_shp_5070 %>%
      mutate(centroid = st_centroid(geometry)) %>% # compute centroids
      mutate(lon = st_coordinates(centroid)[,1], # extract centroid longitude
             lat = st_coordinates(centroid)[,2]) %>% # extract centroid latitude
      st_drop_geometry() %>% # drop geometry for join
      dplyr::select(atlas_block, lon, lat), # keep only needed cols
    by = "atlas_block"
  ) %>%
  mutate(
    lon_z = as.numeric(scale(lon)), # z-standardized longitude
    lat_z = as.numeric(scale(lat)) # z-standardized latitude
  )


### --- DAYMET DATA --- ###

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

