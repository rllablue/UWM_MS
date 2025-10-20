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


### --- CRS ALIGNMENT --- ###

# Set up: crs integration of daymet, atlas geodata
daymet_crs <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" # Daymet custom crs (Lambert Conformal Conic)
blocks_comp_daymet_shp <- st_transform(blocks_comp_shp, crs = daymet_crs)


###################
### DAYMET DATA ###
###################

### --- EXTRACT CLIMATOLOGY --- ###

### Flexible function to extract various climate metrics from Daymet 12-band rasters
# for summer months, May-August. Workflow extracts both the long-term averages of a 
# metric across full reference period (38 years) as well as the averages between subset
# reference periods (18 years prior to the end of each Atlas period) for analysis
# of climate anomalies. 

Extract_climate <- function(files_list, blocks_shp, wibba_df, metrics = c("tmax", "tmin", "precip"),
                                     summer_months = 5:8,
                                     period1 = 1982:2000, period2 = 2001:2019) {
  
  wibba_out <- wibba_df
  
  for(metric in metrics) {
    files <- files_list[[metric]]
    
    Extract_annual <- function(file_path) { 
      year <- as.numeric(str_extract(basename(file_path), "\\d{4}"))
      stack_r <- raster::stack(file_path)
      
      if(nlayers(stack_r) < max(summer_months)) {
        warning(paste("Skipping", basename(file_path), "- not enough bands"))
        return(NULL)
      }
      
      summer_avg <- raster::calc(stack_r[[summer_months]], mean, na.rm = TRUE)
      values <- exact_extract(summer_avg, blocks_shp, 'mean')
      
      # Return tibble with dynamic column name
      tibble(
        atlas_block = blocks_shp$atlas_block,
        year = year,
        !!metric := values
      )
    }
    
    annual_df <- map_dfr(files, Extract_annual) # extract annual values
    
    longterm_df <- annual_df %>% # calc long-term raw, z-standardized avg per block
      group_by(atlas_block) %>%
      summarise(
        !!paste0(metric, "_38yr") := mean(.data[[metric]], na.rm = TRUE),
        n_years = n(),
        .groups = "drop"
      ) %>%
      mutate(
        !!paste0(metric, "_38yr_z") := as.numeric(scale(.data[[paste0(metric, "_38yr")]]))
      )
    
    period_df <- annual_df %>% # calc period raw, z-standardized avgs per block
      mutate(
        period = case_when(
          year %in% period1 ~ "t1",
          year %in% period2 ~ "t2",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(period)) %>%
      group_by(atlas_block, period) %>%
      summarise(
        !!paste0(metric, "_avg") := mean(.data[[metric]], na.rm = TRUE),
        n_years = n(),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = period,
        values_from = !!paste0(metric, "_avg"),
        names_glue = paste0(metric, "_{period}")
      ) %>%
      mutate(
        !!paste0(metric, "_diff") := .data[[paste0(metric, "_t2")]] - .data[[paste0(metric, "_t1")]],
        !!paste0(metric, "_diff_z") := as.numeric(scale(.data[[paste0(metric, "_diff")]]))
      )

    wibba_out <- wibba_out %>% # join values to modeling df
      left_join(dplyr::select(longterm_df, atlas_block, !!paste0(metric, "_38yr_z")), by = "atlas_block") %>%
      left_join(dplyr::select(period_df, atlas_block, !!paste0(metric, "_diff_z")), by = "atlas_block")
  }
  
  return(wibba_out)
}

# Set locations for files
tmax_dir <- "data/maps/daymet/38yr_tmax"
tmin_dir <- "data/maps/daymet/38yr_tmin"
precip_dir <- "data/maps/daymet/38yr_precip"

tmax_files <- list.files(path = tmax_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()
tmin_files <- list.files(path = tmin_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()
precip_files <- list.files(path = precip_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()

files_list <- list( # put in list for function to apply
  tmax = tmax_files,
  tmin = tmin_files,
  precip = precip_files
)

# Apply function, join to modeling df
wibba_modeling_comp <- Extract_climate(
  files_list = files_list,
  blocks_shp = blocks_comp_daymet_shp,
  wibba_df = wibba_modeling_comp,
  metrics = c("tmax", "tmin", "precip")
)


















### --- TEMPERATURE MAX (TMAX) --- ###

# Load daymet geoTifs
tmax_dir <- "data/maps/daymet/38yr_tmax"
tmax_files <- list.files(path = tmax_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()

# Extract summer (May-Aug) avg tmax per year per atlas block
Extract_tmax <- function(file_path) {
  year <- as.numeric(str_extract(basename(file_path), "\\d{4}")) # derive year from file name
  tmax_stack <- raster::stack(file_path)
  tmax_summer_avg <- raster::calc(tmax_stack[[5:8]], fun = mean, na.rm = TRUE) # subset May-Aug from 12 monthly bands
  tmax_vals <- exact_extract(tmax_summer_avg, blocks_comp_daymet_shp, 'mean') # extract per block
  
  data.frame(
    atlas_block = blocks_comp_daymet_shp$atlas_block,
    year = year,
    tmax_summer = tmax_vals
  )
}

# Apply across 38 year reference period
tmax_avg_annual <- map_dfr(tmax_files, Extract_tmax) # Max summer temps by year, block


# Long-term tmax avg
tmax_avg_longterm <- tmax_avg_annual %>% 
  group_by(atlas_block) %>%
  summarise(
    tmax_summer_38yr = mean(tmax_summer, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  ) %>%
  mutate(
    tmax_summer_38yr_z = as.numeric(scale(tmax_summer_38yr))
  )

# Between-period tmax avgs (19 year reference periods for each Atlas)
tmax_block_periods <- tmax_annual %>%
  mutate(
    atlas_tperiod = case_when(
      year >= 1982 & year <= 2000 ~ "t1_Atlas1",
      year >= 2001 & year <= 2019 ~ "t2_Atlas2",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(atlas_block, atlas_tperiod) %>%
  summarise(
    tmax_summer_avg = mean(tmax_summer, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

tmax_block_wide <- tmax_block_periods %>% # pivot wide, compute tmax differences between periods
  pivot_wider(
    names_from = atlas_period,
    values_from = tmax_summer_avg,
    names_glue = "tmax_summer_{atlas_period}"
  ) %>%
  mutate(
    tmax_summer_diff = tmax_summer_Atlas2 - tmax_summer_Atlas1,
    tmax_summer_diff_z = as.numeric(scale(tmax_summer_diff))
  )

# Join to summary dfs
wibba_modeling_comp <- wibba_modeling_comp %>%
  left_join(
    dplyr::select(tmax_avg_longterm, atlas_block, tmax_summer_38yr_z),
    by = "atlas_block"    
  )



### --- PRECIPITATION --- ###

# load daymet geoTifs
prcp_dir <- "data/maps/daymet/38yr_prcp"
prcp_files <- list.files(path = prcp_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()


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
