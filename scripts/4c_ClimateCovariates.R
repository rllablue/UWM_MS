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
library(furrr)
library(future)


### --- CRS ALIGNMENT --- ###

# Set up: crs integration of daymet, atlas geodata
# Daymet custom crs (Lambert Conformal Conic)
blocks_comp_daymet_shp <- st_transform(blocks_comp_shp, crs = daymet_crs)
daymet_crs <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"


###################
### DAYMET DATA ###
###################

### --- EXTRACT CLIMATOLOGY --- ###

### Flexible function to extract various climate metrics from Daymet 12-band rasters
# for summer months, May-August. Workflow extracts both the long-term averages of a 
# metric across full reference period (38 years) as well as the averages between subset
# reference periods (18 years prior to the end of each Atlas period) for analysis
# of climate anomalies. 

Extract_climate <- function(files_list, 
                            blocks_shp, 
                            wibba_df, 
                            metrics = c("tmax", "tmin", "prcp"),
                            summer_months = 5:8,
                            period1 = 1982:2000, 
                            period2 = 2001:2019) {
  
  wibba_out <- dplyr::select(wibba_df, atlas_block)
  
  for (metric in metrics) {
    files <- files_list[[metric]]
    
    # Safeguard: drop missing or invalid file paths
    files <- purrr::compact(files)
    files <- files[file.exists(files)]
    if (length(files) == 0) {
      warning(paste("No valid files found for", metric))
      next
    }
    
    Extract_annual <- function(file_path) {
      year <- as.numeric(stringr::str_extract(basename(file_path), "\\d{4}"))
      stack_r <- raster::stack(file_path)
      
      if (nlayers(stack_r) < max(summer_months)) {
        warning(paste("Skipping", basename(file_path), "- not enough bands"))
        return(NULL)
      }
      
      # Mean of May–Aug bands
      summer_avg <- raster::calc(stack_r[[summer_months]], mean, na.rm = TRUE)
      
      # Extract mean per polygon
      values <- exactextractr::exact_extract(summer_avg, blocks_shp, 'mean')
      values <- unlist(values)
      
      tibble::tibble(
        atlas_block = blocks_shp$atlas_block,
        year = year,
        !!metric := values
      )
    }
    
    # Annual summer values per year per block
    annual_df <- purrr::map_dfr(files, Extract_annual)
    
    if (!"atlas_block" %in% colnames(annual_df)) {
      stop("atlas_block column is missing in annual_df.")
    }
    
    # Long-term mean (38-year reference)
    longterm_df <- annual_df %>%
      dplyr::group_by(atlas_block) %>%
      dplyr::summarise(
        !!paste0(metric, "_38yr") := mean(.data[[metric]], na.rm = TRUE),
        n_years = dplyr::n(),
        .groups = "drop"
      )
    
    # Period means (18 years each)
    period_df <- annual_df %>%
      dplyr::mutate(
        period = dplyr::case_when(
          year %in% period1 ~ "t1",
          year %in% period2 ~ "t2",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(period)) %>%
      dplyr::group_by(atlas_block, period) %>%
      dplyr::summarise(
        !!paste0(metric, "_avg") := mean(.data[[metric]], na.rm = TRUE),
        n_years = dplyr::n(),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = period,
        values_from = !!paste0(metric, "_avg"),
        names_glue = paste0(metric, "_{period}")
      ) %>%
      dplyr::mutate(
        !!paste0(metric, "_diff") := .data[[paste0(metric, "_t2")]] - .data[[paste0(metric, "_t1")]],
      )
    
    # Join to modeling df
    wibba_out <- wibba_out %>%
      dplyr::left_join(longterm_df, by = "atlas_block") %>%
      dplyr::left_join(period_df, by = "atlas_block")
  }
  
  return(wibba_out)
}


### --- FILE SETUP --- ###
tmax_dir <- "data/maps/daymet/38yr_tmax"
tmin_dir <- "data/maps/daymet/38yr_tmin"
prcp_dir <- "data/maps/daymet/38yr_prcp"

tmax_files <- list.files(path = tmax_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()
tmin_files <- list.files(path = tmin_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()
prcp_files <- list.files(path = prcp_dir, pattern = "\\.tif$", full.names = TRUE) %>% sort()

files_list <- list(  # list for function
  tmax = tmax_files,
  tmin = tmin_files,
  prcp = prcp_files
)


### --- APPLY --- ###

# PARALLELIZE PROCESS #
# Set up
future::plan(multisession, workers = 10)

# Assign process to run in parallel
Parallel_climate <- future({
  Extract_climate(
    files_list = files_list,
    blocks_shp = blocks_comp_daymet_shp,
    wibba_df = wibba_modeling_covars,
    metrics = c("tmax", "tmin", "prcp")
  )
})

# Check completion
resolved(Parallel_climate)


# SUMMARIZE #
# Join all metrics to summary df
climate_summary_covars <- value(Parallel_climate)

# Join to modeling df
wibba_modeling_covars <- wibba_modeling_covars %>%
  left_join(climate_summary_covars, by = "atlas_block")



###############
### PROXIES ###
###############

### --- CLIMATE PROXY --- ###
# Longitude, Latitude centroids from atlas blocks as climate proxies
# Proxies to be used in SR residual regression as climate/geographic gradient proxy

wibba_modeling_covars <- wibba_modeling_covars %>%
  left_join(
    blocks_comp_shp_5070 %>%
      st_centroid() %>%  # compute centroids
      dplyr::mutate(
        lon = st_coordinates(.)[,1], # raw longitude
        lat = st_coordinates(.)[,2]  # raw latitude
      ) %>%
      st_drop_geometry() %>%
      dplyr::select(atlas_block, lon, lat),  # keep raw coords only
    by = "atlas_block"
  )