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



############################
### PROTECTED AREA (PAD) ###
############################

### --- UPLOAD FILES --- ### 
blocks_comp_shp_5070 <- blocks_comp_shp %>% # reproject atlas blocks
  st_transform(5070)  # NAD83, Conus Albers EPSG:5070 (ideal for area calcs)

pad_wi_5070 <- st_read("data/maps/pad wi/pad_wi_polygons.shp") %>% # align pad crs
  st_transform(5070)


### --- FULL PAD --- ###

# Calculate PA/block w/ dissolving ie. "full" PA cover per block
pad_wi_diss <- pad_wi_5070 %>% 
  st_union() %>% 
  st_sf(geometry = .) %>% # attempt to remedy inaccurate pa %s over 100
  st_set_crs(5070)

# Summarize PA by blocks
blocks_pad_summary_diss <- st_intersection(blocks_comp_shp_5070, pad_wi_diss) %>%
  mutate(overlap_area = st_area(.)) %>%
  group_by(atlas_block) %>%
  summarise(pa_area = sum(overlap_area), .groups = "drop") %>%
  right_join(
    blocks_comp_shp_5070 %>%
      mutate(block_area = st_area(.)) %>%
      st_drop_geometry() %>%
      dplyr::select(atlas_block, block_area),
    by = "atlas_block"
  ) %>%
  mutate(
    pa_area = if_else(is.na(pa_area), units::set_units(0, m^2), pa_area),
    block_area_km2 = units::set_units(block_area, "m^2") %>% units::set_units("km^2"),
    pa_area_km2 = units::set_units(pa_area, "m^2") %>% units::set_units("km^2"),
    pa_percent = as.numeric(pa_area_km2 / block_area_km2 * 100)
  ) %>%
  dplyr::select(atlas_block, pa_area_km2, block_area_km2, pa_percent)


# Join to modeling df
wibba_covars_raw <- wibba_covars_raw %>% 
  left_join(
    blocks_pad_summary_diss %>% dplyr::select(atlas_block, pa_percent),
    by = "atlas_block"
  ) %>%
  dplyr::select(-geometry)













### -- FILTERED PAD --- ###
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