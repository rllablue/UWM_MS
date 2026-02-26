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
# Atlas blocks geometry
blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>%
  rename(atlas_block = BLOCK_ID)
blocks_all_sf <- st_transform(blocks_all_sf, 5070)
st_crs(blocks_all_sf)

# PAD geometry
st_layers("data/maps/wipad/PADUS4_1_State_WI_GDB_KMZ/PADUS4_1_StateWI.gdb")
wipad_sf <- st_read(
  "data/maps/wipad/PADUS4_1_State_WI_GDB_KMZ/PADUS4_1_StateWI.gdb",
  layer = "PADUS4_1Fee_State_WI"
)
names(wipad_sf)

wipad_sf <- st_transform(wipad_sf, 5070)
st_crs(wipad_sf)


# Filtering for counts
wipad_sf %>%
  filter(Own_Name == "NGO") %>%
  count(Des_Tp, sort = TRUE)

wipad_sf %>% # w/in NGO
  filter(
    Own_Name == "NGO",
    Loc_Own == "National Audubon Society" # or National Audubon Society
  ) %>%
  count(Des_Tp, sort = TRUE)








### --- TOTAL PA --- ###
### Calculate PA/block w/ dissolving ie. "full" PA cover per block

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