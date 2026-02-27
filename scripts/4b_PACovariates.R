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
blocks_all_sf <- st_transform(blocks_all_sf, 5070) # EPSG:5070 uses m as unit
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
  filter(Own_Name == "SDNR") %>%
  count(Des_Tp, sort = TRUE)

wipad_sf %>% # w/in NGO
  filter(
    Own_Name == "NGO",
    Loc_Own == "National Audubon Society" # for NGO: The Nature Conservancy or National Audubon Society
  ) %>%
  count(Des_Tp, sort = TRUE)





### --- TOTAL PA --- ###
### Calculate PA/block w/ dissolving ie. "full" PA cover per block

wipad_diss_sf <- wipad_sf %>% 
  st_union() %>% 
  st_as_sf() %>%
  st_set_crs(5070)
st_crs(wipad_diss_sf)


# Compute PA per block
blocks_all_sf <- blocks_all_sf %>%
  mutate(
    block_area_km2 = as.numeric(st_area(.) / 1e6) # m2 / 1e6 = km2
  )


# Combine, summarize 
wipad_wibba_summary <- blocks_all_sf %>%
  st_intersection(wipad_diss_sf) %>%
  mutate(overlap_km2 = as.numeric(st_area(.) / 1e6)) %>%
  group_by(atlas_block) %>%
  summarise(pa_area_km2 = sum(overlap_km2), .groups = "drop") %>%
  right_join(
    blocks_all_sf %>%
      st_drop_geometry() %>%
      dplyr::select(atlas_block, block_area_km2),
    by = "atlas_block"
  )  %>%
  mutate(
    pa_area_km2 = replace_na(pa_area_km2, 0),
    pa_percent = (pa_area_km2 / block_area_km2) * 100
  ) %>%
  dplyr::select(atlas_block, pa_area_km2, block_area_km2, pa_percent)
  





### -- FILTERED PAD --- ###
# Calculate PA/block w/ no dissolving (allows for land owner filtering)


# Define owner logic once
wipad_sf <- wipad_sf %>%
  mutate(
    owner_id = case_when(
      Own_Name == "SDNR" ~ "sdnr",
      Own_Name == "USFS" ~ "usfs",
      Own_Name == "NGO" & Loc_Own == "The Nature Conservancy" ~ "tnc",
      Own_Name == "NGO" & Loc_Own == "National Audubon Society" ~ "nas",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(owner_id))


# Compute per-owner PA per block
wipad_own_summary <- blocks_all_sf %>%
 
   # keep only blocks of interest
  filter(atlas_block %in% blocks_rll) %>%
  st_intersection(wipad_sf) %>%
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  group_by(atlas_block, owner_id) %>%
  summarise(area_km2 = sum(area_km2), .groups = "drop") %>%
  st_drop_geometry() %>%
 
   # pivot so each owner is a separate column
  pivot_wider(
    names_from = owner_id,
    values_from = area_km2,
    names_glue = "{owner_id}_area_km2",
    values_fill = 0
  ) %>%
  
  # ensure **all blocks_rll** are included
  right_join(
    tibble(atlas_block = blocks_rll),
    by = "atlas_block"
  ) %>%
  
  # replace NAs with 0 for any owner
  mutate(
    sdnr_area_km2 = replace_na(sdnr_area_km2, 0),
    usfs_area_km2 = replace_na(usfs_area_km2, 0),
    tnc_area_km2  = replace_na(tnc_area_km2, 0),
    nas_area_km2  = replace_na(nas_area_km2, 0)
  ) %>%
  
  # add block total area
  left_join(
    blocks_all_sf %>%
      st_drop_geometry() %>%
      dplyr::select(atlas_block, block_area_km2),
    by = "atlas_block"
  ) %>%
  
  # compute percentages
  mutate(
    sdnr_percent = sdnr_area_km2 / block_area_km2 * 100,
    usfs_percent = usfs_area_km2 / block_area_km2 * 100,
    tnc_percent  = tnc_area_km2  / block_area_km2 * 100,
    nas_percent  = nas_area_km2  / block_area_km2 * 100
  ) %>%
  dplyr::select(
    atlas_block,
    sdnr_area_km2, sdnr_percent,
    usfs_area_km2, usfs_percent,
    tnc_area_km2, tnc_percent,
    nas_area_km2, nas_percent
  )


### --- JOIN TO DATA FRAMES --- ###

# Join to full wibba, wipad summary
wipad_wibba_summary <- wipad_wibba_summary %>%
  left_join(wipad_own_summary, by = "atlas_block") %>%
  mutate(across(ends_with(c("_area_km2", "_percent")), ~replace_na(.x, 0)))


# Join % values to raw covariates df
# Add to raw modeling data
covars_raw_rll <- covars_raw_rll %>%
  left_join(
    wipad_wibba_summary %>%
      dplyr::select(atlas_block, pa_percent, sdnr_percent, usfs_percent, nas_percent, tnc_percent),
    by = "atlas_block"
  )

# write.csv(covars_raw_rll, "data/summaries/covars_raw_rll.csv", row.names = FALSE)
