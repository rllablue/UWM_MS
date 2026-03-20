#############
### SETUP ###
#############

### --- LOAD PACKAGES --- ###

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,
  tidyr,
  lubridate,
  purrr,
  ggplot2,
  sf,
  raster,
  terra,
  exactextractr)


### --- GEOMETRY FILES --- ### 

# Atlas Block Geometry #
# Native CRS: WI Transverse Mercator, EPSG:3071

blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>%
  rename(atlas_block = BLOCK_ID)
crs(blocks_all_sf)

blocks_rll_sf <- blocks_all_sf %>%
  filter(atlas_block %in% blocks_rll)
crs(blocks_rll_sf)

blocks_dnr_sf <- blocks_all_sf %>%
  filter(atlas_block %in% blocks_dnr)
crs(blocks_dnr_sf)


# transformation, CRS: Conus Albers, NAD83, EPSG:5070
blocks_all_sf <- st_transform(wipad_sf, 5070)
blocks_rll_sf <- st_transform(wipad_sf, 5070)
blocks_dnr_sf <- st_transform(wipad_sf, 5070)

st_crs(blocks_all_sf)
st_crs(blocks_rll_sf)
st_crs(blocks_dnr_sf)



# PAD Geometry #
# Native CRS: USA Contiguous Albers Equal Area Conic [USGS], NAD83

st_layers("data/maps/wipad/PADUS4_1_State_WI_GDB_KMZ/PADUS4_1_StateWI.gdb")
wipad_sf <- st_read(
  "data/maps/wipad/PADUS4_1_State_WI_GDB_KMZ/PADUS4_1_StateWI.gdb",
  layer = "PADUS4_1Fee_State_WI"
)
names(wipad_sf)
st_crs(wipad_sf)


# transformation, CRS: Conus Albers, NAD83, EPSG:5070
wipad_sf <- st_transform(wipad_sf, 5070)
st_crs(wipad_sf)





### --- SUMMARIES --- ###

# Counts #
wipad_sf %>%
  filter(Own_Name == "SDNR") %>%
  count(Des_Tp, sort = TRUE)

wipad_sf %>%
  filter(
    Own_Name == "USFS",
    GAP_Sts == "4"
    ) %>%
  count(Des_Tp, sort = TRUE)

wipad_sf %>% # w/in NGO
  filter(
    Own_Name == "NGO",
    Loc_Own == "National Audubon Society" # for NGO: The Nature Conservancy or National Audubon Society
  ) %>%
  count(Des_Tp, sort = TRUE)





### --- PA CALCULATIONS --- ###

# Total PA Coverage #
# w/ dissolving, i.e. total regardless of ownership

# geometry: dissolve all PA boundaries
wipad_diss_sf <- wipad_sf %>% 
  st_union() %>% 
  st_as_sf() %>%
  st_set_crs(5070)
st_crs(wipad_diss_sf)

# compute: PA per block
blocks_all_sf <- blocks_all_sf %>%
  mutate(
    block_area_km2 = as.numeric(st_area(.) / 1e6) # m2 / 1e6 = km2
  )

# collect: PA values  
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
    pa_prop = pa_area_km2 / block_area_km2 # total pa proportion of block
  ) %>%
  dplyr::select(atlas_block, pa_area_km2, block_area_km2, pa_prop)

  



# PA Coverage by Sub-type #

# PA Ownership #
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


# Compute owner area per block
wipad_own_summary <- blocks_all_sf %>%
  filter(atlas_block %in% blocks_rll) %>%
  st_intersection(wipad_sf) %>%
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  group_by(atlas_block, owner_id) %>%
  summarise(area_km2 = sum(area_km2), .groups = "drop") %>%
  st_drop_geometry() %>%
  pivot_wider(
    names_from = owner_id,
    values_from = area_km2,
    names_glue = "{owner_id}_area_km2",
    values_fill = 0
  ) %>%
  right_join(
    tibble(atlas_block = blocks_rll),
    by = "atlas_block"
  ) %>%
  mutate(
    sdnr_area_km2 = replace_na(sdnr_area_km2, 0),
    usfs_area_km2 = replace_na(usfs_area_km2, 0),
    tnc_area_km2  = replace_na(tnc_area_km2, 0),
    nas_area_km2  = replace_na(nas_area_km2, 0)
  ) %>%
  left_join(
    blocks_all_sf %>%
      st_drop_geometry() %>%
      dplyr::select(atlas_block, block_area_km2),
    by = "atlas_block"
  ) %>%
  
  # additive intensity = area / block area
  mutate(
    sdnr_prop = sdnr_area_km2 / block_area_km2,
    usfs_prop = usfs_area_km2 / block_area_km2,
    tnc_prop  = tnc_area_km2  / block_area_km2,
    nas_prop  = nas_area_km2  / block_area_km2
  ) %>%
  dplyr::select(
    atlas_block,
    sdnr_prop, usfs_prop, tnc_prop, nas_prop
  )





# PA Protection Stringency #
# GAP Status

wipad_gap_summary <- blocks_all_sf %>%
  filter(atlas_block %in% blocks_rll) %>%
  st_intersection(wipad_sf) %>%
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  group_by(atlas_block, GAP_Sts) %>%
  summarise(area_km2 = sum(area_km2), .groups = "drop") %>%
  st_drop_geometry() %>%
  pivot_wider(
    names_from = GAP_Sts,
    values_from = area_km2,
    names_glue = "gap{GAP_Sts}_area_km2",
    values_fill = 0
  ) %>%
  right_join(
    tibble(atlas_block = blocks_rll),
    by = "atlas_block"
  ) %>%
  mutate(
    gap1_area_km2 = replace_na(gap1_area_km2, 0),
    gap2_area_km2 = replace_na(gap2_area_km2, 0),
    gap3_area_km2 = replace_na(gap3_area_km2, 0),
    gap4_area_km2 = replace_na(gap4_area_km2, 0)
  ) %>%
  left_join(
    blocks_all_sf %>%
      st_drop_geometry() %>%
      dplyr::select(atlas_block, block_area_km2),
    by = "atlas_block"
  ) %>%
  
  # additive intensity = area / block area
  mutate(
    gap1_prop = gap1_area_km2 / block_area_km2,
    gap2_prop = gap2_area_km2 / block_area_km2,
    gap3_prop = gap3_area_km2 / block_area_km2,
    gap4_prop = gap4_area_km2 / block_area_km2
  ) %>%
  dplyr::select(
    atlas_block,
    gap1_prop, gap2_prop, gap3_prop, gap4_prop
  )



# Summarize, Join to Modeling Data
wipad_wibba_summary <- wipad_wibba_summary %>%
  left_join(wipad_own_summary, by = "atlas_block") %>%
  left_join(wipad_gap_summary, by = "atlas_block") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

covars_raw_rll <- covars_raw_rll %>%
  left_join(
    wipad_wibba_summary %>%
      dplyr::select(atlas_block, pa_prop,
                    sdnr_prop, usfs_prop, tnc_prop, nas_prop,
                    gap1_prop, gap2_prop, gap3_prop, gap4_prop),
    by = "atlas_block"
  ) %>%
  dplyr::select(-geometry)

# write.csv(covars_raw_rll, "data/summaries/covars_raw_rll.csv", row.names = FALSE)









############################### COMPOSITIONAL STRUCTURE ########################

### --- FILTERED PAD --- ###

### LAND OWNER ###

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

# Compute owner area per block
wipad_own_summary <- blocks_all_sf %>%
  filter(atlas_block %in% blocks_rll) %>%
  st_intersection(wipad_sf) %>%
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  group_by(atlas_block, owner_id) %>%
  summarise(area_km2 = sum(area_km2), .groups = "drop") %>%
  st_drop_geometry() %>%
  pivot_wider(
    names_from = owner_id,
    values_from = area_km2,
    names_glue = "{owner_id}_area_km2",
    values_fill = 0
  ) %>%
  right_join(
    tibble(atlas_block = blocks_rll),
    by = "atlas_block"
  ) %>%
  mutate(
    sdnr_area_km2 = replace_na(sdnr_area_km2, 0),
    usfs_area_km2 = replace_na(usfs_area_km2, 0),
    tnc_area_km2  = replace_na(tnc_area_km2, 0),
    nas_area_km2  = replace_na(nas_area_km2, 0)
  ) %>%
  left_join(
    wipad_wibba_summary %>%
      dplyr::select(atlas_block, pa_area_km2, pa_prop),
    by = "atlas_block"
  ) %>%
  
  # manager proportions of total pa per block
  mutate(
    sdnr_prop = ifelse(pa_prop > 0, sdnr_area_km2 / pa_area_km2, 0),
    usfs_prop = ifelse(pa_prop > 0, usfs_area_km2 / pa_area_km2, 0),
    tnc_prop  = ifelse(pa_prop > 0, tnc_area_km2  / pa_area_km2, 0),
    nas_prop  = ifelse(pa_prop > 0, nas_area_km2  / pa_area_km2, 0)
  ) %>%
  dplyr::select(
    atlas_block,
    sdnr_area_km2, sdnr_prop,
    usfs_area_km2, usfs_prop,
    tnc_area_km2, tnc_prop,
    nas_area_km2, nas_prop
  )





### --- FILTERED PAD: GAP STATUS --- ###

# Compute per-GAP PA per block
wipad_gap_summary <- blocks_all_sf %>%
  
  filter(atlas_block %in% blocks_rll) %>%
  st_intersection(wipad_sf) %>%
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  group_by(atlas_block, GAP_Sts) %>%
  summarise(area_km2 = sum(area_km2), .groups = "drop") %>%
  st_drop_geometry() %>%
  
  # pivot so each GAP status is its own column
  pivot_wider(
    names_from = GAP_Sts,
    values_from = area_km2,
    names_glue = "gap{GAP_Sts}_area_km2",
    values_fill = 0
  ) %>%
  
  # ensure all blocks_rll included
  right_join(
    tibble(atlas_block = blocks_rll),
    by = "atlas_block"
  ) %>%
  
  # replace NAs with 0
  mutate(
    gap1_area_km2 = replace_na(gap1_area_km2, 0),
    gap2_area_km2 = replace_na(gap2_area_km2, 0),
    gap3_area_km2 = replace_na(gap3_area_km2, 0),
    gap4_area_km2 = replace_na(gap4_area_km2, 0)
  ) %>%
  
  # add block total area
  left_join(
    wipad_wibba_summary %>%
      dplyr::select(atlas_block, pa_area_km2, pa_prop),
    by = "atlas_block"
  ) %>%
  
  # compute proportions of total pa 
  mutate(
    gap1_prop = ifelse(pa_prop > 0, gap1_area_km2 / pa_area_km2, 0),
    gap2_prop = ifelse(pa_prop > 0, gap2_area_km2 / pa_area_km2, 0),
    gap3_prop = ifelse(pa_prop > 0, gap3_area_km2 / pa_area_km2, 0)
    # gap4 = reference category
  ) %>%
  dplyr::select(atlas_block, gap1_prop, gap2_prop, gap3_prop)





### --- JOIN TO DATA FRAMES --- ###

# Join to summary df
wipad_wibba_summary <- wipad_wibba_summary %>%
  left_join(wipad_own_summary, by = "atlas_block") %>%
  left_join(wipad_gap_summary, by = "atlas_block") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

# Join covars for modeling into pool
covars_raw_rll <- covars_raw_rll %>%
  left_join(
    wipad_wibba_summary %>%
      dplyr::select(atlas_block, pa_prop,
                    sdnr_prop, usfs_prop, tnc_prop, nas_prop,
                    gap1_prop, gap2_prop, gap3_prop),
    by = "atlas_block"
  ) %>%
  dplyr::select(-geometry)

# write.csv(covars_raw_rll, "data/summaries/covars_raw_rll.csv", row.names = FALSE)
