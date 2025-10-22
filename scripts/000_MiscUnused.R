########################
### MISC CODE CHUNKS ###
########################

### --- CLIMATE PROXY --- ###
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



### --- ELEVATION DATA (NED) --- ###
# Mean elevation per block done in ArcPro (1/3 arc-sec DEMs,~10m)

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






## GEOGRAPHIC FILTER ## (UNNEEDED FOR WIBBA)
# geographic filtering, clean-up
wibba_grid <- read_sf("maps/wibba blocks/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # read in atlas grid shape file
  st_transform(4269) # convert to NAD83, crs of shp file

checklists_sf <- checklists %>%
  select(checklist_id, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% # convert checklist locations to points
  st_transform(4269) # match atlas block crs

checklists_joined <- st_join(checklists_sf, wibba_grid[, c("BLOCK_ID")], left = TRUE) # create df with checklist id, lat/lon, and BLOCK_ID

in_region <- checklists_joined %>%
  filter(!is.na(BLOCK_ID)) # keep only checklists in joined atlas block shp

checklists_wip <- checklists %>%
  semi_join(st_drop_geometry(in_region), by = "checklist_id") %>% # remove checklists outside atlas blocks
  left_join(st_drop_geometry(in_region[, c("checklist_id", "BLOCK_ID")]), by = "checklist_id") %>% # adds BLOCK_ID from shp to checklist data w/o overwriting
  mutate(
    atlas_block = if_else(is.na(atlas_block), BLOCK_ID, atlas_block) # fills in missing atlas blocks w/o overwriting
  ) %>%
  select(-BLOCK_ID) #remove BLOCK_ID column

observations_wip <- observations %>%
  semi_join(checklists, by = "checklist_id") %>%  # include only matching checklist_id from checklists df
  left_join(
    checklists %>% 
      select(checklist_id, atlas_block) %>%  # join atlas_block to observations using checklist_id
      rename(joined_block = atlas_block),
    by = "checklist_id"
  ) %>%
  mutate(
    old_block = atlas_block, # preserve existing
    atlas_block = joined_block # overwrite with cleaned value
  ) %>%
  select(-joined_block) # drop temp join column
