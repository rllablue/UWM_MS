#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

library(dplyr)
library(tidyr)
library(magrittr)
library(lubridate)
library(purrr)
library(ggplot2)
library(sf)
library(raster)
library(terra)
library(readxl)
library(exactextractr)


################################
### FILTER COMPARABLE BLOCKS ###
################################

# Filtering map-based objects

### --- ALL BLOCKS --- ###

# Load, create maps
blocks_all_shp <- st_read("data/maps/wibba blocks/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # load full block map
  rename(atlas_block = BLOCK_ID)


# --- RLL COMP BLOCKS --- #

# BLOCK MAP #
blocks_rll <- wibba_summary_comp$atlas_block
blocks_rll <- data.frame(atlas_block = blocks_comp) # ?

blocks_rll_shp <- blocks_all_shp %>% # create sf object of filtered comparable blocks
  filter(atlas_block %in% blocks_rll)
st_write(blocks_rll_shp, "outputs/maps/blocks_rll.shp", delete_dsn = TRUE) # create shp file for comp blocks

# BLOCKS OUTLINES #
atlas_outline_crs3071 <- blocks_all_shp %>% # natural crs
  st_union() %>%      # dissolve all blocks into one geometry
  st_as_sf()          # convert back to sf object
st_write(atlas_outline_crs3071, "outputs/maps/atlas_outline_crs3071.shp", append = FALSE)

# PLOTS #
# Visualize all Atlas blocks
ggplot(blocks_all_shp) + # all blocks
  geom_sf(fill = "orange", color = "white") +
  ggtitle("All Atlas Blocks")

# Visualize surveyed v un-surveyed blocks
blocks_surveyed <- checklist_summary %>% 
  dplyr::select(atlas_block, checklists_Atlas1, checklists_Atlas2) %>%
  pivot_longer(
    cols = starts_with("checklists_"),
    names_to = "atlas_period",
    names_prefix = "checklists_",
    values_to = "n_checklists"
  ) %>%
  mutate(surveyed = n_checklists > 0)

blocks_surveyed_shp <- blocks_all_shp %>%
  left_join(blocks_surveyed, by = "atlas_block")

ggplot(blocks_surveyed_shp) +
  geom_sf(aes(fill = surveyed), color = NA) +
  facet_wrap(~atlas_period) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgray")) +
  labs(
    title = "Surveyed Atlas Blocks by Atlas Period",
    fill = "Surveyed"
  ) +
  theme_minimal()

# Visualize comparable Atlas blocks (surveyed in both Atlas periods) (N = 3337)
ggplot(blocks_rll_shp) + 
  geom_sf(fill = "orange", color = "white") +
  ggtitle("RLL Comparable Atlas Blocks")



### --- DNR COMP BLOCKS --- ###

blocks_dnr <- read_xlsx("data/summaries/CompBlocks_DNR2023.xlsx") # df
blocks_dnr <- blocks_dnr$atlas_block # vector

all(blocks_dnr$atlas_block %in% blocks_rll$atlas_block) # ?



# Shp object
blocks_dnr_shp <- blocks_all_shp %>% # create sf object of filtered DNR comparable blocks
  filter(atlas_block %in% blocks_dnr)
st_write(blocks_dnr_shp, "outputs/maps/blocks_dnr.shp", delete_dsn = TRUE) # create shp file for comp blocks

# Visualize DNR comparable Atlas blocks (N = 858)
ggplot(blocks_dnr_shp) + 
  geom_sf(fill = "orange", color = "white") +
  ggtitle("DNR Comparable Atlas Blocks")
