library(sf)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


### --- BASE FILES --- ###

# Atlas blocks geometry
blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # load full block map
  rename(atlas_block = BLOCK_ID)
crs(blocks_all_sf)

blocks_background <- blocks_all_sf %>%
  filter(atlas_block %in% union(blocks_rll, blocks_dnr))


# PAD geometry
wipad_sf <- st_read("data/maps/wipad/pad_wi_polygons.shp") ######### NO NO NO NO
crs(wipad_sf)
wipad_sf <- st_transform(wipad_sf, 3071)



# PAD geometry
st_layers("data/maps/wipad/PADUS4_1_State_WI_GDB_KMZ/PADUS4_1_StateWI.gdb")
wipad_sf <- st_read(
  "data/maps/wipad/PADUS4_1_State_WI_GDB_KMZ/PADUS4_1_StateWI.gdb",
  layer = "PADUS4_1Fee_State_WI"
)
names(wipad_sf)

wipad_sf <- st_transform(wipad_sf, 5070)
st_crs(wipad_sf)

















# Join RAW PA values and atlas blocks
pa_raw_df <- covars_raw_rll %>%
  dplyr::select(
    atlas_block, pa_prop, 
    usfs_prop, nas_prop, sdnr_prop, tnc_prop,
    gap1_prop, gap2_prop, gap3_prop, gap4_prop
  ) %>% mutate(
    pa_percent = round(pa_prop * 100, 6),
    pa_percent = pmin(pa_percent, 100)
  )

# %>% mutate(pa_percent = pa_prop * 100)
  

blocks_pa_sf <- blocks_all_sf %>%
  left_join(pa_raw_df, by = "atlas_block")

summary(blocks_pa_sf$pa_prop)




### --- Visualizations --- ###

## WIBBA ##

# All blocks
ggplot() +
  geom_sf(data = blocks_all_sf,
          fill = NA,
          color = "grey40",
          linewidth = 0.2) +
  theme_void() +
  labs(title = "Wisconsin Breeding Bird Atlas Blocks")


# RLL blocks
ggplot() +
  geom_sf(data = blocks_rll_sf,
          fill = NA,
          color = "grey40",
          linewidth = 0.2) +
  theme_void() +
  labs(title = "Wisconsin Breeding Bird Atlas Blocks")


# DNR blocks
ggplot() +
  geom_sf(data = blocks_dnr_sf,
          fill = NA,
          color = "grey40",
          linewidth = 0.2) +
  theme_void() +
  labs(title = "Wisconsin Breeding Bird Atlas Blocks")






## PROTECTED AREA ##

# PAD map all
ggplot() +
  geom_sf(data = wipad_sf,
          fill = "darkgreen",
          color = NA,
          alpha = 0.6) +
  theme_void() +
  labs(title = "Protected Areas in Wisconsin")


ggplot() +
  geom_sf(data = wipad_sf,
          fill = "orange",
          color = "black",
          linewidth = 0.2) +
  theme_void() +
  labs(title = "Protected Area Boundaries")


# PAD map by ownership
# Filter for the main ownerships
wipad_selected <- wipad_sf %>%
  filter(
    Own_Name %in% c("SDNR", "USFS") |           # SDNR and USFS directly
      (Own_Name == "NGO" & Loc_Own %in% c("The Nature Conservancy", "National Audubon Society"))
  ) %>%
  mutate(
    Owner_Group = case_when(
      Own_Name == "SDNR" ~ "SDNR",
      Own_Name == "USFS" ~ "USFS",
      Own_Name == "NGO" & Loc_Own == "The Nature Conservancy" ~ "TNC",
      Own_Name == "NGO" & Loc_Own == "National Audubon Society" ~ "NAS",
      TRUE ~ "Other"
    ),
    Owner_Group = factor(Owner_Group, levels = c("SDNR", "TNC", "NAS", "USFS"))
  )

# Define colors for each group
owner_colors <- c(
  "SDNR" = "firebrick",
  "TNC"  = "cornflowerblue",
  "NAS"  = "plum",
  "USFS" = "lightgreen"
)

# Plot
ggplot() +
  geom_sf(data = wipad_selected,
          aes(fill = Owner_Group),
          color = NA,
          alpha = 0.8) +
  geom_sf(data = wipad_sf,
          fill = NA,
          color = NA,
          linewidth = 0.2) +
  scale_fill_manual(
    values = owner_colors,
    name = "Ownership",
    guide = guide_legend(
      label.position = "left",   # moves the text to the left of the swatch
      keywidth = 1.5,            # adjust swatch width
      keyheight = 0.8,
      reverse = FALSE
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 11)
  ) +
  labs(
    title = "Wisconsin Protected Areas",
    subtitle = "Selected Ownerships Highlighted"
  )





# PAD map by gap status
wipad_gap <- wipad_sf %>%
  filter(GAP_Sts %in% c(1, 2, 3, 4)) %>%
  mutate(
    GAP_Sts = factor(GAP_Sts, levels = c(1, 2, 3, 4))
  )

# Define colors (ecologically intuitive gradient)
gap_colors <- c(
  "1" = "#006400",  # dark green (highest protection)
  "2" = "#41ab5d",
  "3" = "#fe9929",
  "4" = "#de2d26"   # red (lowest protection)
)

# Plot
ggplot() +
  geom_sf(data = wipad_gap,
          aes(fill = GAP_Sts),
          color = NA,
          alpha = 0.85) +
  
  # Outline all PAD lands
  geom_sf(data = wipad_sf,
          fill = NA,
          color = NA,
          linewidth = 0.2) +
  
  scale_fill_manual(
    values = gap_colors,
    name = "GAP Status",
    guide = guide_legend(
      label.position = "left",   # moves the text to the left of the swatch
      keywidth = 1.5,            # adjust swatch width
      keyheight = 0.8,
      reverse = FALSE
    )
  ) +
  
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 11)
  ) +
  
  labs(
    title = "Wisconsin Protected Areas",
    subtitle = "Colored by GAP Status (1 = Highest Protection)"
  )








# PAD-WIBBA Overlap
ggplot() +
  geom_sf(data = blocks_all_sf,
          fill = NA,
          color = "grey70",
          linewidth = 0.1) +
  geom_sf(data = wipad_sf,
          fill = "forestgreen",
          color = NA,
          alpha = 0.6) +
  theme_void() +
  labs(title = "Protected Areas Over Atlas Blocks")



# WIBBA-PAD subset overlap
wipad_sdnr <- wipad_sf %>%
  filter(Own_Name == "SDNR")

ggplot() +
  geom_sf(data = wipad_sdnr,
          fill = "darkgreen",
          color = NA,
          alpha = 0.7) +
  theme_void() +
  labs



wipad_sdnr <- wipad_sf %>%
  filter(Own_Name == "SDNR")

ggplot() +
  # Fill SDNR lands
  geom_sf(data = wipad_sdnr,
          fill = "forestgreen",
          color = NA,
          alpha = 0.8) +
  
  # Outline all PAD lands
  geom_sf(data = wipad_sf,
          fill = NA,
          color = "black",
          linewidth = 0.2) +
  
  theme_void() +
  labs(
    title = "Wisconsin Protected Areas",
    subtitle = "SDNR Lands Highlighted"
  )






# Histogram: Distribution of PA across blocks
pa_raw_df <- pa_raw_df %>%
  mutate(
    pa_bin = case_when(
      pa_percent == 0 ~ "0%",
      pa_percent > 0 & pa_percent <= 25 ~ ">0–25%",
      pa_percent > 25 & pa_percent <= 50 ~ ">25–50%",
      pa_percent > 50 & pa_percent <= 75 ~ ">50–75%",
      pa_percent > 75 ~ ">75–100%"
    )
  )

pa_raw_df$pa_bin <- factor(
  pa_raw_df$pa_bin,
  levels = c("0%", ">0–25%", ">25–50%", ">50–75%", ">75–100%")
)

# Text sizes
title_size  <- 16
axis_title_size <- 13
axis_text_size  <- 12

# PA counts per bin
pa_counts <- pa_raw_df %>%
  group_by(pa_bin) %>%
  summarise(count = n(), .groups = "drop")

# Plot with counts on top
ggplot(pa_counts, aes(x = pa_bin, y = count)) +
  geom_col(fill = "orange", color = "black") +
  geom_text(aes(label = count), 
            vjust = -0.5,  # moves text above bar
            size = 4,       # text size
            fontface = "bold") +
  labs(
    title = "Distribution of Protected Area Coverage Across Atlas Blocks",
    x = "Protected Area Coverage (% per block)",
    y = "Number of Atlas Blocks"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = title_size,
      margin = margin(b = 12)
    ),
    axis.title.x = element_text(
      face = "bold",
      size = axis_title_size,
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = axis_title_size,
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(size = axis_text_size, angle = 0, vjust = 0.5),
    axis.text.y = element_text(size = axis_text_size)
  ) +
  ylim(0, max(pa_counts$count) * 1.15)  # add space for labels above bars






### PREDICTED PROBABILITIES w PA ###

# PA z-scaled 
blocks_pa_z_sf <- blocks_all_sf %>%
  left_join(mod_data_all %>%
              dplyr::select(atlas_block, pa_prop),
            by = "atlas_block")

ggplot() +
  geom_sf(data = blocks_pa_z_sf,
          aes(fill = pa_prop),
          color = NA) +
  geom_sf(data = wipad_sf,
          fill = NA,
          color = "black",
          linewidth = 0.2) +
  scale_fill_viridis_c(
    option = "C",
    name = "Protected Area\n(Z-score)"
  ) +
  theme_minimal() +
  labs(
    title = "Z-Scaled Protected Area Coverage",
    subtitle = "Used in Occupancy Modeling"
  )


# Predicted response x pa
predict_by_block <- function(model_name) {
  
  model_obj <- global_glm_models[[model_name]]
  model_data <- data_dir[[model_name]]
  
  # Predict probabilities
  model_data$pred_prob <- predict(
    model_obj,
    type = "response"
  )
  
  # Aggregate to block level
  block_preds <- model_data %>%
    group_by(atlas_block) %>%
    summarise(pred_prob = mean(pred_prob, na.rm = TRUE),
              .groups = "drop")
  
  block_preds$model_name <- model_name
  
  return(block_preds)
}


block_predictions <- lapply(names(pa_models), # predictions 
                            predict_by_block)
block_predictions_df <- bind_rows(block_predictions)

blocks_pred_sf <- blocks_all_sf %>% # join predictions to geometry
  left_join(block_predictions_df,
            by = "atlas_block")


plot_prediction_map <- function(mod_name) {
  
  ggplot(blocks_pred_sf %>% 
           filter(model_name == mod_name)) +
    geom_sf(aes(fill = pred_prob),
            color = NA) +
    geom_sf(data = wipad_sf,
            fill = NA,
            color = "black",
            linewidth = 0.2) +
    scale_fill_viridis_c(
      option = "C",
      name = "Predicted\nProbability"
    ) +
    theme_void() +
    labs(
      title = paste("Predicted Probability -", mod_name)
    )
}

plot_prediction_map("RLL_col")
plot_prediction_map("RLL_ext")
plot_prediction_map("DNR_col")
plot_prediction_map("DNR_ext")


# w/o PA outlines
plot_prediction_map <- function(mod_name) {
  
  ggplot(blocks_pred_sf %>% 
           filter(model_name == mod_name)) +
    geom_sf(aes(fill = pred_prob),
            color = NA) +
    scale_fill_viridis_c(
      option = "C",
      name = "Predicted\nProbability"
    ) +
    theme_void() +
    labs(
      title = paste("Predicted Probability -", mod_name)
    )
}

plot_prediction_map("RLL_col")
plot_prediction_map("RLL_ext")
plot_prediction_map("DNR_col")
plot_prediction_map("DNR_ext")





















### SPECIES-RESPONSE MAPS ###

blocks_spp_sf <- blocks_all_sf %>%
  left_join(spp_blocks, by = "atlas_block")

blocks_rll_sf <- blocks_spp_sf %>%
  filter(atlas_block %in% blocks_rll)

# Count blocks per response type
response_counts <- blocks_rll_sf %>%
  st_set_geometry(NULL) %>%  # drop geometry for counting
  group_by(transition_state) %>%
  summarise(n_blocks = n(), .groups = "drop")

# Create labels with counts
response_labels <- paste0(response_counts$transition_state, " (", response_counts$n_blocks, ")")
names(response_labels) <- response_counts$transition_state

# Colors
state_colors <- c(
  "Colonization" = "darkorchid",
  "Persistence"  = "darkslategray3",
  "Extinction"   = "orange",
  "Absence"      = "#bdbdbd"
)


# Text size
title_size       <- 16
subtitle_size    <- 14
legend_title_size <- 12
legend_text_size  <- 11


# Plot
ggplot(blocks_rll_sf) +
  geom_sf(aes(fill = transition_state), color = NA, linewidth = 0.2) +
  scale_fill_manual(
    values = state_colors,
    labels = response_labels,
    na.value = "white",
    name = "Occurence State",
    guide = guide_legend(
      label.position = "left",   # moves the text to the left of the swatch
      keywidth = 1.5,            # adjust swatch width
      keyheight = 0.8,
      reverse = FALSE
    )
  ) +
  labs(
    title = paste0("Atlas Occurence of ", spp_name),
    subtitle = "Atlas 1 (1995-2000) → Atlas 2 (2015-2019)"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

