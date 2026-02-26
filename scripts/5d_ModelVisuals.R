library(sf)
library(dplyr)
library(ggplot2)
library(RColorBrewer)



# Atlas blocks geometry
blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>% # load full block map
  rename(atlas_block = BLOCK_ID)
crs(blocks_all_sf)

blocks_background <- blocks_all_sf %>%
  filter(atlas_block %in% union(blocks_rll, blocks_dnr))


# PAD geometry
wipad_sf <- st_read("data/maps/wipad/pad_wi_polygons.shp")
crs(wipad_sf)
wipad_sf <- st_transform(wipad_sf, 3071)


# Join RAW PA values and atlas blocks
pa_raw_df <- covars_raw_rll %>%
  dplyr::select(atlas_block, pa_percent)

blocks_pa_sf <- blocks_all_sf %>%
  left_join(pa_raw_df, by = "atlas_block")

summary(blocks_pa_sf$pa_percent)




# Visualizations

## PROTECTED AREA ##
# Raw PA Coverage by Atlas Block
ggplot() +
  geom_sf(data = blocks_pa_sf,
          aes(fill = pa_percent),
          color = NA) +
  geom_sf(data = wipad_sf,
          fill = NA,
          color = "black",
          linewidth = 0.2) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    name = "Protected Area (%)"
  ) +
  theme_minimal() +
  labs(
    title = "Protected Area Coverage by Atlas Block",
    subtitle = "Raw Percent Coverage",
    caption = "Values represent raw % of block under protection"
  )



# Binned RAW PA Coverage by Atlas Block
breaks <- quantile(blocks_pa_sf$pa_percent,
                   probs = seq(0, 1, 0.2),
                   na.rm = TRUE)

breaks <- unique(breaks)  # prevents non-unique break error

labels <- paste0(
  round(breaks[-length(breaks)], 1),
  "–",
  round(breaks[-1], 1),
  "%"
)

blocks_pa_sf <- blocks_pa_sf %>%
  mutate(pa_bin = cut(pa_percent,
                      breaks = breaks,
                      include.lowest = TRUE,
                      labels = labels)) %>%
  mutate(pa_bin = factor(pa_bin,
                         levels = rev(labels)))

ggplot() +
  geom_sf(data = blocks_pa_sf,
          aes(fill = pa_bin),
          color = NA) +
  geom_sf(data = wipad_sf,
          fill = NA,
          color = "black",
          linewidth = 0.2) +
  scale_fill_brewer(
    palette = "YlGnBu",
    name = "Protected Area (%)\nQuintiles"
  ) +
  theme_minimal() +
  labs(
    title = "Protected Area Coverage by Atlas Block",
    subtitle = "Raw Percent (Quintiles)"
  )




### PREDICTED PROBABILITIES w PA ###

# PA z-scaled 
blocks_pa_z_sf <- blocks_all_sf %>%
  left_join(mod_data_all %>%
              dplyr::select(atlas_block, pa_percent),
            by = "atlas_block")

ggplot() +
  geom_sf(data = blocks_pa_z_sf,
          aes(fill = pa_percent),
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


block_predictions <- lapply(names(global_glm_models), # predictions 
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
    theme_minimal() +
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
    theme_minimal() +
    labs(
      title = paste("Predicted Probability -", mod_name)
    )
}

plot_prediction_map("RLL_col")
plot_prediction_map("RLL_ext")
plot_prediction_map("DNR_col")
plot_prediction_map("DNR_ext")