#############
### SETUP ###
#############

### --- LOAD PACKAGES --- ###

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  sf,
  dplyr,
  ggplot2,
  viridis,
  viridisLite,
  RColorBrewer,
  grid,
  cowplot)



### --- FLEXIBLE SPECS --- ###

# Species to map #
spp_name <- "Canada Warbler"



### --- GEOMETRY FILES --- ###

# Atlas Block Geometry #
# Native CRS: WI Transverse Mercator, EPSG:3071

blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_Breeding_Bird_Atlas_Blocks.shp") %>%
  rename(atlas_block = BLOCK_ID)
crs(blocks_all_sf)

blocks_rll_sf <- blocks_all_sf %>%
  filter(atlas_block %in% blocks_rll)
crs(blocks_rll_sf)


# transformation, CRS: Conus Albers, NAD83, EPSG:5070
blocks_all_sf <- st_transform(blocks_all_sf, 5070)
blocks_rll_sf <- st_transform(blocks_rll_sf, 5070)

st_crs(blocks_all_sf)
st_crs(blocks_rll_sf)




# WI State outline
state_outline_sf <- st_union(blocks_all_sf)
state_outline_sf <- st_sf(geometry = state_outline_sf)

# Plot to check
ggplot() +
  geom_sf(data = state_outline_sf,
          fill = NA,        # no fill
          color = "black",  # outline color
          linewidth = 0.5) + 
  theme_void() +
  labs(title = "Wisconsin State Outline")


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







### --- VISUALIZATIONS --- ###

# BASE MAPS #

# Atlas Blocks #
ggplot() +
  geom_sf(data = blocks_all_sf, # blocks_rll_sf, blocks_dnr_sf
          fill = "white",
          color = "grey40",
          linewidth = 0.2) +
  theme_void() +
  labs(title = "Wisconsin Breeding Bird Atlas Blocks")


# Protected Area (all) #
ggplot() +
  geom_sf(data = wipad_sf,
          color = "turquoise",
          fill = "turquoise", # boundaries, add arg linewidth = 
          alpha = 0.6) +
  theme_void() +
  labs(title = "Protected Areas in Wisconsin")



# Protected Area (ownership) #
owner_colors <- c(
  "SDNR" = "firebrick",
  "TNC"  = "cornflowerblue",
  "NAS"  = "plum",
  "USFS" = "lightgreen"
)

wipad_partners <- wipad_sf %>%
  filter(
    Own_Name %in% c("SDNR", "USFS") | 
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


# Plot
ggplot() +
  geom_sf(data = wipad_partners,
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
    subtitle = "Partner Ownerships Highlighted"
  )




# Protected Area (gap status) #
gap_colors <- c(
  "1" = "#006400", 
  "2" = "#41ab5d",
  "3" = "#fe9929",
  "4" = "#de2d26"
)

wipad_gap <- wipad_sf %>%
  filter(GAP_Sts %in% c(1, 2, 3, 4)) %>%
  mutate(
    GAP_Sts = factor(GAP_Sts, levels = c(1, 2, 3, 4))
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











### SPECIES-RESPONSE MAPS ###

# Join two species response dfs, mod_col_rll and mod_ext_rll


# Colonization / Absence
col_states <- mod_col_rll %>%
  dplyr::select(atlas_block, col) %>%
  mutate(transition_state = ifelse(col == 1, "Colonization", "Absence")) %>%
  dplyr::select(atlas_block, transition_state)

# Extinction / Persistence
ext_states <- mod_ext_rll %>%
  dplyr::select(atlas_block, ext) %>%
  mutate(transition_state = ifelse(ext == 1, "Extinction", "Persistence")) %>%
  dplyr::select(atlas_block, transition_state)

# Combine
blocks_species <- bind_rows(col_states, ext_states)


blocks_rll_sf <- blocks_rll_sf %>%
  left_join(blocks_species, by = "atlas_block")


response_counts <- blocks_rll_sf %>%
  st_set_geometry(NULL) %>%
  group_by(transition_state) %>%
  summarise(n_blocks = n(), .groups = "drop")


blocks_rll_sf$transition_state <- factor(
  blocks_rll_sf$transition_state,
  levels = c("Colonization", "Extinction", "Persistence", "Absence")
)


# Create labels with counts
response_labels <- paste0(response_counts$transition_state, " (", response_counts$n_blocks, ")")
names(response_labels) <- response_counts$transition_state

# Colors
state_colors <- c(
  "Colonization" = "darkorchid",
  "Persistence"  = "turquoise",
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
  geom_sf(aes(fill = transition_state), color = NA) +
  
  
  scale_fill_manual(
    values = state_colors,
    labels = response_labels,
    na.value = "white",
    name = paste0(spp_name, " Occurence State"),
    guide = guide_legend(
      label.position = "left",   # moves the text to the left of the swatch
      keywidth = 1,            # adjust swatch width
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



### Data Pool Maps
### --- SPECIES RESPONSE MAPS: POOLED MODEL BINS --- ###

# Colonization model universe = Colonization + Absence
col_bin <- mod_col_rll %>%
  dplyr::select(atlas_block) %>%
  mutate(model_bin = "Colonization / Absence")

# Extinction model universe = Extinction + Persistence
ext_bin <- mod_ext_rll %>%
  dplyr::select(atlas_block) %>%
  mutate(model_bin = "Extinction / Persistence")

# Combine pooled bins
blocks_species_bins <- bind_rows(col_bin, ext_bin)

# Join to sf
blocks_bins_sf <- blocks_rll_sf %>%
  dplyr::select(atlas_block, geometry) %>%
  left_join(blocks_species_bins, by = "atlas_block")

# Count blocks in each pooled bin
bin_counts <- blocks_bins_sf %>%
  st_set_geometry(NULL) %>%
  count(model_bin)

# Set factor order
blocks_bins_sf$model_bin <- factor(
  blocks_bins_sf$model_bin,
  levels = c(
    "Colonization / Absence",
    "Extinction / Persistence"
  )
)

# Legend labels with total counts
bin_labels <- paste0(bin_counts$model_bin, " (", bin_counts$n, ")")
names(bin_labels) <- bin_counts$model_bin

# Colors
bin_colors <- c(
  "Colonization / Absence"   = "darkorchid",
  "Extinction / Persistence" = "orange"
)

# Plot
ggplot(blocks_bins_sf) +
  geom_sf(aes(fill = model_bin), color = NA) +
  scale_fill_manual(
    values = bin_colors,
    labels = bin_labels,
    na.value = "white",
    name = "Occurence State Model Bins",
    guide = guide_legend(
      label.position = "left",
      keywidth = 1.5,
      keyheight = 0.8
    )
  ) +
  labs(
    title = paste0("Modeling State Space of ", spp_name),
    subtitle = "Colonization–Absence vs Extinction–Persistence blocks"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )







### --- MODEL RESPONSE MAPS --- ###
# Predicted Probabilities w/ PA

# ENV AND PA
predict_by_block <- function(model_name){
  
  model_obj  <- pa_glm_models[[model_name]]
  model_data <- data_dir[[model_name]]
  
  model_data$pred_prob <- predict(
    model_obj,
    type = "response"
  )
  
  block_preds <- model_data %>%
    group_by(atlas_block) %>%
    summarise(
      pred_prob = mean(pred_prob, na.rm = TRUE),
      .groups = "drop"
    )
  
  block_preds$model_name <- model_name
  
  block_preds
}


block_predictions <- lapply(
  c("RLL_col","RLL_ext"),
  predict_by_block
)

block_predictions_df <- bind_rows(block_predictions)


blocks_pred_sf <- blocks_rll_sf %>%
  left_join(block_predictions_df,
            by = "atlas_block")



plot_prediction_map <- function(mod_name){
  
  resp <- ifelse(grepl("col", mod_name),
                 "P(Colonization)",
                 "P(Extinction)")
  
  ggplot(
    blocks_pred_sf %>% 
      filter(model_name == mod_name)
  ) +
    geom_sf(aes(fill = pred_prob),
            color = NA) +
    scale_fill_viridis_c(
      option = "C",
      direction = 1,
      limits = c(0,1),
      name = resp
    ) +
    theme_void() +
    labs(
      title = paste("Predicted Probability -", mod_name)
    )
}


plot_prediction_map("RLL_col")
plot_prediction_map("RLL_ext")



### PA EFFECT PLOT

### --- FACETED PA OCCUPANCY EFFECT MAP --- ###


# ==============================
# 1) FUNCTION: partial PA effect by block
# ==============================
predict_pa_effect_by_block <- function(model_name){
  
  model_obj  <- pa_glm_models[[model_name]]
  model_data <- data_dir[[model_name]]
  
  # Prediction with observed PA
  pred_obs <- predict(
    model_obj,
    newdata = model_data,
    type = "response"
  )
  
  # Counterfactual: remove PA
  no_pa_data <- model_data
  no_pa_data$pa_prop <- 0
  
  pred_no_pa <- predict(
    model_obj,
    newdata = no_pa_data,
    type = "response"
  )
  
  # Partial PA contribution
  model_data$pa_effect <- pred_obs - pred_no_pa
  
  # Average to atlas block
  model_data %>%
    group_by(atlas_block) %>%
    summarise(
      pa_effect = mean(pa_effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(model_name = model_name)
}


# ==============================
# 2) GET COL + EXT BLOCK EFFECTS
# ==============================
col_pa <- predict_pa_effect_by_block("RLL_col") %>%
  mutate(
    model_bin = "Colonization / Absence",
    pa_effect_signed = pa_effect
  )

ext_pa <- predict_pa_effect_by_block("RLL_ext") %>%
  mutate(
    model_bin = "Extirpation / Persistence",
    pa_effect_signed = -pa_effect
  )


# ==============================
# 3) COMBINE FOR FACETING
# ==============================
pa_turnover_facet_df <- bind_rows(col_pa, ext_pa)


# ==============================
# 4) SHARED SYMMETRIC COLOR SCALE
# ==============================
max_abs <- max(
  abs(pa_turnover_facet_df$pa_effect_signed),
  na.rm = TRUE
)


# ==============================
# 5) JOIN TO SPATIAL BLOCKS
# ==============================
blocks_pa_facet_sf <- blocks_rll_sf %>%
  dplyr::select(atlas_block, geometry) %>%
  left_join(pa_turnover_facet_df, by = "atlas_block")


# ==============================
# 6) FINAL FACETED MAP
# ==============================

ggplot() +
  
  # block-level PA effect
  geom_sf(
    data = blocks_pa_facet_sf,
    aes(fill = pa_effect_signed),
    color = NA
  ) +
  
  # subtle PAD overlay
  geom_sf(
    data = wipad_sf,
    fill = "grey40",
    alpha = 0.05
  ) +
  
  facet_wrap(~model_bin) +
  
  scale_fill_gradient2(
    low = "orange",
    mid = "white",
    high = "darkorchid",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    oob = scales::squish,
    name = "PA Effect",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "cm"),
      barheight = unit(0.4, "cm")
    )
  ) +
  
  coord_sf(expand = FALSE) +
  theme_void() +
  
  labs(
    title = paste("Protected Area Contribution to", spp_name, "Occurrence"),
    subtitle = "Purple = supports occurrence, Orange = undermines occurrence"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 12
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.position = "bottom",   # puts legend below plot
    legend.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )










### --- MODEL RESPONSE MAPS: COMBINED COL/EXT --- ###

### Conditional Partial Effect Plot ###
ScalePA <- function(pa_raw, raw_df = covars_raw_rll) {
  mu <- mean(raw_df$pa_prop, na.rm = TRUE)
  sdv <- sd(raw_df$pa_prop, na.rm = TRUE)
  (pa_raw - mu) / sdv
}

MakePAResponseCurve <- function(model_col,
                                model_ext,
                                col_data,
                                ext_data,
                                raw_df,
                                pa_seq = seq(0, 1, length.out = 100)) {
  
  pa_scaled <- ScalePA(pa_seq, raw_df)
  
  all_data <- bind_rows(col_data, ext_data)
  
  vars_col <- all.vars(formula(model_col))[-1]
  vars_ext <- all.vars(formula(model_ext))[-1]
  all_vars <- unique(c(vars_col, vars_ext))
  
  baseline <- lapply(all_vars, function(v) {
    if (v == "pa_prop") return(0)
    median(all_data[[v]], na.rm = TRUE)
  })
  names(baseline) <- all_vars
  
  pred_data <- as.data.frame(baseline)
  pred_data <- pred_data[rep(1, length(pa_seq)), ]
  pred_data$pa_prop <- pa_scaled
  
  col_pred <- predict(model_col, newdata = pred_data,
                      type = "link", se.fit = TRUE)
  
  ext_pred <- predict(model_ext, newdata = pred_data,
                      type = "link", se.fit = TRUE)
  
  tibble(
    pa_raw = pa_seq,
    col_fit = plogis(col_pred$fit),
    col_lwr = plogis(col_pred$fit - 1.96 * col_pred$se.fit),
    col_upr = plogis(col_pred$fit + 1.96 * col_pred$se.fit),
    ext_fit = plogis(ext_pred$fit),
    ext_lwr = plogis(ext_pred$fit - 1.96 * ext_pred$se.fit),
    ext_upr = plogis(ext_pred$fit + 1.96 * ext_pred$se.fit)
  )
}


MakeObservedPA <- function(col_data,
                           ext_data,
                           raw_df,
                           bins = 12) {
  
  raw_lookup <- raw_df %>%
    dplyr::select(atlas_block, pa_raw = pa_prop) %>%
    distinct()
  
  col_obs <- col_data %>%
    left_join(raw_lookup, by = "atlas_block") %>%
    mutate(bin = cut(pa_raw, breaks = bins)) %>%
    group_by(bin) %>%
    summarise(
      pa_raw = mean(pa_raw, na.rm = TRUE),
      prob = mean(col, na.rm = TRUE),
      process = "Colonization",
      .groups = "drop"
    )
  
  ext_obs <- ext_data %>%
    left_join(raw_lookup, by = "atlas_block") %>%
    mutate(bin = cut(pa_raw, breaks = bins)) %>%
    group_by(bin) %>%
    summarise(
      pa_raw = mean(pa_raw, na.rm = TRUE),
      prob = mean(ext, na.rm = TRUE),
      process = "Extinction",
      .groups = "drop"
    )
  
  bind_rows(col_obs, ext_obs)
}

curve_pa <- MakePAResponseCurve(
  pa_glm_models$RLL_col,
  pa_glm_models$RLL_ext,
  mod_col_rll,
  mod_ext_rll,
  covars_raw_rll
)

obs_pa <- MakeObservedPA(
  mod_col_rll,
  mod_ext_rll,
  covars_raw_rll
)


ggplot(curve_pa, aes(x = pa_raw)) +
  geom_ribbon(aes(ymin = col_lwr, ymax = col_upr, fill = "Colonization"),
              alpha = 0.2) +
  geom_line(aes(y = col_fit, color = "Colonization"), linewidth = 1) +
  
  geom_ribbon(aes(ymin = ext_lwr, ymax = ext_upr, fill = "Extinction"),
              alpha = 0.2) +
  geom_line(aes(y = ext_fit, color = "Extinction"), linewidth = 1) +
  
  geom_point(data = obs_pa,
             aes(y = prob, color = process),
             size = 2.5) +
  
  scale_color_manual(values = c(
    Colonization = "darkorchid",
    Extinction = "orange"
  )) +
  scale_fill_manual(values = c(
    Colonization = "darkorchid",
    Extinction = "orange"
  )) +
  scale_x_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Protected Area (%)",
    y = "Transition Probability",
    title = paste(spp_name, "Conditional Partial PA Response")
  ) +
  theme_minimal()





### Avg Marginal Predictions Plot ###
ScalePA <- function(pa_raw, raw_df = covars_raw_rll) {
  mu <- mean(raw_df$pa_prop, na.rm = TRUE)
  sdv <- sd(raw_df$pa_prop, na.rm = TRUE)
  (pa_raw - mu) / sdv
}


MakeAvgMarginalPA_CI <- function(model_col,
                                 model_ext,
                                 col_data,
                                 ext_data,
                                 raw_df,
                                 pa_seq = seq(0, 1, length.out = 100),
                                 n_sims = 500) {
  
  col_sims <- arm::sim(model_col, n.sims = n_sims)
  ext_sims <- arm::sim(model_ext, n.sims = n_sims)
  
  curve_list <- lapply(pa_seq, function(pa_val) {
    
    pa_scaled <- ScalePA(pa_val, raw_df)
    
    col_new <- col_data
    ext_new <- ext_data
    
    col_new$pa_prop <- pa_scaled
    ext_new$pa_prop <- pa_scaled
    
    # model matrix
    X_col <- model.matrix(formula(model_col), col_new)
    X_ext <- model.matrix(formula(model_ext), ext_new)
    
    # simulate predictions
    col_draws <- apply(col_sims@coef, 1, function(b) {
      mean(plogis(X_col %*% b))
    })
    
    ext_draws <- apply(ext_sims@coef, 1, function(b) {
      mean(plogis(X_ext %*% b))
    })
    
    tibble(
      pa_raw = pa_val,
      col_fit = mean(col_draws),
      col_lwr = quantile(col_draws, 0.025),
      col_upr = quantile(col_draws, 0.975),
      ext_fit = mean(ext_draws),
      ext_lwr = quantile(ext_draws, 0.025),
      ext_upr = quantile(ext_draws, 0.975)
    )
  })
  
  bind_rows(curve_list)
}


curve_avg_pa <- MakeAvgMarginalPA_CI(
  pa_glm_models$RLL_col,
  pa_glm_models$RLL_ext,
  mod_col_rll,
  mod_ext_rll,
  covars_raw_rll
)


obs_pa <- MakeObservedPA(
  mod_col_rll,
  mod_ext_rll,
  covars_raw_rll,
  bins = 12
)

ggplot(curve_avg_pa, aes(x = pa_raw)) +
  
  geom_ribbon(
    aes(ymin = col_lwr, ymax = col_upr, fill = "Colonization"),
    alpha = 0.2
  ) +
  geom_line(
    aes(y = col_fit, color = "Colonization"),
    linewidth = 1.2
  ) +
  
  geom_ribbon(
    aes(ymin = ext_lwr, ymax = ext_upr, fill = "Extinction"),
    alpha = 0.2
  ) +
  geom_line(
    aes(y = ext_fit, color = "Extinction"),
    linewidth = 1.2
  ) +
  
  geom_point(
    data = obs_pa,
    aes(y = prob, color = process),
    size = 2.5,
    alpha = 0.85
  ) +
  
  scale_color_manual(values = c(
    Colonization = "darkorchid",
    Extinction = "orange"
  )) +
  scale_fill_manual(values = c(
    Colonization = "darkorchid",
    Extinction = "orange"
  )) +
  
  scale_x_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  
  labs(
    x = "Protected Area (%)",
    y = "Average Predicted Transition Probability",
    title = paste(spp_name, "Average Marginal PA Response"),
    color = "Process",
    fill = "Process"
  ) +
  theme_minimal()






### --- MODEL RESPONSE MAPS: COMBINED COL/EXT W PA --- ###


# --- 1. Pull colonization & extinction predictions from block_predictions_df --- #
col_preds <- block_predictions_df %>%
  filter(model_name == "RLL_col") %>%
  dplyr::select(atlas_block, col_prob = pred_prob)

ext_preds <- block_predictions_df %>%
  filter(model_name == "RLL_ext") %>%
  dplyr::select(atlas_block, ext_prob = pred_prob)

# --- 2. Build full block table --- #
full_blocks <- blocks_rll_sf %>%
  st_set_geometry(NULL) %>%
  dplyr::select(atlas_block) %>%
  left_join(col_preds, by = "atlas_block") %>%
  left_join(ext_preds, by = "atlas_block")

# --- 3. Z-score standardization --- #
full_blocks <- full_blocks %>%
  mutate(
    col_scaled = ifelse(!is.na(col_prob),
                        (col_prob - mean(col_prob, na.rm = TRUE)) / sd(col_prob, na.rm = TRUE),
                        NA),
    ext_scaled = ifelse(!is.na(ext_prob),
                        (ext_prob - mean(ext_prob, na.rm = TRUE)) / sd(ext_prob, na.rm = TRUE),
                        NA),
    # Combined response surface: positive = colonization, negative = extinction
    response_surface = case_when(
      !is.na(col_scaled) & !is.na(ext_scaled) ~ col_scaled - ext_scaled,
      !is.na(col_scaled) & is.na(ext_scaled)  ~ col_scaled,
      is.na(col_scaled) & !is.na(ext_scaled)  ~ -ext_scaled,
      TRUE ~ NA_real_
    )
  )

# --- 4. Rescale response_surface to -1 -> 1 --- #
max_abs <- max(abs(full_blocks$response_surface), na.rm = TRUE)
full_blocks <- full_blocks %>%
  mutate(response_scaled = response_surface / max_abs)

# --- 5. Join back to spatial blocks --- #
blocks_combined_sf <- blocks_rll_sf %>%
  left_join(full_blocks %>% dplyr::select(atlas_block, response_scaled), by = "atlas_block")

# --- 6. Define legend ticks --- #
numeric_ticks <- c(-1, -0.5, 0, 0.5, 1)

# --- 7. Map with PA overlay --- #
map_pa_plot <- ggplot() +
  geom_sf(data = blocks_combined_sf,
          aes(fill = response_scaled),
          color = NA) +
  geom_sf(data = wipad_sf,
          fill = "grey40",
          alpha = 0.2,
          color = "black",
          linewidth = 0.03) +
  scale_fill_gradient2(
    low = "orange",
    mid = "white",
    high = "darkorchid",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90",
    name = "Pressure",
    guide = "none"
  ) +
  coord_sf(expand = FALSE) +
  theme_void() +
  labs(title = paste("Colonization vs Extinction Pressure -", spp_name)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# --- 8. Legend --- #
legend_plot <- ggplot(blocks_combined_sf) +
  geom_tile(aes(x = response_scaled, y = 1, fill = response_scaled)) +
  scale_fill_gradient2(
    low = "orange",
    mid = "white",
    high = "darkorchid",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = numeric_ticks,
    labels = numeric_ticks,
    guide = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(8, "cm"),
      barheight = unit(0.35, "cm"),
      ticks.colour = "black",
      title.position = "top",
      label.position = "bottom"
    ),
    name = NULL
  ) +
  annotate(
    "text", x = c(-1, 0, 1), y = 1.05,
    label = c("Extinction", "Stable", "Colonization"),
    size = 3.5,
    fontface = "bold"
  ) +
  theme_void() +
  theme(legend.position = "bottom")

# --- 9. Combine map + legend --- #
plot_grid(
  map_pa_plot,
  legend_plot,
  ncol = 1,
  rel_heights = c(0.9, 0.1)
)





### --- SPECIES PA CATERPILLAR PLOT --- ### 


write.csv(rll_results_df,"outputs/data/basepa_effects_rll.csv")
rll_results_df <- read.csv("outputs/data/basepa_effects_rll.csv")


# Caterpillar Plot

caterpillar_df <- rll_results_df %>%
  group_by(species, response) %>%
  mutate(strength = max(abs(estimate), na.rm = TRUE)) %>%
  ungroup()


# Plot

PlotEffects <- function(df, species_keep = NULL) {
  
  # Species filter
  if (!is.null(species_keep)) {
    df <- df %>%
      filter(species %in% species_keep)
  }
  
  # Rank order
  df <- df %>%
    group_by(species) %>%
    mutate(strength = max(abs(estimate), na.rm = TRUE)) %>%
    ungroup()
  
  dodge <- position_dodge(width = 0.5)
  
  ggplot(
    df,
    aes(
      x = estimate,
      y = reorder(species, strength),
      color = response,
      group = response
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      width = 0.2,
      position = dodge
    ) +
    
    geom_point(
      aes(shape = sig),
      size = 2.8,
      stroke = 1,
      fill = "white",
      position = dodge
    ) +
    
    scale_shape_manual(values = c("yes" = 16, "no" = 1)) +
    
    scale_color_manual(
      values = c("col" = "#1f78b4", "ext" = "tomato")
    ) +
    
    labs(
      x = "PA Effect (log-odds)",
      y = "Species",
      color = "Response",
      shape = "Significant"
    ) +
    
    theme_minimal()
}




species_subset <- c(
  "Canada Warbler",
  "Hooded Warbler",
  "Cerulean Warbler",
  "Winter Wren",
  "Olive-sided Flycatcher",
  "Henslow's Sparrow"
)


species_subset <- c(
  "Dickcissel",
  "Northern Cardinal",
  "Orchard Oriole",
  "Red-bellied Woodpecker",
  "Tufted Titmouse"
)


species_subset <- c(
  "Dickcissel",
  "Bobolink",
  "Eastern Meadowlark",
  "Grasshopper Sparrow",
  "Henslow's Sparrow",
  "Upland Sandpiper"
)


species_subset <- c(
  "Ruby-crowned Kinglet",
  "Canada Jay",
  "Olive-sided Flycatcher",
  "Canada Warbler",
  "Winter Wren",
  "Cape May Warbler",
  "Evening Grosbeak"
)



plot_rll <- PlotEffects(
  caterpillar_df,
  species_keep = species_subset
)

plot_rll







### BOLDED NAMES COMP CATERPILLAR 

# -----------------------------
# Define the two species groups
# -----------------------------

# Group you want emphasized (bold labels)
species_focus <- c(
  "Dickcissel",
  "Northern Cardinal",
  "Orchard Oriole",
  "Red-bellied Woodpecker"
)

# Comparison group
species_compare <- c(
  "Ruby-crowned Kinglet",
  "Olive-sided Flycatcher",
  "Canada Warbler",
  "Winter Wren"
)

# Combine both groups for plotting
species_subset <- c(species_focus, species_compare)


# -----------------------------
# Plot function
# -----------------------------

PlotEffects <- function(df,
                        species_keep = NULL,
                        species_bold = NULL) {
  
  # Filter to only chosen species
  if (!is.null(species_keep)) {
    df <- df %>%
      filter(species %in% species_keep)
  }
  
  # Rank species by strongest absolute effect
  df <- df %>%
    group_by(species) %>%
    mutate(strength = max(abs(estimate), na.rm = TRUE)) %>%
    ungroup()
  
  # Create axis labels (bold only selected species)
  species_order <- df %>%
    distinct(species, strength) %>%
    arrange(strength) %>%
    pull(species)
  
  axis_labels <- species_order
  
  if (!is.null(species_bold)) {
    axis_labels <- ifelse(
      species_order %in% species_bold,
      paste0("bold('", species_order, "')"),
      paste0("'", species_order, "'")
    )
  }
  
  names(axis_labels) <- species_order
  
  dodge <- position_dodge(width = 0.5)
  
  ggplot(
    df,
    aes(
      x = estimate,
      y = reorder(species, strength),
      color = response,
      group = response
    )
  ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    
    geom_errorbar(
      aes(
        xmin = conf.low,
        xmax = conf.high
      ),
      orientation = "y",
      width = 0.2,
      position = dodge
    ) +
    
    geom_point(
      aes(shape = sig),
      size = 2.8,
      stroke = 1,
      fill = "white",
      position = dodge
    ) +
    
    scale_shape_manual(
      values = c(
        "yes" = 16,
        "no" = 1
      )
    ) +
    
    scale_color_manual(
      values = c(
        "col" = "#1f78b4",
        "ext" = "tomato"
      )
    ) +
    
    scale_y_discrete(
      labels = parse(text = axis_labels)
    ) +
    
    labs(
      x = "Protected Area Effect (log-odds)",
      color = "Response",
      shape = "Significant"
    ) +
    
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 11)
    )
}


# -----------------------------
# Final plot
# -----------------------------

plot_rll <- PlotEffects(
  caterpillar_df,
  species_keep = species_subset,
  species_bold = species_focus   # <- bold only Dickcissel group
)

plot_rll
