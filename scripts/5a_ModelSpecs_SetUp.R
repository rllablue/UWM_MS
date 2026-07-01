### Script: 5a_ModelSpecs_SetUp.R
### Purpose:


#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  
  # Core
  dplyr,
  tidyverse,
  magrittr,
  purrr,
  stringr,
  readxl,
  
  # Diagnostics
  nnet,
  stats,
  broom,
  corrplot,
  car,
  performance,
  DescTools,
  Metrics,
  pROC,
  rsample,
  lme4,
  pscl,
  AICcmodavg,
  MuMIn,
  arm,
  ncf,
  DHARMa,
  psych, 
  usdm,
  glmmTMB,
  
  # Visualization
  ggplot2,
  viridis,
  gt,
  sf,
  gpkg,
  webshot2
  
)

options(na.action = "na.fail") # required for MuMIn



### --- DATA --- ###

### --- Bird Data

#spp_zf_rll <- read.csv("data/summaries/spp_zf_rll.csv") %>% # df
#  dplyr::select(-X)
#spp_list <- unique(spp_zf_rll$common_name) # vector

#spp_guilds <- read.csv("data/summaries/species_guilds.csv")

#spp_zf_rll <- spp_zf_rll %>%
#  left_join(spp_guilds, by = "common_name")



### --- Atlas Block Data

# Native: EPSG:3071, Wisconsin Transverse Mercator
blocks_all_sf <- st_read("data/maps/wibba/Wisconsin_BBA_blocks.gpkg") %>%
  mutate(priority_block = factor(priority_block, levels = c("nonpriority", "priority")))
blocks_all <- blocks_all_sf$atlas_block # vector

blocks_rll_sf <- st_read("data/maps/wibba/Wisconsin_BBA_blocks_rll.gpkg") %>%
  mutate(priority_block = factor(priority_block, levels = c("nonpriority", "priority")))
blocks_rll <- blocks_rll_sf$atlas_block # vector


#st_write(
#  blocks_all_sf,
#  "data/maps/wibba/Wisconsin_BBA_blocks.gpkg",
#  delete_layer = TRUE
#)

#st_write(
#  blocks_rll_sf,
#  "data/maps/wibba/Wisconsin_BBA_blocks_rll.gpkg",
#  delete_layer = TRUE
#)



### --- Covariate Data

covs_raw_rll_df <- read.csv("data/summaries/covs_raw_rll_df.csv") # df

#covs_raw_rll_df <- covs_raw_rll_df %>%
#write.csv(covs_raw_rll_df, "data/summaries/covs_raw_rll_df.csv", row.names = FALSE)



### --- Consolidated Modeling Data 

#moddata_raw_rll_df <- spp_zf_rll %>% 
#  left_join(covs_raw_rll_df, 
#    by = "atlas_block") %>%
#  left_join(blocks_rll_sf %>%
#    st_drop_geometry() %>%
#    dplyr::select(atlas_block, priority_block),
#    by = "atlas_block"
#  )

#write.csv(moddata_raw_rll_df, "outputs/data/moddata_raw_rll_df.csv", row.names = FALSE)



#####################
### MODELING DATA ###
#####################

# FULL DATASET #
moddata_raw_rll_df <- read.csv("outputs/data/moddata_raw_rll_df.csv") %>%
  mutate(priority_block = factor(priority_block, levels = c("nonpriority", "priority")))


## SPECIES ##
spp_alpha <- "CAWA"
spp_name <- "Canada Warbler"



## SPECIES COUNTS ##
rll_spp_counts <- moddata_raw_rll_df %>%
  filter(atlas_block %in% blocks_rll) %>%
  count(common_name, transition_state) %>%
  pivot_wider(
    names_from  = transition_state,
    values_from = n,
    values_fill = 0
  )

# Counts for spp with sufficient sample size 
# flexible sample size thresholds
min_total <- 30
min_events <- 10

rll_valid_counts <- moddata_raw_rll_df %>%
  filter(atlas_block %in% blocks_rll) %>%
  count(common_name, transition_state) %>%
  pivot_wider(
    names_from   = transition_state,
    values_from  = n,
    values_fill  = 0
  ) %>%
  mutate(
    
    # modeling totals
    col_abs_total = Absence + Colonization,
    ext_per_total = Persistence + Extirpation,
    
    # event structure
    col_events     = Colonization,
    col_nonevents  = Absence,
    ext_events     = Extirpation,
    ext_nonevents  = Persistence,
    
    # validity thresholds
    col_valid = (
      col_abs_total   >= min_total &
        col_events    >= min_events &
        col_nonevents >= min_events
    ),
    ext_valid = (
      ext_per_total   >= min_total &
        ext_events    >= min_events &
        ext_nonevents >= min_events
    )
  ) %>%
  filter(col_valid | ext_valid) %>%
  dplyr::select(
    common_name,
    Absence,
    Colonization,
    Extirpation,
    Persistence,
    col_abs_total,
    ext_per_total,
    col_valid,
    ext_valid
  )




# --- SPECIES-RESPONSE MOD DATA SUBSETTING --- #

### Separate data into bins where A2 detection is either 1 or 0; not necessarily
# in precise probabilities between col, per, abs, ext, but more what promotes
# 'new' v. 'continued' colonization, ie. what promotes det = 1, as opposed to 
# what promotes det = 0. Colonization = col + abs data, Extirpation = ext + pre data
### Scaling of covars w/in subsets for relevant normalized values


# Helper: Build filtered modeling dfs
BuildSppRespDfs <- function(data,
                            species,
                            block_vector = NULL,
                            state_pairs,
                            response_name,
                            response_one) {
  
  df <- data %>%
    filter(alpha_code == species)
  
  if (!is.null(block_vector)) {
    df <- df %>% filter(atlas_block %in% block_vector)
  }
  
  df %>%
    filter(transition_state %in% state_pairs) %>%
    mutate(
      !!response_name := as.integer(
        ifelse(transition_state == response_one, 1, 0)
      )
    )
}


# Apply: Colonization Model (Col + Abs)
modcol_rll_raw_df <- BuildSppRespDfs(
  data = moddata_raw_rll_df,
  species = spp_alpha,
  block_vector = blocks_rll,
  state_pairs = c("Colonization", "Absence"),
  response_name = "col",
  response_one = "Colonization"
)


# Apply: Extirpation Model (Ext + Per)
modext_rll_raw_df <- BuildSppRespDfs(
  data = moddata_raw_rll_df,
  species = spp_alpha,
  block_vector = blocks_rll,
  state_pairs = c("Extirpation", "Persistence"),
  response_name = "ext",
  response_one = "Extirpation"
)


# Helper: Scale covariates within species-response subsets
ZScaleNumericData <- function(df,
                              exclude = c("col", "ext", "det_Atlas1", "det_Atlas2")) {
  
  nums <- names(df)[sapply(df, is.numeric)] 
  nums <- setdiff(nums, exclude)
  df[nums] <- lapply(df[nums], \(x) as.numeric(scale(x)))
  df
}

# Apply: independently z-scaled modeling data subsets
modcol_rll_z_df <- modcol_rll_raw_df %>% ZScaleNumericData()
modext_rll_z_df <- modext_rll_raw_df %>% ZScaleNumericData()


# Store in list
mod_dfs_all <- list(
  colmod = modcol_rll_z_df,
  extmod = modext_rll_z_df
)





##################
### COVARIATES ###
##################

# FULL DATASET #
chr_covs_all <- c("atlas_block", 
                  "common_name", "alpha_code", "guild", 
                  "transition_state")

stable_covs_all <- c("lon", "lat", 
                     "sr_diff", "priority_block", # effort control terms
                     "sac_kernel", # sac term
                     
                     "pa_prop", # total prop pa per block 
                     "sdnr_prop", "usfs_prop", "tnc_prop", "nas_prop", # manager prop of total pa per block
                     "gap1_prop", "gap2_prop", "gap3_prop") # gap status prop of total pa per block

land_covs_all <- c("water_open_base", "barren_land_base", "shrub_scrub_base", # base year values
                   "developed_open_base", "developed_low_base", "developed_med_base", "developed_high_base", 
                   "developed_lower_base", "developed_upper_base", "developed_total_base", 
                   "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base", "forest_total_base", 
                   "pasture_base", "cropland_base", "grassland_base", "grass_pasture_base",
                   "wetlands_woody_base", "wetlands_herb_base", "wetlands_total_base",
                   
                   "water_open_diff", "barren_land_diff", "shrub_scrub_diff", # change values
                   "developed_open_diff", "developed_low_diff", "developed_med_diff", "developed_high_diff", 
                   "developed_lower_diff", "developed_upper_diff", "developed_total_diff", 
                   "forest_deciduous_diff", "forest_evergreen_diff", "forest_mixed_diff", "forest_total_diff", 
                   "pasture_diff", "cropland_diff", "grassland_diff", "grass_pasture_diff",
                   "wetlands_woody_diff", "wetlands_herb_diff", "wetlands_total_diff")

climate_covs_all <- c("tmax_38yr", "tmin_38yr", "prcp_38yr", # base year values
                      
                      "tmax_diff", "tmin_diff", "prcp_diff") # change values



### --- COVARIATE SUBSETS --- ###

## ECOLOGICAL FILTERS ##
### Habitat guilds to bin land covariates

# Inputs
guildcovs_map <- list(
  
  forest = c("developed_total_base", 
             "forest_deciduous_base", "forest_mixed_base", "forest_evergreen_base", 
             "wetlands_woody_base", "wetlands_herb_base",
             
             "forest_total_diff", "wetlands_total_diff"),
  
  grass = c("developed_total_base","forest_total_base", "cropland_base",
            "grassland_base", "pasture_base",
            
            "grassland_diff", "pasture_diff"),
  
  marsh = c("water_open_base", "grass_pasture_base",
            "developed_total_base", "forest_total_base", 
            "wetlands_woody_base", "wetlands_herb_base",
            
            "wetlands_total_diff"), 
  
  water = c("water_open_base", "developed_total_base", "grass_pasture_base", 
            "forest_deciduous_base", "forest_mixed_base", "forest_evergreen_base", 
            "wetlands_woody_base", "wetlands_herb_base",
            
            "water_open_diff", "forest_total_diff", 
            "grass_pasture_diff", "wetlands_total_diff"),
  
  
  general = c("developed_total_base",
              "forest_total_base", "grass_pasture_base", "wetlands_total_base",
              
              "forest_total_diff", "grass_pasture_diff", "wetlands_total_diff") # "water_open_base"
  
)


# lookup table
sppguild_lookup <- moddata_raw_rll_df %>%
  distinct(alpha_code, guild)


# Helper: pull spp guild-covariate combinations
PullGuildCovs <- function(species, guild_lookup, guild_map) {
  
  guild_value <- guild_lookup %>%
    filter(alpha_code == species) %>%
    pull(guild)
  
  if (length(guild_value) == 0) {
    stop("Species lacks guild assignment.")
  }
  
  if (!guild_value %in% names(guild_map)) {
    stop("Guild not found in key.")
  }
  
  guild_map[[guild_value]]
}


# Apply: Species-Specific Covariate Sets
chr_covs_reduced <- c("atlas_block", "transition_state")
stable_covs_reduced <- c("sr_diff", "priority_block", "sac_kernel") # no PA in this script
land_covs_reduced <- PullGuildCovs(spp_alpha, sppguild_lookup, guildcovs_map)
climate_covs_reduced <- c("tmax_38yr", "prcp_38yr", # base year values
                          "tmax_diff", "tmin_diff", "prcp_diff") # change values

numeric_covs_reduced <- c(stable_covs_reduced, land_covs_reduced, climate_covs_reduced)



### STATISTICAL FILTERS ### 

# CORRELATIONS #
### Pair-wise, linear relationships (i.e. bivariate)
stable_covs_reduced <- c("sr_diff") # no priority_block, corrs cannot handle factors
numeric_covs_reduced <- c(stable_covs_reduced, land_covs_reduced, climate_covs_reduced)

covnames <- unique(numeric_covs_reduced)


VisualizeCorrelations <- function(data, covs, label = "") {
  
  covnames <- covs[covs %in% names(data)]
  
  # Pairwise Plots
  psych::pairs.panels(
    data[, covnames],  
    method = "pearson", 
    hist.col = "#00AFBB",
    density = FALSE,
    ellipses = FALSE,
    lm = TRUE,
    main = paste("Visualizing Predictor Relationships:", label)
  )
  
  # Histograms 
  print(
    data %>%
      dplyr::select(all_of(covnames)) %>%
      pivot_longer(
        cols = everything(), 
        names_to = "Variable", 
        values_to = "Value"
      ) %>%
      ggplot(aes(x = Value)) +
      geom_histogram(bins = 20, fill = "steelblue", color = "white") +
      facet_wrap(~Variable, scales = "free") +
      theme_minimal() +
      labs(
        title = paste("Predictor Distributions:", label),
        y = "Count", 
        x = "Value"
      )
  )
  
  # Correlation Matrix
  cor_mat <- cor(data[, covnames], method = 'spearman')
  
  corrplot::corrplot.mixed(
    cor_mat, 
    tl.pos = 'lt', 
    tl.cex = 0.8, 
    number.cex = 0.7, 
    addCoefasPercent = TRUE,
    title = paste("Spearman Correlation Matrix:", label),
    mar = c(0,0,2,0)
  )
  
  invisible(cor_mat)
}

# Apply: aseess pairwise corrs
corr_results <- lapply(
  names(mod_dfs_all),
  function(nm) {
    VisualizeCorrelations(
      data = mod_dfs_all[[nm]],
      covs = numeric_covs_reduced,
      label = nm
    )
  }
)

names(corr_results) <- names(mod_dfs_all)


# Helper: numerically identify correlated predictors
GetHighCorrs <- function(data, covs, threshold = 0.7) {
  
  covs <- covs[sapply(data[, covs, drop = FALSE], function(x)
    sd(x, na.rm = TRUE) > 0)]
  
  pc <- cor(data[, covs], use = "pairwise.complete.obs")
  
  # only upper triangle, excluding diagonal
  idx <- which(abs(pc) > threshold & upper.tri(pc), arr.ind = TRUE)
  
  if (nrow(idx) == 0) {
    message("No corrs above threshold.")
  } else {
    apply(idx, 1, function(i)
      cat(rownames(pc)[i[1]], "-", colnames(pc)[i[2]],
          "r =", pc[i[1], i[2]], "\n")
    )
  }
  
  invisible(pc)
}

corrs1 <- GetHighCorrs(modcol_rll_z_df, numeric_covs_reduced)
corrs2 <- GetHighCorrs(modext_rll_z_df, numeric_covs_reduced)



# Output covariates
#stable_covars_reduced

#guild_key
#land_covs_reduced <- c()

#climate_covs_reduced # "tmax_38yr", "tmax_diff", "prcp_38yr", "tmin_diff", "prcp_diff"
#climate_covs_reduced <- c()

#numeric_covs_reduced <- c(land_covs_reduced, climate_covs_reduced, stable_covs_reduced)



# VARIANCE INFLATION FACTOR #
### Single variable relationship to all others as a group (i.e. multivariate); 
# multicollinearity redundancy measure, i.e. too similar to uniquely est coefficients (ideally VIF < 5)
stable_covs_reduced <- c("sr_diff", "priority_block") # VIF can handle factors
numeric_covs_reduced <- c(stable_covs_reduced, land_covs_reduced, climate_covs_reduced)


responses <- c(
  colmod = "col",
  extmod = "ext"
)


GetVIFs <- function(data, response, covs, family = "binomial") {
  
  formula_str <- as.formula(
    paste(response, "~ .")
  )
  
  df <- data[, c(response, covs)]
  
  mod <- glm(formula_str,
             data = df, 
             family = family)
  
  vifs <- car::vif(mod)
  
  vifs_sorted <- sort(vifs, decreasing = TRUE)
  cat("\nVIF (sorted):\n")
  print(vifs_sorted)
  
  cat("\nAliased coefficients:\n")
  print(alias(mod))
  
  invisible(list(model = mod,
                 vif = vifs_sorted))
}


vif_results <- lapply(
  names(mod_dfs_all),
  function(nm) {
    
    cat("VIF for", nm, "\n")
    
    GetVIFs(
      data = mod_dfs_all[[nm]],
      response = responses[nm],
      covs = numeric_covs_reduced
    )
  }
)

names(vif_results) <- names(mod_dfs_all)


# covariate revisions
#land_covs_reduced <- c("developed_total_base", "grass_pasture_base",
#"forest_mixed_base", "forest_evergreen_base", 
#"wetlands_woody_base", "wetlands_herb_base",

#"forest_total_diff", "wetlands_total_diff")

#climate_covs_reduced # "tmax_38yr", "tmax_diff", "prcp_38yr", "tmin_diff", "prcp_diff"
#climate_covs_reduced <- c()

#numeric_covs_reduced <- c(land_covs_reduced, climate_covs_reduced, stable_covs_reduced)



#####################################
### MODEL CONSTRUCTION, SELECTION ### 
#####################################

### --- PROCESS, SET UP --- ###

### Steps ###

# 1) Partitioned AICc Model Selection
# Construct, fit, rank separate models partitioning land and climate covariates;
# carry over best covariates for each set into full ENV models.

# 2) Additive ENV Model Selection
# Construct, fit, rank top models AICc < 2 w/ additive land, climate, effort terms;
# extract reference model for each response subset.

# 3) Model Diagnostics
# Assess top candidate models for each response for correlations, dispersion, fit, uninformative parameters, etc.


### DATA ###


### Create grid data directory of blocks x response subsets for fx lookup

data_dir <- list(
  RLL_col = modcol_rll_z_df,
  RLL_ext = modext_rll_z_df
)

covar_dir <- list(
  climate = climate_covs_reduced,
  land    = land_covs_reduced
)
