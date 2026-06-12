### Script: 5a_ModelSpecs_ENV.R
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
#  left_join(
#    block_centroids_df %>%
#      dplyr::select(atlas_block, lat, lon, x_km, y_km),
#    by = "atlas_block"
#    )
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


# STORING #
pa_results_df <- data.frame()

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
chr_covs_all <- c("atlas_block", "common_name", "alpha_code", "transition_state", "guild")

stable_covs_all <- c("lon", "lat", "sr_diff", "priority_block", "pa_prop", # total prop pa per block 
                     
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
  
  forest = c("developed_total_base", "grass_pasture_base",
             "forest_deciduous_base", "forest_mixed_base", "forest_evergreen_base", 
             "wetlands_woody_base", "wetlands_herb_base",
             
             "forest_total_diff", "wetlands_total_diff"),
  
  grass = c("developed_total_base","forest_total_base", "cropland_base",
            "grassland_base", "pasture_base",
            
            "grassland_diff", "pasture_diff"),
  
  marsh = c("water_open_base", "grass_pasture_base",
            "developed_total_base", "forest_total_base", 
            "wetlands_woody_base", "wetlands_herb_base",
            
            "grass_pasture_diff", "wetlands_total_diff"), 
  
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


# Helper: pull spp-guild-covars combinations
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


# Apply: Species Specific Covariate Sets
chr_covs_reduced <- c("atlas_block", "transition_state")
stable_covs_reduced <- c("sr_diff", "priority_block") # no PA in this script
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



# Output Covariates
#guild_key
#land_covs_reduced <- c()

#climate_covs_reduced # "tmax_38yr", "tmax_diff", "prcp_38yr", "tmin_diff", "prcp_diff"
#climate_covs_reduced <- c()

#numeric_covs_reduced <- c(land_covs_reduced, climate_covs_reduced, stable_covs_reduced)

### repeat above if needed 




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

effort_covars <- c("sr_diff", "priority_block")



### Create grid data directory of blocks x response subsets for fx lookup
partition_grid <- expand.grid(
  response = c("col", "ext"),
  blocks = "RLL",
  partition = c("climate", "land"),
  stringsAsFactors = FALSE
)

data_dir <- list(
  RLL_col = modcol_rll_z_df,
  RLL_ext = modext_rll_z_df
)

covar_dir <- list(
  climate = climate_covs_reduced,
  land    = land_covs_reduced
)



### --- STEP 1: ENV PARTITIONED MODELS --- ###

### Full, additive candidate model sets for climate and land cover covariates
# w/ null included. AICc fitting, ranking; extract covars to pass on from "reference" 
# model (ie. delta = 0) to global model selection steps.


### FUNCTIONS ###

### Helper: generate all additive predictor combos w/in partitioned covariate sets
# required = structural effort control variable, PA% as primary effect of interest
BuildPartitionedRHS <- function(covariates, required) {
  
  covariates <- unique(covariates)
  required   <- unique(required)
  
  optional <- setdiff(covariates, required)
  
  # Case: required only
  rhs <- paste(required, collapse = " + ")
  
  # Add combinations of optional covariates
  if (length(optional) > 0) {
    rhs_optional <- unlist(
      lapply(seq_along(optional), function(k) {
        combn(optional, k, FUN = function(x) {
          paste(c(required, x), collapse = " + ")
        })
      }),
      
      use.names = FALSE
    )
    
    rhs <- c(rhs, rhs_optional)
  }
  
  unique(rhs)
}


### Helper: Fit, rank partitioned candidate models
FitPartitionedModels <- function(response, 
                                 data, 
                                 covariates, 
                                 required = effort_covars,
                                 family = binomial,
                                 include_null = TRUE) {
  
  # Predictor structure
  rhs_terms <- BuildPartitionedRHS(
    covariates = covariates,
    required = required
  )
  
  formulas <- lapply(rhs_terms, function(rhs) {
    as.formula(paste(response, "~", rhs))
  })
  
  modnames <- rhs_terms
  
  # Null structure
  if (include_null) {
    formulas <- c(list(as.formula(paste(response, "~ 1"))), formulas)
    modnames <- c("NULL", rhs_terms)
  }
  
  # Fit models
  models <- lapply(formulas, function(f) {
    glm(f, data = data, family = family)
  })
  
  names(models) <- modnames
  
  MuMIn::model.sel(models, rank = "AICc")
}


# Apply: Index, fit partitions in single loop function
partitioned_add_models <- lapply(seq_len(nrow(partition_grid)), function(i) {
  
  row <- partition_grid[i, ]
  key <- paste(row$blocks, row$response, sep = "_")
  
  
  FitPartitionedModels(
    response     = row$response,
    data         = data_dir[[key]],
    covariates   = covar_dir[[row$partition]],
    family       = binomial,
    include_null = TRUE
  )
})

# Names in AICc table
names(partitioned_add_models) <- with(
  partition_grid,
  paste(response, blocks, partition, sep = "_")
)


### Helper: Extract top models (threshold: delta < 2)
ExtractTopModels <- function(model_sel, delta = 2, drop_null = TRUE) {
  
  df <- as.data.frame(model_sel)
  df$Modnames <- rownames(df)
  df <- df[df$delta <= delta, , drop = FALSE]
  
  if (drop_null) {
    df <- df[df$Modnames != "NULL", , drop = FALSE]
  }
  
  df
}


# Apply: threshold to all candidate models 
top_part_models <- lapply(partitioned_add_models, ExtractTopModels, delta = 2)



### Helper: Extract [REFERENCE] model [delta = lowest value]
ExtractReferenceModel <- function(model_sel_df) {
  
  if (nrow(model_sel_df) == 0) return(model_sel_df)
  
  min_delta <- min(model_sel_df$delta)
  ref <- model_sel_df[model_sel_df$delta == min_delta, , drop = FALSE]
  
  if (nrow(ref) != 1) {
    warning(
      "Expected 1 reference model, found ", nrow(ref),
      ". Using first."
    )
    ref <- ref[1, , drop = FALSE]
  }
  
  ref
}


# Apply: reference model extraction
reference_part_models <- lapply(top_part_models, ExtractReferenceModel)

### Helper: Extract [REFERENCE] model covariates
ExtractReferenceCovariates <- function(model_string) {
  
  strsplit(model_string, "\\+")[[1]] %>%
    trimws
}

# Apply: covariate extraction from ref model
reference_part_covariates <- lapply(reference_part_models, function(df) {
  
  if (nrow(df) == 0) return(character(0))
  
  ExtractReferenceCovariates(df$Modnames[1])
})



### Helper: Merge partitioned ref mod climate, land covar lists
MergePartitionedCovariates <- function(ref_covariate_list) {
  
  # Full index names
  names_dir <- data.frame(
    full_name = names(ref_covariate_list),
    stringsAsFactors = FALSE
  )
  
  # Extract index components by name (response_blocks_partition)
  components_dir <- do.call(
    rbind,
    strsplit(names_dir$full_name, "_")
  ) 
  
  colnames(components_dir) <- c("response", "blocks", "partition")
  names_dir <- cbind(names_dir, components_dir)
  combos <- unique(names_dir[, c("blocks", "response")])
  
  merged <- lapply(seq_len(nrow(combos)), function(i) {
    
    blocks <- combos$blocks[i]
    resp <- combos$response[i]
    
    idx <- names_dir$blocks   == blocks &
      names_dir$response == resp
    
    covs <- unlist(ref_covariate_list[idx])
    
    if (length(covs) == 0) { # fail-safe
      warning("No covariates for ", blocks, "_", resp)
      return(character(0))
    }
    unique(covs)
  })
  names(merged) <- paste(combos$blocks, combos$response, sep = "_")
  
  merged
}


# Apply: to partitioned covar sets
merged_ref_covariates <- MergePartitionedCovariates(reference_part_covariates)




### --- STEP 2: ENV MODEL SELECTION --- ###

### Combine climate, land covars lists from partitioned models into single covar list/pool for 
# new additive mod selection process; find new ref model to carry into interaction step (2B).

effort_covars <- c("sr_diff", "priority_block")


merged_ref_covariates <- lapply(
  merged_ref_covariates,
  function(covs) unique(c(covs))
)



BuildGlobalRHS <- function(covariates,
                           forced = effort_covars) {
  
  covariates <- unique(covariates)
  
  # Ensure forced is present
  forced <- intersect(forced, covariates)
  
  # Optional = everything except forced
  optional <- setdiff(covariates, forced)
  
  rhs <- character(0)
  
  # Forced-only model
  rhs <- c(rhs, paste(forced, collapse = " + "))
  
  # Forced + subsets of optional
  if (length(optional) > 0) {
    rhs_optional <- unlist(
      lapply(seq_along(optional), function(k) {
        combn(optional, k, FUN = function(x) {
          paste(c(forced, x), collapse = " + ")
        })
      }),
      use.names = FALSE
    )
    rhs <- c(rhs, rhs_optional)
  }
  
  unique(rhs)
}



FitGlobalModels <- function(response,
                            data,
                            covariates,
                            family = binomial,
                            include_null = TRUE) {
  
  rhs_terms <- BuildGlobalRHS(
    covariates = covariates,
    forced     = effort_covars
  )
  
  models <- list()
  
  # NULL model
  if (include_null) {
    models[["NULL"]] <- glm(
      as.formula(paste(response, "~ 1")),
      data = data,
      family = family
    )
  }
  
  # Additive models
  for (rhs in rhs_terms) {
    models[[rhs]] <- glm(
      as.formula(paste(response, "~", rhs)),
      data = data,
      family = family
    )
  }
  
  model.sel(models)
}


# Filter out response x block models w 0 covars
valid_keys <- names(merged_ref_covariates)[
  lengths(merged_ref_covariates) > 0
]


global_models <- lapply(valid_keys, function(key) { # or lapply(merged_left_covariates){}
  
  FitGlobalModels(
    response   = strsplit(key, "_")[[1]][2],
    data       = data_dir[[key]],
    covariates = merged_ref_covariates[[key]],
    family     = binomial
  )
})

names(global_models) <- valid_keys

top_global_models <- lapply(global_models, ExtractTopModels, delta = 2)
reference_global_models <- lapply(top_global_models, ExtractReferenceModel) 





### --- STEP 3: ENV TOP CANDIDATE MODEL FITTING --- ###

### Helper: Pull formula for ref mod
ExtractRefFormula <- function(ref_df, response) {
  
  if (is.null(ref_df) || nrow(ref_df) == 0) {
    stop("Empty reference model for response: ", response)
  }
  
  as.formula(
    paste(response, "~", ref_df$Modnames[1])
  )
}



global_glm_models <- lapply(names(reference_global_models), function(nm) { 
  
  # nm format: "_col" v "_ext"
  response <- strsplit(nm, "_")[[1]][2]
  
  glm(
    formula = ExtractRefFormula(reference_global_models[[nm]], response),
    data    = data_dir[[nm]],
    family  = binomial
  )
})

names(global_glm_models) <- names(reference_global_models)

global_glm_summaries <- lapply(global_glm_models, summary)
global_glm_summaries





########################
### MODEL ASSESSMENT ###
########################

models <- global_glm_models


### --- 1: DATA STATS CHECKS --- ###

### --- 1A: Correlated Pairs
GetCorrelatedPairs <- function(model, cutoff = 0.7, digits = 2) {
  
  mm <- model.matrix(model)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  
  cor_mat <- cor(mm, use = "complete.obs")
  
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("var1", "var2", "correlation")
  
  cor_df %>%
    filter(var1 != var2) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(var1, var2)), collapse = "__")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    dplyr::select(-pair) %>%
    filter(abs(correlation) >= cutoff) %>%
    mutate(correlation = round(correlation, digits)) %>%
    arrange(desc(abs(correlation)))
}

corr_pairs <- lapply(models, GetCorrelatedPairs, cutoff = 0.7)
#corr_pairs



### --- 1B: VIF
vif_mods <- lapply(models, car::vif)
#vif_mods



### --- 2: MODEL PERFORMANCE --- ###

### --- 2A: AUC
auc_mods <- bind_rows(lapply(names(models), function(nm) {
  mod <- models[[nm]]
  
  data.frame(
    model = nm,
    AUC = as.numeric(
      pROC::auc(
        mod$model[[1]],
        predict(mod, type = "response")
      )
    )
  )
}))
#auc_mods


### --- 2B: McFadden's Pseudo-R2
r2_mods <- bind_rows(lapply(names(models), function(nm) {
  r2 <- pscl::pR2(models[[nm]])
  data.frame(
    model = nm,
    McFadden_R2 = r2["McFadden"]
  )
}))
#r2_mods


### --- 2C: Over-dispersion
### disp. ratio approx. = 1 (none); disp. ration > 1.2 (possible); p < 0.05 (sig. over-dispersion)
overdisp_mods <- bind_rows(lapply(names(models), function(nm) {
  chk <- performance::check_overdispersion(models[[nm]])
  
  tibble(
    model = nm,
    dispersion_ratio = chk$dispersion_ratio,
    p_value = chk$p_value
  )
}))
#overdisp_mods



### --- 3: PARAMETER INFERENCE --- ###

# Arnold (2010) uninformative parameters:
# 1) 95% CI overlaps 0
# 2) LRT/Wald p-value nonsig.
# 3) delta_AICc negligible with removal of term (approx. < 2)


### --- 3A: Coefficient Summary
coefs_summary <- bind_rows(lapply(names(models), function(nm) {
  broom::tidy(models[[nm]], conf.int = TRUE) %>%
    mutate(model = nm)
}))
#coefs_summary


### --- 3B: Flag Potentially Uninformative Parameters
# flag parameters where 95% CI overlaps 0, p > 0.05
uninf_pars <- coefs_summary %>%
  filter(term != "(Intercept)") %>%
  mutate(
    ci_overlaps_zero = conf.low <= 0 & conf.high >= 0,
    non_significant = p.value > 0.05,
    potentially_uninformative = ci_overlaps_zero & non_significant
  ) %>%
  arrange(model, desc(potentially_uninformative))
#uninf_pars


### --- 3C: Summarize 
uninf_summary <- uninf_pars %>%
  group_by(model) %>%
  summarise(
    n_predictors = n(),
    n_uninformative = sum(potentially_uninformative, na.rm = TRUE),
    pct_uninformative = 100 * n_uninformative / n_predictors,
    .groups = "drop"
  )
#uninf_summary


### --- 3D: Likelihood Ratio Tests
# evaluate whether removing parameter worsens model fit
lrt_results <- bind_rows(lapply(names(models), function(nm) {
  drop1(models[[nm]], test = "Chisq") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    mutate(model = nm)
}))
#lrt_results


### --- 3E: Visualize Coefficient Estimates
coefs_summary %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), width = 0.15) +
  facet_wrap(~ model, scales = "free_y") +
  labs(
    x = "Coefficient Estimate",
    y = NULL,
    title = "Global Model Effect Estimates"
  ) +
  theme_minimal(base_size = 13)



### --- 4: RESIDUAL DIAGNOSTICS --- ###

set.seed(123)

resids_dharma <- lapply(models, DHARMa::simulateResiduals, plot = FALSE) # DHARMa std quantile sim residuals


### --- 4A: Residual Uniformity
### Checks whether residuals follow expected uniform distribution;
# p > 0.05 (acceptable distribution), p < 0.05 (possible model misspecification) 
resids_uniformity <- bind_rows(lapply(names(resids_dharma), function(nm) {
  tst <- DHARMa::testUniformity(resids_dharma[[nm]], plot = FALSE)
  data.frame(model = nm, statistic = tst$statistic, p_value = tst$p.value)
}))
#resids_uniformity



### --- 4C: Residual Dispersion
### p > 0.05 (dispersion acceptable), p < 0.005 (possible over- or under-dispersion)
resids_dispersion <- bind_rows(lapply(names(resids_dharma), function(nm) {
  tst <- DHARMa::testDispersion(resids_dharma[[nm]], plot = FALSE)
  data.frame(model = nm, statistic = tst$statistic, p_value = tst$p.value)
}))
#resids_dispersion


### --- 4D: Residual Outlier Test
resids_outlier <- bind_rows(lapply(names(resids_dharma), function(nm) {
  tst <- DHARMa::testOutliers(resids_dharma[[nm]], plot = FALSE)
  
  tibble(
    model = nm,
    outlier_frequency = if (!is.null(tst$estimate)) as.numeric(tst$estimate) else NA_real_,
    p_value = tst$p.value %||% NA_real_,
    method = tst$method %||% NA_character_
  )
}))
#resids_outlier



### --- 4E: Diagnostics Results
corr_pairs
vif_mods
auc_mods
r2_mods
overdisp_mods
coefs_summary
uninf_pars
uninf_summary
lrt_results
resids_uniformity
resids_dispersion
resids_outlier



### --- 5: Visualizations --- ###

labels <- c(
  RLL_col = "Colonization (RLL_col)",
  RLL_ext = "Extirpation (RLL_ext)"
)


### --- 5A: DHARMa resids QQ, pred probs

# Residual summary plots
par(mfrow = c(1, 2),
    mar = c(4, 4, 3, 1))

for (nm in names(resids_dharma)) {
  plot(
    resids_dharma[[nm]],
    main = labels[[nm]]
  )
}

par(mfrow = c(1, 1))


# QQ plots
par(mfrow = c(1, 2),
    mar = c(4, 4, 4, 1))

for (nm in names(resids_dharma)) {
  DHARMa::plotQQunif(
    resids_dharma[[nm]],
    main = labels[[nm]]
  )
}

par(mfrow = c(1, 1))


# Dispersion
par(mfrow = c(1, 2),
    mar = c(4, 4, 4, 1))

for (nm in names(resids_dharma)) {
  DHARMa::testDispersion(
    resids_dharma[[nm]],
    plot = TRUE
  )
  title(main = labels[[nm]], line = 2)
}

par(mfrow = c(1, 1))


# Outliers
par(mfrow = c(1, 2),
    mar = c(4, 4, 4, 1))

for (nm in names(resids_dharma)) {
  DHARMa::testOutliers(
    resids_dharma[[nm]],
    plot = TRUE
  )
  title(main = labels[[nm]], line = 2)
}

par(mfrow = c(1, 1))






### --- 5B: ROC Curves
pred_probs <- lapply(models, predict, type = "response")

roc_col <- pROC::roc(models$RLL_col$model[[1]], pred_probs$RLL_col)
roc_ext <- pROC::roc(models$RLL_ext$model[[1]], pred_probs$RLL_ext)

roc_df <- bind_rows(
  data.frame(
    FPR = 1 - roc_col$specificities,
    TPR = roc_col$sensitivities,
    Model = "Colonization"
  ),
  data.frame(
    FPR = 1 - roc_ext$specificities,
    TPR = roc_ext$sensitivities,
    Model = "Extirpation"
  )
)

auc_col <- as.numeric(pROC::auc(roc_col))
auc_ext <- as.numeric(pROC::auc(roc_ext))


# plot
ggplot(roc_df, aes(FPR, TPR, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Colonization" = "#2C7BB6",
                                "Extirpation" = "#D7191C"),
                     labels = c(
                       paste0("Colonization (AUC = ", round(auc_col, 3), ")"),
                       paste0("Extirpation (AUC = ", round(auc_ext, 3), ")")
                     )) +
  labs(
    title = "ROC Curves for ENV Logistic GLMs",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))
