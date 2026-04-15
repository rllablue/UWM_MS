##################### Script with only RLL block modeling ######################

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
  psych, 
  usdm,
  
  # Visualization
  ggplot2,
  viridis,
  gt,
  webshot2
  
)

options(na.action = "na.fail") # required for MuMIn




### --- DATA --- ###

# Upload, Wrangle Data #
spp_zf_rll <- read.csv("data/summaries/spp_zf_rll.csv") # df
spp_list <- unique(spp_zf_rll$common_name) # vector

covars_raw_rll <- read.csv("data/summaries/covars_raw_rll.csv") # df
covars_z_rll <- covars_raw_rll %>% # z-standardized covariate values
  mutate(across(
    .cols = -atlas_block,
    .fns = ~ as.numeric(scale(.))
  ))

wibba_summary_rll <- read.csv("data/summaries/wibba_summary_rll.csv") # df
blocks_rll <- wibba_summary_rll$atlas_block # vector



# Create Guilds #
# spp_guilds <- tibble(
#  common_name = unique(spp_zf_rll$common_name),
#  guild = NA_character_
#)

# write.csv(spp_guilds, "data/summaries/species_guilds.csv", row.names = FALSE)
spp_guilds <- read.csv("data/summaries/species_guilds.csv")

spp_zf_rll <- spp_zf_rll %>%
  left_join(spp_guilds, by = "common_name")



# FULL MODELING DF #
# all zero-filled species, z-scaled covariate data
mod_data_all <- spp_zf_rll %>%
  left_join(covars_z_rll, by = "atlas_block")
# write.csv(mod_data_all, "outputs/data/mod_data_all.csv", row.names = FALSE)



#####################
### MODELING DATA ###
#####################

# FULL DATASET #
mod_data_all <- read.csv("outputs/data/mod_data_all.csv")

# STORING #
pa_results_df <- data.frame()


# --- SPECIES-RESPONSE SUBSETTING --- #

### Separate data into bins where A2 detection is either 1 or 0; not necessarily
# in precise probabilities between col, per, abs, ext, but more what promotes
# 'new' v. 'continued' colonization, ie. what promotes det = 1, as opposed to 
# what promotes det = 0. Colonization = col + abs data, Extinction = ext + pre data
### Scaling of covars w/in subsets for relevant normalized values

# Species to model
spp_alpha <- "CAWA"
spp_name <- "Canada Warbler"


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
      !!response_name := ifelse(transition_state == response_one, 1, 0)
    )
}




# Apply: create species-state-specific dfs

# Col-RLL
mod_col_rll <- BuildSppRespDfs(
  data = mod_data_all,
  species = spp_alpha,
  block_vector = blocks_rll,
  state_pairs = c("Colonization", "Absence"),
  response_name = "col",
  response_one = "Colonization"
)

# Ext-RLL
mod_ext_rll <- BuildSppRespDfs(
  data = mod_data_all,
  species = spp_alpha,
  block_vector = blocks_rll,
  state_pairs = c("Extinction", "Persistence"),
  response_name = "ext",
  response_one = "Extinction"
)



# Store in list
mod_dfs_all <- list(
  col_rll = mod_col_rll,
  ext_rll = mod_ext_rll
)


# DATA SUMMARIES #
### Counts of each response type per species (247) and block set (3337, 858)
rll_spp_counts <- mod_data_all %>%
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

rll_valid_counts <- mod_data_all %>%
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
    ext_per_total = Persistence + Extinction,
    
    # event structure
    col_events     = Colonization,
    col_nonevents  = Absence,
    ext_events     = Extinction,
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
    Extinction,
    Persistence,
    col_abs_total,
    ext_per_total,
    col_valid,
    ext_valid
  )




##################
### COVARIATES ###
##################

# FULL DATASET #
factor_covs_all <- c("atlas_block", "common_name", "alpha_code", "transition_state", "guild")

stable_covs_all <- c("lon", "lat", "sr_diff", "pa_prop", # total prop pa per block 
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

### ECOLOGICAL FILTERS ###
### Habitat guilds to bin land cover covariates

# Inputs
guild_key <- list(
  
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


GetGuildCovs <- function(species, data, guild_map) {
  guild_value <- data %>%
    filter(alpha_code == species) %>%
    distinct(guild) %>%
    pull(guild)
  
  if (length(guild_value) == 0) {
    stop("Species lacks guild assignment.")
  }
  
  if (!guild_value %in% names(guild_map)) {
    stop("Guild not found in key.")
  }
  
  guild_map[[guild_value]]
}


# Subset outputs
factor_covs_reduced <- c("atlas_block", "transition_state")
stable_covs_reduced <- c("sr_diff") # no pa in this stage
land_covs_reduced <- GetGuildCovs(spp_alpha, mod_data_all, guild_key)
climate_covs_reduced <- c("tmax_38yr", "prcp_38yr", # base year values
                          
                          "tmax_diff", "tmin_diff", "prcp_diff") # change values

numeric_covs_reduced <- c(stable_covs_reduced, land_covs_reduced, climate_covs_reduced)



### STATISTICAL FILTERS ### 

# CORRELATIONS #
### Pair-wise, linear relationships (i.e. bivariate)

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

corrs1 <- GetHighCorrs(mod_col_rll, numeric_covs_reduced)
corrs2 <- GetHighCorrs(mod_ext_rll, numeric_covs_reduced)



# Output covariates
# stable_covars_reduced

guild_key
land_covs_reduced <- c("developed_total_base", "grass_pasture_base",
                       "forest_deciduous_base", "forest_mixed_base", "forest_evergreen_base", 
                       "wetlands_woody_base", "wetlands_herb_base",
                       
                       "forest_total_diff", "wetlands_total_diff")

climate_covs_reduced # "tmax_38yr", "tmax_diff", "prcp_38yr", "tmin_diff", "prcp_diff"
climate_covs_reduced <- c("tmax_38yr", "tmax_diff", "tmin_diff", "prcp_diff")

numeric_covs_reduced <- c(land_covs_reduced, climate_covs_reduced, stable_covs_reduced)



# VARIANCE INFLATION FACTOR #
### Single variable relationship to all others as a group (i.e. multivariate); 
# multicollinearity redundancy measure, i.e. too similar to uniquely est coefficients (ideally VIF < 5)
responses <- c(
  col_rll = "col",
  ext_rll = "ext"
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
guild_key
land_covs_reduced <- c("developed_total_base", 
                       "forest_deciduous_base", "forest_mixed_base",
                       "wetlands_woody_base", "wetlands_herb_base",
                       
                       "forest_total_diff", "wetlands_total_diff")

climate_covs_reduced # "tmax_38yr", "tmax_diff", "prcp_38yr", "tmin_diff", "prcp_diff"
climate_covs_reduced <- c("tmax_38yr")

numeric_covs_reduced <- c(land_covs_reduced, climate_covs_reduced, stable_covs_reduced)



# Repeat
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





#######################
### MODEL SELECTION ### 
#######################

### --- PROCESS, SET UP --- ###

### Steps ###

# 1) Partitioned AICc Model Selection
# Construct, rank separate models partitioning land and climate covariates
# Carry over best covariates for each set into interactions, global model


# 2) Global AICc Model Selection

# 2A) Additive only models
# Construct, fit, rank top models AICc < 2 w/ addiitve terms only
# Extract reference model for each block/response subset to carry into interactions

# 2B) Interaction Terms
# Generate plausible 2-way interactions w/ protected area using terms from additive
# reference model for each block/response subset


# 3) Model Diagnostics
# AICc Table visualization
# Assess reference models for each block/response for uninformative parameters, etc.


### DATA ###

effort_covar <- "sr_diff"



### Create grid data directory of blocks x response subsets for fx lookup
partition_grid <- expand.grid(
  response = c("col", "ext"),
  blocks = "RLL",
  partition = c("climate", "land"),
  stringsAsFactors = FALSE
)

data_dir <- list(
  RLL_col = mod_col_rll,
  RLL_ext = mod_ext_rll
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
                                 required = effort_covar,
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




### --- STEP 2: GLOBAL ENV MODELS --- ###

### Combine climate, land covars lists from partitioned models into single covar list/pool for 
# new additive mod selection process; find new ref model to carry into interaction step (2B).

effort_covar <- "sr_diff"


merged_ref_covariates <- lapply(
  merged_ref_covariates,
  function(covs) unique(c(covs)) # WIP: removed pa_covar
)



BuildGlobalRHS <- function(covariates,
                           forced = effort_covar) { # WIP: removed pa_covar
  
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
    forced     = effort_covar
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





### --- STEP 3: ENV MODELING --- ###

# Model data (responses: col, abs)
data_dir <- list(
  RLL_col = mod_col_rll,
  RLL_ext = mod_ext_rll
)


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

vif_global_models <- lapply(global_glm_models, car::vif)
vif_global_models





### --- 4: GLOBAL PA MODELS, MODELING --- ###

pa_covar <- c("pa_prop")



### BASE PA ###

ExtractPAFormula <- function(ref_df, response) {
  
  if (is.null(ref_df) || nrow(ref_df) == 0) {
    stop("Empty reference model for response: ", response)
  }
  
  as.formula(
    paste(response, "~", ref_df$Modnames[1], "+ pa_prop")
  )
}


pa_glm_models <- lapply(names(reference_global_models), function(nm) {
  
  response <- strsplit(nm, "_")[[1]][2]
  
  glm(
    formula = ExtractPAFormula(reference_global_models[[nm]], response),
    data    = data_dir[[nm]],
    family  = binomial
  )
})

names(pa_glm_models) <- names(reference_global_models)


pa_glm_summaries <- lapply(pa_glm_models, summary)
pa_glm_summaries

vif_pa_models <- lapply(pa_glm_models, car::vif)
vif_pa_models




# PA: Interactions
mod_col_pa <- pa_glm_models[["RLL_col"]]
mod_ext_pa <- pa_glm_models[["RLL_ext"]]

form_col_base <- formula(mod_col_pa)
form_ext_base <- formula(mod_ext_pa)

form_col_int <- update(form_col_base, . ~ . + pa_prop:forest_total_base) # :tmax_38yr, :forest_total_base, :grassland_base, :wetlands_total_base
form_ext_int <- update(form_ext_base, . ~ . + pa_prop:forest_total_base)


mod_col_pa_int <- glm(
  formula = form_col_int,
  data    = mod_col_rll,
  family  = binomial
)
summary(mod_col_pa_int)
car::vif(mod_col_pa_int, type = "terms")


mod_ext_pa_int <- glm(
  formula = form_ext_int,
  data    = mod_ext_rll,
  family  = binomial
)
summary(mod_ext_pa_int)
car::vif(mod_ext_pa_int, type = "terms")






# Extract, Collect PA effects 
ExtractPACoefficients <- function(model, model_name, spp_name) {
  
  parts <- strsplit(model_name, "_")[[1]]
  block    <- parts[1]
  response <- parts[2]
  
  df <- broom::tidy(model)
  
  if (!"pa_prop" %in% df$term) {
    return(NULL)
  }
  
  df %>%
    filter(term == "pa_prop") %>%
    mutate(
      species   = spp_name,
      response  = response,
      block     = block,
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      sig       = ifelse(conf.low > 0 | conf.high < 0, "yes", "no")
    ) %>%
    dplyr::select(species, block, response, estimate, std.error, conf.low, conf.high, p.value, sig)
}


species_results <- dplyr::bind_rows(
  lapply(names(pa_glm_models), function(nm) {
    ExtractPACoefficients(pa_glm_models[[nm]], nm, spp_name)
  })
)


if (!"species" %in% names(rll_results_df)) {
  
  rll_results_df <- species_results
  
} else {
  
  rll_results_df <- rll_results_df %>%
    dplyr::filter(species != spp_name) %>%
    dplyr::bind_rows(species_results)
  
}




write.csv(rll_results_df,"outputs/data/basepa_effects_rll.csv")
rll_results_df <- read.csv("outputs/data/basepa_effects_rll.csv")




# Caterpillar Plot

caterpillar_df <- rll_results_df %>%
  group_by(species, response) %>%
  mutate(strength = max(abs(estimate), na.rm = TRUE)) %>%
  ungroup()


PlotEffects <- function(df) {
  
  dodge <- position_dodge(width = 0.5)
  
  ggplot(
    df,
    aes(x = estimate,
        y = reorder(species, strength),
        color = response,
        group = response)
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
    
    scale_color_manual(values = c("col" = "#1f78b4", "ext" = "#e31a1c")) +
    
    labs(
      title = "Effect of PA (RLL only)",
      x = "PA Effect (log-odds)",
      y = "Species",
      color = "Response",
      shape = "Significant"
    ) +
    
    theme_minimal()
}


plot_rll <- PlotEffects(caterpillar_df)
plot_rll









### GAP STATUS #################################################################

gap_covars <- c("gap1_prop", "gap2_prop", "gap3_prop")




MakeSubsets <- function(vars) {
  
  unlist(
    lapply(0:length(vars), function(k) {
      combn(vars, k, simplify = FALSE)
    }),
    recursive = FALSE
  )
}


FitPAGapModels <- function(ref_df, response, data) {
  
  # Base environmental structure
  env_rhs <- ref_df$Modnames[1]
  
  # Base PA model (always included)
  base_rhs <- paste(env_rhs, "+ pa_prop")
  
  # All GAP combinations
  gap_sets <- MakeSubsets(gap_covars)
  
  # Build RHS formulas
  rhs_list <- lapply(gap_sets, function(gaps) {
    if (length(gaps) == 0) {
      base_rhs
    } else {
      paste(base_rhs, "+", paste(gaps, collapse = " + "))
    }
  })
  
  # Name models nicely
  names(rhs_list) <- sapply(gap_sets, function(gaps) {
    if (length(gaps) == 0) return("PA_only")
    paste("PA", paste(gaps, collapse = "_"), sep = "_")
  })
  
  # Fit models
  models <- lapply(rhs_list, function(rhs) {
    glm(as.formula(paste(response, "~", rhs)),
        data = data,
        family = binomial)
  })
  
  # Model selection
  ms_table <- MuMIn::model.sel(models, rank = "AICc")
  
  list(models = models, sel_table = ms_table)
}



pa_sig_models <- c("RLL_col")

pa_gap_models <- lapply(pa_sig_models, function(nm) {
  
  response <- strsplit(nm, "_")[[1]][2]
  
  FitPAGapModels(
    ref_df = reference_global_models[[nm]],
    response = response,
    data = data_dir[[nm]]
  )
})

names(pa_gap_models) <- pa_sig_models



RefitTopModel <- function(model_obj) {
  
  ms_table <- model_obj$sel_table
  models   <- model_obj$models
  
  df <- as.data.frame(ms_table)
  df$Modnames <- rownames(df)
  
  top_name <- df$Modnames[which.min(df$delta)][1]
  
  models[[top_name]]
}




top_pa_gap_models <- lapply(pa_gap_models, RefitTopModel)

# Check results
summary(top_pa_gap_models$RLL_col)
car::vif(top_pa_gap_models$RLL_col)


pa_vars <- c("pa_prop", "gap1_prop", "gap2_prop", "gap3_prop")

cor_mat <- cor(
  data_dir$RLL_col[, pa_vars],
  use = "complete.obs"
)

round(cor_mat, 3)

corrplot::corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.cex = 0.8,
  addCoef.col = "black"
)



pa_only <- subset(data_dir$RLL_col, pa_prop > 0)

cor_pa_only <- cor(
  pa_only[, pa_vars],
  use = "complete.obs"
)

round(cor_pa_only, 3)

corrplot::corrplot(
  cor_pa_only,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.cex = 0.8,
  addCoef.col = "black"
)


pairs(
  data_dir$RLL_col[, pa_vars],
  pch = 16,
  cex = 0.5
)




# GAP-only Mod
gap_only_mod <- glm(
  formula = col ~ sr_diff + tmax_38yr + developed_total_base + forest_deciduous_base + forest_mixed_base + 
    gap1_prop + gap2_prop + gap3_prop,
  data    = mod_col_rll,
  family  = binomial)

summary(gap_only_mod)






mod_vars <- c("sr_diff", "tmax_38yr", "developed_total_base", "forest_deciduous_base", "forest_mixed_base", "pa_prop",
              "gap1_prop", "gap2_prop", "gap3_prop")

cor_mat_full <- cor(
  data_dir$RLL_col[, mod_vars],
  use = "complete.obs"
)

round(cor_mat_full, 3)

corrplot::corrplot(
  cor_mat_full,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.cex = 0.8,
  addCoef.col = "black"
)







### MANAGER ####################################################################

own_covars <- c("sdnr_prop", "usfs_prop", "tnc_prop") # "nas_prop"


FitPAModels <- function(ref_df, response, data) {
  
  # Extract environmental RHS from reference model
  env_rhs <- ref_df$Modnames[1]
  
  # Build candidate RHS strings
  rhs_list <- list(
    ENV = env_rhs,
    ENV_PA = paste(env_rhs, "+", paste(pa_covars, collapse = " + ")),
    ENV_PA_OWN = paste(env_rhs, "+",
                       paste(c(pa_covars, own_covars), collapse = " + "))
  )
  
  # Fit models
  models <- lapply(rhs_list, function(rhs) {
    glm(as.formula(paste(response, "~", rhs)),
        data = data,
        family = binomial)
  })
  
  # Model selection table
  ms_table <- MuMIn::model.sel(models, rank = "AICc")
  
  list(models = models, sel_table = ms_table)
}


RefitTopModel <- function(model_obj) {
  
  # Extract selection table
  ms_table <- model_obj$sel_table
  models   <- model_obj$models
  
  # Determine top model name (Δ = 0)
  df <- as.data.frame(ms_table)
  df$Modnames <- rownames(df)
  top_name <- df$Modnames[df$delta == min(df$delta)][1]
  
  # Return the **fitted model object** directly
  models[[top_name]]
}



pa_models <- lapply(names(reference_global_models), function(nm) {
  response <- strsplit(nm, "_")[[1]][2]
  
  FitPAModels(
    ref_df = reference_global_models[[nm]],
    response = response,
    data = data_dir[[nm]]
  )
})

names(pa_models) <- names(reference_global_models)

# Extract top Δ=0 models
top_pa_models <- lapply(pa_models, RefitTopModel)
names(top_pa_models) <- names(pa_models)


# Check VIFs
vif_results_pa <- lapply(top_pa_models, car::vif)

# Quick summary example
summary(top_pa_models$DNR_col)
summary(top_pa_models$DNR_ext)
summary(top_pa_models$RLL_col)
summary(top_pa_models$RLL_ext)





########################## VISUALIZATIONS ######################################

MakeResponseCurve_fixed <- function(model_col,
                                    model_ext,
                                    focal_var,
                                    data,
                                    focal_seq = seq(0, 1, length.out = 100),
                                    custom_baseline = NULL) {
  
  vars_col <- all.vars(formula(model_col))[-1]
  vars_ext <- all.vars(formula(model_ext))[-1]
  all_vars <- unique(c(vars_col, vars_ext))
  
  base_vals <- list()
  
  for (v in all_vars) {
    
    if (!is.null(custom_baseline) && v %in% names(custom_baseline)) {
      base_vals[[v]] <- custom_baseline[[v]]
      
    } else {
      if (is.numeric(data[[v]])) {
        base_vals[[v]] <- median(data[[v]], na.rm = TRUE)
      } else {
        base_vals[[v]] <- names(sort(table(data[[v]]), decreasing = TRUE))[1]
      }
    }
  }
  
  pred_data <- do.call(expand.grid, base_vals)
  pred_data <- pred_data[rep(1, length(focal_seq)), ]
  pred_data[[focal_var]] <- focal_seq
  
  # --- COLONIZATION ---
  col_link <- predict(model_col, newdata = pred_data, type = "link", se.fit = TRUE)
  
  col_fit  <- plogis(col_link$fit)
  col_lwr  <- plogis(col_link$fit - 1.96 * col_link$se.fit)
  col_upr  <- plogis(col_link$fit + 1.96 * col_link$se.fit)
  
  # --- EXTINCTION ---
  ext_link <- predict(model_ext, newdata = pred_data, type = "link", se.fit = TRUE)
  
  ext_fit  <- plogis(ext_link$fit)
  ext_lwr  <- plogis(ext_link$fit - 1.96 * ext_link$se.fit)
  ext_upr  <- plogis(ext_link$fit + 1.96 * ext_link$se.fit)
  
  data.frame(
    x = focal_seq,
    col_fit = col_fit,
    col_lwr = col_lwr,
    col_upr = col_upr,
    ext_fit = ext_fit,
    ext_lwr = ext_lwr,
    ext_upr = ext_upr
  )
}






# Start with model variables
vars_needed_dnr <- unique(c(
  all.vars(formula(top_pa_models$DNR_col))[-1],
  all.vars(formula(top_pa_models$DNR_ext))[-1]
))

# Add gap & manager variables you want to explore
extra_covs <- c("gap1_prop", "gap2_prop", "gap3_prop", "sdnr_prop", "tnc_prop", "usfs_prop")

# Combine, keeping only columns that exist in your dataset
vars_needed_dnr <- intersect(c(vars_needed_dnr, extra_covs), names(mod_col_dnr))

dnr_all <- rbind(
  mod_col_dnr[, vars_needed_dnr, drop = FALSE],
  mod_ext_dnr[, vars_needed_dnr, drop = FALSE]
)




vars_needed_rll <- unique(c(
  all.vars(formula(top_pa_models$RLL_col))[-1],
  all.vars(formula(top_pa_models$RLL_ext))[-1]
))

vars_needed_rll <- intersect(c(vars_needed_rll, extra_covs), names(mod_col_rll))

rll_all <- rbind(
  mod_col_rll[, vars_needed_rll, drop = FALSE],
  mod_ext_rll[, vars_needed_rll, drop = FALSE]
)








curve_dnr_pa <- MakeResponseCurve_fixed(
  model_col = top_pa_models$DNR_col,
  model_ext = top_pa_models$DNR_ext,
  focal_var = "pa_prop",
  data = dnr_all
)


curve_rll_pa <- MakeResponseCurve_fixed(
  model_col = top_pa_models$RLL_col,
  model_ext = top_pa_models$RLL_ext,
  focal_var = "pa_prop",
  data = rll_all
)



# Text size
title_size       <- 16
subtitle_size    <- 14
legend_title_size <- 12
legend_text_size  <- 11





ggplot(curve_rll_pa, aes(x = x)) +
  
  # Colonization ribbon
  geom_ribbon(aes(ymin = col_lwr, ymax = col_upr, fill = "Colonization"),
              alpha = 0.2) +
  geom_line(aes(y = col_fit, color = "Colonization"),
            size = 1) +
  
  # Extinction ribbon
  geom_ribbon(aes(ymin = ext_lwr, ymax = ext_upr, fill = "Extinction"),
              alpha = 0.2) +
  geom_line(aes(y = ext_fit, color = "Extinction"),
            size = 1) +
  
  scale_color_manual(values = c("Colonization" = "darkorchid",
                                "Extinction" = "orange")) +
  
  scale_fill_manual(values = c("Colonization" = "darkorchid",
                               "Extinction" = "orange")) +
  
  ylim(0, 0.75) +
  labs(x = "Total PA (%)",
       y = "Predicted Probability",
       color = "Process",
       fill  = "Process",
       title = "Red-bellied Woodpecker Response") +
  
  theme_minimal() +
  
  theme(
    
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
    axis.text.y = element_text(size = axis_text_size),
    
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  ) 









