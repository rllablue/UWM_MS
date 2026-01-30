#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

# Core
library(dplyr)
library(tidyverse)
library(magrittr)
library(purrr)
library(stringr)
library(readxl)

# Diagnostics
library(nnet)
library(stats)
library(broom)
library(corrplot)
library(car)
library(performance)
library(DescTools)
library(Metrics)
library(pROC)
library(rsample)
library(lme4)
library(pscl)
library(AICcmodavg)
library(MuMIn)
options(na.action = "na.fail")
library(arm)

# Visualization
library(ggplot2)
library(viridis)
library(gt)
library(webshot2)


### --- DATA --- ###

# Base Data Frames #
spp_zf_rll <- read.csv("data/summaries/spp_zf_rll.csv") # df
spp_list <- unique(spp_zf_rll$common_name) # vector

covars_raw_rll <- read.csv("data/summaries/covars_raw_all.csv") # df
covars_raw_rll <- covars_raw_rll %>% # add spp richness effort proxy
  mutate(sr_diff = sr_Atlas2 - sr_Atlas1,
         grass_pasture_crop_base = grassland_base + pasture_crop_base,
         grass_pasture_crop_diff = grassland_diff + pasture_crop_diff)

wibba_summary_rll <- read.csv("data/summaries/wibba_summary_rll.csv") # df
blocks_rll <- wibba_summary_rll$atlas_block # vector

blocks_dnr <- read_xlsx("data/summaries/CompBlocks_DNR2023.xlsx") # df
blocks_dnr <- blocks_dnr$atlas_block # vector


# Summarize Spp Counts by block
rll_counts <- spp_zf_rll %>%
  filter(atlas_block %in% blocks_rll) %>%
  count(common_name, transition_state) %>%
  pivot_wider(
    names_from  = transition_state,
    values_from = n,
    names_prefix = "rll_",
    values_fill = 0
  )

dnr_counts <- spp_zf_rll %>%
  filter(atlas_block %in% blocks_dnr) %>%
  count(common_name, transition_state) %>%
  pivot_wider(
    names_from  = transition_state,
    values_from = n,
    names_prefix = "dnr_",
    values_fill = 0
  )

spp_list_counts <- rll_counts %>%
  left_join(dnr_counts, by = "common_name") %>%
  mutate(
    total_rll = rowSums(dplyr::select(., starts_with("rll_"))),
    total_dnr = rowSums(dplyr::select(., starts_with("dnr_")))
  )




# Covariates #
factor_covars_all <- c("atlas_block", "common_name", "alpha_code", "transition_state")

stable_covars_all <- c("lon", "lat", "sr_diff", "pa_percent")

land_covars_all <- c("water_open_base", "barren_land_base", "shrub_scrub_base", 
                     "developed_open_base", "developed_low_base", "developed_med_base", "developed_high_base", 
                     "developed_lower_base", "developed_upper_base", "developed_total_base", 
                     "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                     "forest_total_base", "pasture_base", "cropland_base", 
                     "grassland_base", "pasture_crop_base", "grass_pasture_crop_base",
                     "wetlands_woody_base", "wetlands_herb_base", "wetlands_total_base",
                     
                     "water_open_diff", "barren_land_diff", "shrub_scrub_diff", 
                     "developed_open_diff", "developed_low_diff", "developed_med_diff", "developed_high_diff", 
                     "developed_lower_diff", "developed_upper_diff", "developed_total_diff", 
                     "forest_deciduous_diff", "forest_evergreen_diff", "forest_mixed_diff", 
                     "forest_total_diff", "pasture_diff", "cropland_diff", 
                     "grassland_diff", "pasture_crop_diff", "grass_pasture_crop_diff",
                     "wetlands_woody_diff", "wetlands_herb_diff", "wetlands_total_diff")

climate_covars_all <- c("tmax_38yr", "tmin_38yr", "prcp_38yr", 
                        
                        "tmax_diff", "tmin_diff", "prcp_diff")


# New DataFrames #

# Species
spp_name <- "Eastern Meadowlark"


# RLL modeling df
mod_data_rll <- spp_zf_rll %>%
  filter(common_name == spp_name) %>%
  left_join(covars_raw_rll, by = "atlas_block")


# DNR modeling df
mod_data_dnr <- spp_zf_rll %>%
  filter(atlas_block %in% blocks_dnr, 
         common_name == spp_name) %>%  
  left_join(covars_raw_rll, by = "atlas_block")


# write.csv(mod_data_rll, "outputs/data/mod_data_rll.csv", row.names = FALSE)
# write.csv(mod_data_dnr, "outputs/data/mod_data_dnr.csv", row.names = FALSE)



### -- COVARIATE THINNING --- ###

# PRE-SCREENING #
### Screen for biologically relevant covariates (for landcover and climate)
# on a species-, state-specific basis, ie. within-group reduction; run pairwise 
# correlations, VIF to assess multi-collinearity (alt. approach: PCA &/or CLT)
# and further thin predictors


# Species-specific Thinned Covariate Sets

factor_covars_reduced <- c("atlas_block", "transition_state")

stable_covars_all <- c("sr_diff", "pa_percent")

land_covars_reduced <- c("developed_lower_base",
                         "forest_total_base",
                         "grassland_base", "pasture_crop_base",
                         
                         "grassland_diff", "pasture_crop_diff")

covars_numeric_reduced <- c(stable_covars_all, land_covars_reduced, climate_covars_all)


# State-specific Covariate Sets
### Separate data into bins where A2 detection is either 1 or 0; not necessarily
# in precise probabilities between colo, pers, abs, ext, but more what promotes
# 'new' v. 'continued' colonization, ie. what promotes det = 1, as opposed to 
# what promotes det = 0. 
### Scaling of covars w/in subsets for relevant normalized values


# Col-RLL
mod_colabs_rll_z <- mod_data_rll %>%
  filter(transition_state %in% c("Colonization", "Absence")) %>% # filter by state
  dplyr::select(all_of(factor_covars_reduced), # select numeric, factor covars
                all_of(covars_numeric_reduced)) %>%
  mutate(across( # scale only numeric covars
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(col = ifelse(transition_state == "Colonization", 1, 0)) %>% # binomial response variable
  dplyr::select(all_of(factor_covars_reduced), # columns to keep
                ends_with("_z"),
                col)

# Ext-RLL
mod_extper_rll_z <- mod_data_rll %>%
  filter(transition_state %in% c("Extinction", "Persistence")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(ext = ifelse(transition_state == "Extinction", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                ext)

# Per-RLL
mod_perext_rll_z <- mod_data_rll %>%
  filter(transition_state %in% c("Persistence", "Extinction")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(per = ifelse(transition_state == "Persistence", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                per)


# Col-DNR
mod_colabs_dnr_z <- mod_data_dnr %>%
  filter(transition_state %in% c("Colonization", "Absence")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(col = ifelse(transition_state == "Colonization", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                col)

# Ext-DNR
mod_extper_dnr_z <- mod_data_dnr %>%
  filter(transition_state %in% c("Extinction", "Persistence")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(ext = ifelse(transition_state == "Extinction", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                ext)

# Per-DNR
mod_perext_dnr_z <- mod_data_dnr %>%
  filter(transition_state %in% c("Persistence", "Extinction")) %>%
  dplyr::select(all_of(factor_covars_reduced),
                all_of(covars_numeric_reduced)) %>%
  mutate(across(
    .cols = all_of(covars_numeric_reduced),
    .fns = ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  )) %>%
  mutate(per = ifelse(transition_state == "Persistence", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                per)


# CORRELATIONS, COLLINEARITY # 

# Pairwise Correlations
covar_corrs_colabs_rll <- grep("_z$", names(mod_colabs_rll_z), value = TRUE)
covar_corrs_extper_rll <- grep("_z$", names(mod_extper_rll_z), value = TRUE)
covar_corrs_perext_rll <- grep("_z$", names(mod_perext_rll_z), value = TRUE)

covar_corrs_colabs_dnr <- grep("_z$", names(mod_colabs_dnr_z), value = TRUE)
covar_corrs_extper_dnr <- grep("_z$", names(mod_extper_dnr_z), value = TRUE)
covar_corrs_perext_dnr <- grep("_z$", names(mod_perext_rll_z), value = TRUE)


M1 <- cor(mod_colabs_rll_z[, covar_corrs_colabs_rll], use = "pairwise.complete.obs")
M2 <- cor(mod_extper_rll_z[, covar_corrs_extper_rll], use = "pairwise.complete.obs")
M3 <- cor(mod_perext_rll_z[, covar_corrs_perext_rll], use = "pairwise.complete.obs")

M4 <- cor(mod_colabs_dnr_z[, covar_corrs_colabs_dnr], use = "pairwise.complete.obs")
M5 <- cor(mod_extper_dnr_z[, covar_corrs_extper_dnr], use = "pairwise.complete.obs")
M6 <- cor(mod_perext_dnr_z[, covar_corrs_perext_rll], use = "pairwise.complete.obs")


# corrplot::corrplot(M1, method = "color", tl.cex = 0.7, number.cex = 0.6) # visualize correlation plots
# corrplot::corrplot(M2, method = "color", tl.cex = 0.7, number.cex = 0.6)
# corrplot::corrplot(M3, method = "color", tl.cex = 0.7, number.cex = 0.6)

# corrplot::corrplot(M4, method = "color", tl.cex = 0.7, number.cex = 0.6)
# corrplot::corrplot(M5, method = "color", tl.cex = 0.7, number.cex = 0.6)
# corrplot::corrplot(M6, method = "color", tl.cex = 0.7, number.cex = 0.6)


high_corr1 <- which(abs(M1) > 0.7 & abs(M1) < 1, arr.ind = TRUE) # list out high correlation pairs
apply(high_corr1, 1, function(i) cat(rownames(M1)[i[1]], "-", colnames(M1)[i[2]], "r =", M1[i[1],i[2]], "\n"))

high_corr2 <- which(abs(M2) > 0.7 & abs(M2) < 1, arr.ind = TRUE)
apply(high_corr2, 1, function(i) cat(rownames(M2)[i[1]], "-", colnames(M2)[i[2]], "r =", M2[i[1],i[2]], "\n"))

# high_corr3 <- which(abs(M3) > 0.7 & abs(M3) < 1, arr.ind = TRUE)
# apply(high_corr3, 1, function(i) cat(rownames(M3)[i[1]], "-", colnames(M3)[i[2]], "r =", M3[i[1],i[2]], "\n"))

high_corr4 <- which(abs(M4) > 0.7 & abs(M4) < 1, arr.ind = TRUE)
apply(high_corr4, 1, function(i) cat(rownames(M4)[i[1]], "-", colnames(M4)[i[2]], "r =", M4[i[1],i[2]], "\n"))

high_corr5 <- which(abs(M5) > 0.7 & abs(M5) < 1, arr.ind = TRUE)
apply(high_corr5, 1, function(i) cat(rownames(M5)[i[1]], "-", colnames(M5)[i[2]], "r =", M5[i[1],i[2]], "\n"))

# high_corr6 <- which(abs(M6) > 0.7 & abs(M6) < 1, arr.ind = TRUE)
# apply(high_corr6, 1, function(i) cat(rownames(M6)[i[1]], "-", colnames(M6)[i[2]], "r =", M6[i[1],i[2]], "\n"))


# Correlation Thinned Covariate Sets
### Removed one covar from highly correlated pairs when |r| > 0.7
factor_covars_reduced

land_covars_reduced <- c("developed_lower_base",
                         "forest_total_base",
                         "grassland_base", "pasture_crop_base",
                         
                         "grassland_diff", "pasture_crop_diff")

climate_covars_reduced <- c("tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff")

stable_covars_reduced <- c("sr_diff", "pa_percent")

covars_numeric_reduced <- c(land_covars_reduced, climate_covars_reduced, stable_covars_reduced)
covars_numeric_reduced_z <- paste0(covars_numeric_reduced, "_z")



# VIF
vif_model1 <- glm(col ~ ., data = mod_colabs_rll_z[, c("col", covars_numeric_reduced_z)], family = binomial)
vif(vif_model1)
alias(vif_model1)

vif_model2 <- glm(ext ~ ., data = mod_extper_rll_z[, c("ext", covars_numeric_reduced_z)], family = binomial)
vif(vif_model2)
alias(vif_model2)

vif_model3 <- glm(per ~ ., data = mod_perext_rll_z[, c("per", covars_numeric_reduced_z)], family = binomial)
vif(vif_model3)
alias(vif_model3)

vif_model4 <- glm(col ~ ., data = mod_colabs_dnr_z[, c("col", covars_numeric_reduced_z)], family = binomial)
vif(vif_model4)
alias(vif_model4)

vif_model5 <- glm(ext ~ ., data = mod_extper_dnr_z[, c("ext", covars_numeric_reduced_z)], family = binomial)
vif(vif_model5)
alias(vif_model5)

vif_model6 <- glm(per ~ ., data = mod_perext_dnr_z[, c("per", covars_numeric_reduced_z)], family = binomial)
vif(vif_model6)
alias(vif_model6)


# VIF Thinned Covariate Sets
### Keep relatively loose, keep when VIF < 10 to not thin data too much prior
# to first model ranking step

factor_covars_reduced <- c("atlas_block")

land_covars_reduced <- c("developed_lower_base",
                         "forest_total_base",
                         "grassland_base", "pasture_crop_base",
                         
                         "grassland_diff", "pasture_crop_diff")
land_covars_reduced_z <- paste0(land_covars_reduced, "_z")

climate_covars_reduced <- c("tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff")
climate_covars_reduced_z <- paste0(climate_covars_reduced, "_z")

stable_covars_reduced # "sr_diff", "pa_percent"
stable_covars_reduced_z <- paste0(stable_covars_reduced, "_z")

covars_numeric_reduced <- c(land_covars_reduced, climate_covars_reduced, stable_covars_reduced)
covars_numeric_reduced_z <- paste0(covars_numeric_reduced, "_z")



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

effort_covar <- "sr_diff_z"
pa_covar <- "pa_percent_z"


### Create grid data directory of blocks x response subsets for fx lookup
partition_grid <- expand.grid(
  response = c("col", "ext", "per"),
  blocks = c("DNR", "RLL"),
  partition = c("climate", "land"),
  stringsAsFactors = FALSE
)

data_dir <- list(
  RLL_col = mod_colabs_rll_z,
  RLL_ext = mod_extper_rll_z,
  RLL_per = mod_perext_rll_z,
  
  DNR_col = mod_colabs_dnr_z,
  DNR_ext = mod_extper_dnr_z,
  DNR_per = mod_perext_dnr_z
)

covar_dir <- list(
  climate = climate_covars_reduced_z,
  land    = land_covars_reduced_z
)



### --- STEP 1: COVARIATE PARTITIONED MODELS --- ###

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




### --- STEP 2: GLOBAL MODELS --- ###

### Combine climate, land covars lists from partitioned models into single covar list/pool for 
# new additive mod selection process; find new ref model to carry into interaction step (2B).

effort_covar <- "sr_diff_z"
pa_covar <- "pa_percent_z"


merged_ref_covariates <- lapply(
  merged_ref_covariates,
  function(covs) unique(c(covs, pa_covar)) # add in PA to pool in non-forced manner
)



BuildGlobalRHS <- function(covariates,
                           forced = effort_covar,
                           pa = pa_covar) {
  
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
    forced     = effort_covar,
    pa         = pa_covar
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




global_models <- lapply(names(merged_ref_covariates), function(key) {
  
  FitGlobalModels(
    response   = strsplit(key, "_")[[1]][2],
    data       = data_dir[[key]],
    covariates = merged_ref_covariates[[key]],
    family     = binomial
  )
})

names(global_models) <- names(merged_ref_covariates)

top_global_models <- lapply(global_models, ExtractTopModels, delta = 2)
reference_global_models <- lapply(top_global_models, ExtractReferenceModel)







### --- STEP 3: MODELING --- ###

### Use reference model for each blocks x response subsets to obtain effect sizes, etc.

# Model data (responses: col, abs, per)
data_dir <- list(
  RLL_col = mod_colabs_rll_z,
  RLL_ext = mod_extper_rll_z,
  RLL_per = mod_perext_rll_z,
  
  DNR_col = mod_colabs_dnr_z,
  DNR_ext = mod_extper_dnr_z,
  DNR_per = mod_perext_dnr_z
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
  
  # nm format: "DNR_col", "RLL_ext", etc.
  response <- strsplit(nm, "_")[[1]][2]
  
  glm(
    formula = ExtractRefFormula(reference_global_models[[nm]], response),
    data    = data_dir[[nm]],
    family  = binomial
  )
})

names(global_glm_models) <- names(reference_global_models)

global_glm_summaries <- lapply(global_glm_models, summary)