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



##############################  GUILD-BASED COVARIATE SETS ######################


grass_guild_all

forest_guild_all

wetland

riparian

generalist

water






#################################################################################




# MODELING DATA FRAMES #

# Species to Model #
spp_name <- "Ruby-crowned Kinglet"


# RLL blocks modeling df
mod_data_rll <- spp_zf_rll %>%
  filter(common_name == spp_name) %>%
  left_join(covars_raw_rll, by = "atlas_block")

# DNR blocks modeling df
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

# Full covariate sets
factor_covars_all
stable_covars_all
land_covars_all
climate_covars_all



# Species-specific Thinned Covariate Sets

factor_covars_reduced <- c("atlas_block", "transition_state")

stable_covars_all <- c("sr_diff", "pa_percent")

land_covars_reduced <- c("developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_woody_base", "wetlands_herb_base",
                         
                         "developed_total_diff", "forest_total_diff", "wetlands_total_diff")

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


# corrplot::corrplot(M1, method = "color", tl.cex = 0.7, number.cex = 0.6)
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

land_covars_reduced <- c("developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_woody_base", "wetlands_herb_base",
                         
                         "developed_total_diff", "forest_total_diff", "wetlands_total_diff")

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

land_covars_reduced <- c("developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_woody_base", "wetlands_herb_base",
                         
                         "developed_total_diff", "forest_total_diff", "wetlands_total_diff")
land_covars_reduced_z <- paste0(land_covars_reduced, "_z")

climate_covars_reduced <- c("tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff")
climate_covars_reduced_z <- paste0(climate_covars_reduced, "_z")

stable_covars_reduced # "sr_Diff", "pa_percent"
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





### --- Model Comparison --- ###
# McFadden Pseudo-R2 b/c binomial glm

# Refit null models
McFaddenR2 <- function(model) {
  ll_full <- as.numeric(logLik(model))
  null <- glm(
    formula = reformulate("1", response = all.vars(formula(model))[1]),
    data = model$data,
    family = model$family
  )
  
  ll_null <- as.numeric(logLik(null))
  1 - (ll_full / ll_null)
}


glm_r2 <- data.frame(
  model = names(global_glm_models),
  R2_McFadden = vapply(
    global_glm_models,
    McFaddenR2,
    numeric(1)
  )
)


# Response x blocks pairs
glm_r2_comp <- glm_r2 %>%
  separate(model, into = c("source", "response"), sep = "_") %>%
  pivot_wider(
    names_from = source,
    values_from = R2_McFadden
  ) %>%
  mutate(
    delta_R2 = RLL - DNR
  )

# Visualize
ggplot(glm_r2_comp, aes(x = response, y = delta_R2)) +
  geom_col(fill = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    y = expression(Delta~R^2~"(RLL - DNR)"),
    x = "response",
    title = "Difference in Model Fit (McFadden R2)"
  ) +
  theme_bw()


# Likelihood per observation
ll_comp <- data.frame(
  model = names(global_glm_models),
  logLik = sapply(global_glm_models, logLik),
  n = sapply(global_glm_models, nobs)
) %>%
  tidyr::separate(model, into = c("source", "response"), sep = "_") %>%
  dplyr::mutate(
    ll_per_obs = as.numeric(logLik) / n
  ) %>%
  dplyr::select(source, response, ll_per_obs) %>% 
  tidyr::pivot_wider(
    names_from = source,
    values_from = ll_per_obs
  ) %>%
  dplyr::mutate(
    delta_ll = RLL - DNR
  )

# Visualize
ll_comp %>%
  tidyr::pivot_longer(c(DNR, RLL), names_to = "source", values_to = "ll") %>%
  ggplot(aes(source, ll, fill = source)) +
  geom_col() +
  facet_wrap(~ response, scales = "free_y") +
  labs(
    y = "Log-likelihood per observation",
    x = NULL
  ) +
  theme_minimal()

















########################## FOSSILS #############################################

## using AIcmodavg


### Helper: Build main effects formulas
BuildGlobalRHS <- function(selected_covariates,
                           required = effort_covar,
                           optional = pa_covar) {
  
  selected_covariates <- unique(c(selected_covariates, optional))
  required <- unique(required)
  optional <- setdiff(selected_covariates, required)
  
  rhs_list <- list()
  rhs_list[[1]] <- paste(required, collapse = " + ")
  
  if (length(optional) > 0) {
    rhs_optional <- unlist(
      lapply(seq_along(optional), function(k) {
        combn(optional, k , FUN = function(x) {
          paste(c(required, x), collapse = " + ")
        })
      }),
      use.names = FALSE
    )
    rhs_list <- c(rhs_list, rhs_optional)
  }
  
  unique(rhs_list)
}



FitGlobalModels <- function(response,
                            data,
                            covariates,
                            required = effort_covar,
                            optional = pa_covar,
                            family = binomial,
                            include_null = TRUE) {
  
  rhs_terms <- BuildGlobalRHS(
    selected_covariates = covariates,
    required = required,
    optional = optional
  )
  
  if (length(rhs_terms) == 0) {
    rhs_terms <- paste(required, collapse = " + ")
  }
  
  formulas <- lapply(rhs_terms, function(rhs) {
    as.formula(paste(response, "~", rhs))
  })
  
  modnames <- rhs_terms
  
  if (include_null) {
    formulas <- c(list(glm(as.formula(paste(response, "~ 1")), data = data, family = family)), formulas)
    modnames <- c("NULL", rhs_terms)
  } else {
    modnames <- rhs_terms
  }
  
  models <- lapply(formulas, function(f) {
    glm(f, data = data, family = family)
  }) 
  
  names(models) <- modnames
  
  aictab(
    cand.set = models,
    modnames = modnames,
    sort = TRUE
  )
}



# Apply: to all blocks x response subsets
global_models <- lapply(names(merged_ref_covariates), function(key) {
  
  FitGlobalModels(
    response = strsplit(key, "_")[[1]][2],
    data = data_dir[[key]],
    covariates = merged_ref_covariates[[key]],
    required = effort_covar,
    optional = pa_covar,
    family = binomial
  )
})


names(global_models) <- names(merged_ref_covariates)

top_global_models <- lapply(global_models, ExtractTopModels, delta = 2)
reference_global_models <- lapply(top_global_models, ExtractReferenceModel)



##################################








# WIP update to main effects-only global models (BROKEN)

# Function: Build formula for main effects models
BuildGlobalRHS <- function(covariates,
                           required = effort_covar,
                           optional = pa_covar) {
  
  covariates <- unique(c(covariates, optional))
  required <- unique(required)
  optional <- setdiff(covariates, required)
  
  rhs_list <- list()
  
  rhs_list[[1]] <- paste(required, collapse = " + ")
  
  if (length(optional) > 0) {
    rhs_optional <- unlist(
      lapply(seq_along(optional), function(k) {
        combn(optional, k , FUN = function(x) {
          paste(c(required, x), collapse = " + ")
        })
      }),
      use.names = FALSE
    )
    rhs_list <- c(rhs_list, rhs_optional)
  }
  unique(rhs_list)
}



FitGlobalModels <- function(response,
                            data,
                            covariates,
                            required = effort_covar,
                            optional = pa_covar,
                            family = binomial,
                            include_null = TRUE) {
  
  rhs_terms <- BuildGlobalRHS(
    covariates = covariates,
    required = required,
    optional = optional
  )
  
  if (length(rhs_terms) == 0) {
    rhs_terms <- paste(required, collapse = " + ")
  }
  
  formulas <- lapply(rhs_terms, function(rhs) {
    as.formula(paste(response, "~", rhs))
  })
  
  modnames <- rhs_terms
  
  if (include_null) {
    formulas <- c(list(as.formula(paste(response, "~ 1"))), formulas)
    modnames <- c("NULL", rhs_terms)
  }
  
  models <- lapply(formulas, function(f) {
    glm(f, data = data, family = family)
  })
  
  names(models) <- modnames
  
  aictab(
    cand.set = models,
    modnames = modnames,
    sort = TRUE
  )
}



# Apply: to all blocks x response subsets
global_models <- lapply(names(merged_ref_covariates), function(key) {
  
  FitGlobalModels(
    response = strsplit(key, "_")[[1]][2],
    data = data_dir[[key]],
    covariates = merged_ref_covariates[[key]],
    required = effort_covar,
    optional = pa_covar,
    family = binomial
  )
})

names(global_models) <- names(merged_ref_covariates)

top_global_models <- lapply(global_models, ExtractTopModels, delta = 2)

reference_global_models <- lapply(top_global_models, ExtractReferenceModel)













### Model workflow incl interaction terms (1.29.26)
pa_int_covs <- c("tmax_38yr_z",
                 "prcp_38yr_z",
                 "tmax_diff_z",
                 "tmin_diff_z",
                 "prcp_diff_z",
                 "developed_total_base_z",
                 "forest_deciduous_base_z", 
                 "forest_evergreen_base_z", 
                 "forest_mixed_base_z",
                 "wetlands_woody_base_z",
                 "forest_total_diff_z"
)



BuildGlobalIntRHS <- function(selected_covariates,
                           effort = effort_covar,
                           pa = pa_covar,
                           pa_int_covs) {
  
  # Partition-selected covariates for this block x response
  selected_covariates <- unique(selected_covariates)
  
  # Interaction-eligible covariates that were ALSO selected upstream
  int_covs <- intersect(selected_covariates, pa_int_covs)
  
  rhs_list <- character(0)
  
  # ---- 1. Additive (no PA unless selected upstream) ----
  base_main <- unique(c(effort, setdiff(selected_covariates, pa)))
  rhs_list <- c(rhs_list, paste(base_main, collapse = " + "))
  
  # ---- 2. Additive with PA only if PA was selected ----
  if (pa %in% selected_covariates) {
    rhs_list <- c(
      rhs_list,
      paste(unique(c(base_main, pa)), collapse = " + ")
    )
  }
  
  # ---- 3. Interaction models (enforce hierarchy) ----
  if (length(int_covs) > 0) {
    for (x in int_covs) {
      rhs_list <- c(
        rhs_list,
        paste(
          unique(c(base_main, pa, x)),
          collapse = " + "
        ) %>%
          paste(paste0(pa, ":", x), sep = " + ")
      )
    }
  }
  
  unique(rhs_list)
}



FitGlobalModels <- function(response,
                            data,
                            covariates,
                            effort = effort_covar,
                            pa_int_covs,
                            family = binomial,
                            include_null = TRUE) {
  
  rhs_terms <- BuildGlobalRHS(
    selected_covariates = covariates,
    effort = effort,
    pa = pa_covar,
    pa_int_covs = pa_int_covs
  )
  
  formulas <- lapply(rhs_terms, function(rhs) {
    as.formula(paste(response, "~", rhs))
  })
  
  if (include_null) {
    formulas <- c(list(as.formula(paste(response, "~ 1"))), formulas)
    modnames <- c("NULL", rhs_terms)
  } else {
    modnames <- rhs_terms
  }
  
  models <- lapply(formulas, function(f) {
    glm(f, data = data, family = family)
  })
  
  names(models) <- modnames
  
  aictab(
    cand.set = models,
    modnames = modnames,
    sort = TRUE
  )
}


# Apply: to all blocks x response subsets
global_models <- lapply(names(merged_ref_covariates), function(key) {
  
  FitGlobalModels(
    response = strsplit(key, "_")[[1]][2],
    data = data_dir[[key]],
    covariates = merged_ref_covariates[[key]],
    pa_int_covs = pa_int_covs,
    family = binomial
  )
})

names(global_models) <- names(merged_ref_covariates)

top_global_models <- lapply(global_models, ExtractTopModels, delta = 2)

reference_global_models <- lapply(top_global_models, ExtractReferenceModel)

























# Quick diagnostics: number of top models in each set and head()
lapply(Top_global_models, function(df) {
  list(n_models = nrow(df), head = if (nrow(df)>0) head(df, 6) else df)
})



# Identify uninformative parameters #
### Compare K sets among all models w/ delta < 2 

# Helper: Get Reference Model (lowest AICc, smallest K as tie-breaker)
GetReferenceModel <- function(df) {
  df[order(df$AICc, df$K), ][1, ]
}

# Apply
ReferenceModels <- lapply(Top_models_list, GetReferenceModel)
ReferenceModels

# Helper: flag uninformative models 
FlagUninformativeModels <- function(df, aicc_tol = 2) {
  # Reference model: lowest AICc, tie-break with smallest K
  ref <- GetReferenceModel(df)
  
  df$extra_K <- df$K - ref$K
  df$aicc_gain <- ref$AICc - df$AICc
  
  df$uninformative_model <- with(
    df,
    extra_K > 0 & aicc_gain < aicc_tol
  )
  
  df
}


# Helper: extract parameters from automated model construction
ExtractTerms <- function(formula_string) {
  rhs <- gsub(".*~", "", formula_string)
  terms <- trimws(unlist(strsplit(rhs, "\\+")))
  terms <- terms[terms != "(Intercept)" & terms != "1"]
  terms
}


# Helper: build sumamry table comparing candidate moedls to refrence model
GetTermSupport <- function(df, always_keep = c("pa_percent_z", "sr_Diff_z")) {
  
  df <- FlagUninformativeModels(df)
  ref_model <- GetReferenceModel(df)
  ref_terms <- ExtractTerms(ref_model$formula)
  
  term_df <- do.call(
    rbind,
    lapply(seq_len(nrow(df)), function(i) {
      terms <- ExtractTerms(df$formula[i])
      data.frame(
        term = terms,
        model_id = df$model_id[i],
        K = df$K[i],
        AICc = df$AICc[i],
        uninformative_model = df$uninformative_model[i],
        in_reference = terms %in% ref_terms,
        stringsAsFactors = FALSE
      )
    })
  )
  
  # Aggregate counts per term
  agg <- aggregate(
    cbind(n_models = model_id, n_uninf = uninformative_model) ~ term,
    data = term_df,
    FUN = function(x) if(is.logical(x)) sum(x) else length(x)
  )
  
  # Include reference info (does this term appear in reference model?)
  agg$in_reference <- agg$term %in% ref_terms
  
  # Classify as uninfomrative
  agg$category <- "Supported"
  agg$category[!agg$in_reference & agg$n_models == agg$n_uninf] <- "Uninformative"
  
  # force keep PA and SR terms as Supported
  agg$category[agg$term %in% always_keep] <- "Supported"
  
  # Sort for readability: reference terms first
  agg <- agg[order(!agg$in_reference, -agg$n_models), ]
  
  agg
}

# Apply
Supported_terms_clean <- lapply(Top_models_list, GetTermSupport)
lapply(Supported_terms_clean, head, 10)

# Helper: build summary table for best candidate moedl parameters

GetBestCandidateTerms <- function(supported_terms_list, keep_terms = c("pa_percent_z", "sr_Diff_z")) {
  
  lapply(supported_terms_list, function(term_df) {
    # Always keep specified terms
    term_df$category <- as.character(term_df$category)
    term_df$category[term_df$term %in% keep_terms] <- "Supported"
    
    # Keep only supported terms
    best_terms <- term_df[term_df$category == "Supported", , drop = FALSE]
    
    # Flag interactions separately
    best_terms$interaction <- grepl(":", best_terms$term)
    
    # Arrange: interactions last, alphabetically
    best_terms <- best_terms[order(best_terms$interaction, best_terms$term), ]
    
    # Return relevant columns
    best_terms[, c("term", "n_models", "interaction", "category")]
  })
}

# Apply
BestCandidateTerms <- GetBestCandidateTerms(Supported_terms_clean)
lapply(BestCandidateTerms, head, 10)



### --- MODEL EXAMINATION, DIAGNOSTICS --- ###
### Use reference models for each block set, response combo 
# Can't use ggfortify::autoplot for non-Gaussian; insetad, arm::binnedplot() for residual plots,
# performance::check_outliers 

# Run, diagnose reference models
# RLL, ColABs
rll_colabs_ref <- glm(col_abs ~
                        cropland_base_z + grassland_base_z + pa_percent_z + pasture_base_z + 
                        pasture_crop_diff_z + prcp_38yr_z + shrub_scrub_base_z + sr_Diff_z + 
                        tmax_38yr_z + tmin_diff_z + cropland_base_z:pa_percent_z + 
                        pa_percent_z:pasture_base_z + pa_percent_z:tmax_38yr_z,
                      data = mod_colabs_rll_z,
                      family = "binomial"
)

summary(rll_colabs_ref)
binnedplot(
  x = fitted(rll_colabs_ref),
  y = residuals(rll_colabs_ref, type = "response")
)
check_outliers(rll_colabs_ref)


# RLL, ExtPer
rll_extper_ref <- glm(ext_per ~
                        cropland_base_z + grassland_base_z + pa_percent_z + 
                        pasture_base_z + prcp_38yr_z + shrub_scrub_base_z + 
                        sr_Diff_z + tmax_38yr_z + tmax_diff_z + 
                        grassland_base_z:pa_percent_z + pa_percent_z:sr_Diff_z + 
                        pa_percent_z:tmax_38yr_z,
                      data = mod_extper_rll_z,
                      family = "binomial"
)

summary(rll_extper_ref)
binnedplot(
  x = fitted(rll_extper_ref),
  y = residuals(rll_extper_ref, type = "response")
)
check_outliers(rll_extper_ref)






# DNR, ColAbs
dnr_colabs_ref <- glm(col_abs ~
                        pa_percent_z +
                        cropland_base_z + developed_total_base_z + forest_total_base_z + 
                        pasture_base_z + pasture_crop_diff_z + prcp_38yr_z + sr_Diff_z + tmin_diff_z,
                      data = mod_colabs_dnr_z,
                      family = "binomial"                 
)

summary(dnr_colabs_ref)
binnedplot(
  x = fitted(dnr_colabs_ref),
  y = residuals(dnr_colabs_ref, type = "response")
)
check_outliers(dnr_colabs_ref)


# DNR, ExtPer
dnr_extper_ref <- glm(ext_per ~
                        cropland_base_z + grassland_base_z + pa_percent_z + 
                        pasture_base_z + shrub_scrub_base_z + shrub_scrub_diff_z + 
                        sr_Diff_z + tmax_38yr_z + tmax_diff_z + cropland_base_z:pa_percent_z + 
                        grassland_base_z:pa_percent_z + pa_percent_z:pasture_base_z + 
                        pa_percent_z:shrub_scrub_base_z + pa_percent_z:sr_Diff_z,
                      data = mod_extper_dnr_z,
                      family = "binomial"                      
)

summary(dnr_extper_ref)
binnedplot(
  x = fitted(dnr_extper_ref),
  y = residuals(dnr_extper_ref, type = "response")
)
check_outliers(dnr_extper_ref)


# Visualize PA Effects #
# Consolidate reference models
ref_models <- list(
  "RLL Colonization" = rll_colabs_ref,
  "RLL Extinction"   = rll_extper_ref,
  "DNR Colonization" = dnr_colabs_ref,
  "DNR Extinction"   = dnr_extper_ref
)

# Extract PA coefficient value, se
pa_effects <- lapply(names(ref_models), function(mod_name){
  tidy_mod <- broom::tidy(ref_models[[mod_name]])
  
  # PA main efffects only
  pa_row <- tidy_mod %>% filter(term == "pa_percent_z") %>%
    mutate(Model = mod_name)
  
  return(pa_row)
}) %>% bind_rows()


# Visualize: Coefficient plot (PA main effects)
ggplot(pa_effects, aes(x = Model, y = estimate, ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error)) +
  geom_pointrange(color = "orange", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Model",
    y = "Standardized Effect of Protected Area",
    caption = "Figure 1. Standardzied main effect of protected area (PA) on apparent colonization and extinction 
    for the Eastern Meadowlark among two survey unit subsets. Points represent effect estimates for each model, 
    with horizantal lines representing 95% confidence intervals."
  ) +
  theme_minimal(base_size = 13) +
  coord_flip() + # flip for horizontal readability
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    plot.caption = element_text(hjust = 0.5, size = 11, margin = margin(t = 15))
  ) 



# Assumption Checking #
# Helper: extract binned residuals to create a single, faceted binned reidual plot for all model results
# Helper function: compute binned residuals
ComputeBinnedResiduals <- function(model, label, bins = 20) {
  df <- tibble(
    fitted = fitted(model),
    resid  = residuals(model, type = "response")
  )
  
  # create equal-sized bins of fitted probabilities
  df <- df %>%
    mutate(bin = cut(fitted,
                     breaks = quantile(fitted, probs = seq(0, 1, length.out = bins + 1)),
                     include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      fitted_mean = mean(fitted, na.rm = TRUE),
      resid_mean  = mean(resid, na.rm = TRUE),
      resid_se    = sd(resid, na.rm = TRUE)/sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(model = label)
  
  return(df)
}


binned_all <- bind_rows(
  ComputeBinnedResiduals(rll_colabs_ref, "RLL – Colonization–Absence"),
  ComputeBinnedResiduals(rll_extper_ref, "RLL – Extinction–Persistence"),
  ComputeBinnedResiduals(dnr_colabs_ref, "DNR – Colonization–Absence"),
  ComputeBinnedResiduals(dnr_extper_ref, "DNR – Extinction–Persistence")
)


# construct faceted plot
ggplot(binned_all, aes(x = fitted_mean, y = resid_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = resid_mean - resid_se,
                    ymax = resid_mean + resid_se), width = 0.02) +
  facet_wrap(~ model, ncol = 2) +
  labs(
    x = "Mean fitted probability",
    y = "Mean binned residual",
    caption = "Figure S2. Binned residual diagnostic plots for each block set x response logistic model combination.
    Each point represents the mean residual within a group of fitted probabilities, with vertical error bars representing
    a +/- 1 standard error of the mean residual."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 15))
  ) 


### --- OTHER SUPPLEMENTAL PLOTS --- ###

# Species Richness #
### Histogram of SR diff between RLL and NDR comp block sets

# Data
srdiff_rll <- covars_raw_rll %>% # RLL blocks
  filter(atlas_block %in% blocks_rll) %>%
  dplyr::select(sr_Diff) %>%
  mutate(source = "rll")

srdiff_dnr <- covars_raw_rll %>% # DNR blocks
  filter(atlas_block %in% blocks_dnr) %>%
  dplyr::select(sr_Diff) %>%
  mutate(source = "dnr")

srdiff_hist_data <- bind_rows(srdiff_rll, srdiff_dnr) # combine

# Plot
ggplot(srdiff_hist_data, aes(x = sr_Diff, fill = source)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(
    values = c("rll" = "orange", "dnr" = "steelblue"),
    labels = c("DNR", "RLL")
  ) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Species Richness Difference (Atlas 2 - Atlas 1)",
    y = "Count",
    fill = "Block Set",
    caption = "Figure S1. Distribution of difference in species richness between Wisconsin Breeding Bird Atlas 1 (1995-2000) and 2 (2015-2019) between two
    different survey block subsets (All Blocks, N = 7056; RLL, n = 2535; DNR, n = 858). Overall, blocks in the second Atlas were surveyed more comprehensively than 
    in Atlas 1, resulting in higher-per-block species richness in Atlas 2 generally. The minimally filtered RLL block set retains a significant number of 
    survey units with high species richness differences compared to the heavily filtered DNR block set, which largely controlled for these coverage discrepancies."
  ) +
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 15))
  )