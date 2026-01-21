#############
### SETUP ###
#############

## --- LOAD PACKAGES --- ##

# Core
library(dplyr)
library(tidyr)
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
library(arm)

# Visualization
library(ggplot2)
library(viridis)
library(gt)
library(webshot2)


### --- DATAFRAMES --- ###

# Carry-over DataFrames #
spp_zf_rll <- read.csv("data/summaries/spp_zf_rll.csv") 
covars_raw_rll <- read.csv("data/summaries/covars_raw_all.csv")

wibba_summary_rll <- read.csv("data/summaries/wibba_summary_rll.csv") # df
blocks_rll <- wibba_summary_rll$atlas_block # vector

blocks_dnr <- read_xlsx("data/summaries/CompBlocks_DNR2023.xlsx") # df
blocks_dnr <- blocks_dnr$atlas_block # vector


# SPECIES RICHNESS / EFFORT PROXY 
covars_raw_rll <- covars_raw_rll %>%
  mutate(sr_Diff = sr_Atlas2 - sr_Atlas1,
         grass_pasture_crop_base = grassland_base + pasture_crop_base,
         grass_pasture_crop_diff = grassland_diff + pasture_crop_diff)


# Covariate Sets #
factor_covars_all <- c("atlas_block", "common_name", "alpha_code", "transition_state")

stable_covars_all <- c("lon", "lat", "sr_Diff", "pa_percent")

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

land_covars_base <- c("water_open_base", "barren_land_base", "shrub_scrub_base", "grassland_base",
                      "developed_open_base", "developed_low_base", "developed_med_base", "developed_high_base", 
                      "developed_lower_base", "developed_upper_base", "developed_total_base", 
                      "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                      "forest_total_base", "pasture_base", "cropland_base", "pasture_crop_base", 
                      "wetlands_woody_base", "wetlands_herb_base", "wetlands_total_base")

land_covars_diff <- c("water_open_diff", "barren_land_diff", "shrub_scrub_diff", "grassland_diff",
                      "developed_open_diff", "developed_low_diff", "developed_med_diff", "developed_high_diff", 
                      "developed_lower_diff", "developed_upper_diff", "developed_total_diff", 
                      "forest_deciduous_diff", "forest_evergreen_diff", "forest_mixed_diff", 
                      "forest_total_diff", "pasture_diff", "cropland_diff", "pasture_crop_diff", 
                      "wetlands_woody_diff", "wetlands_herb_diff", "wetlands_total_diff")

climate_covars_all <- c("tmax_38yr", "tmin_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff")

climate_covars_base <- c("tmax_38yr", "tmin_38yr", "prcp_38yr")

climate_covars_diff <- c("tmax_diff", "tmin_diff", "prcp_diff")



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


write.csv(mod_data_rll, "outputs/data/mod_data_rll.csv", row.names = FALSE)
write.csv(mod_data_dnr, "outputs/data/mod_data_dnr.csv", row.names = FALSE)



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
land_covars_base
land_covars_diff
climate_covars_all
climate_covars_base
climate_covars_diff


# Species-specific Thinned Covariate Sets
spp_name <- "Eastern Meadowlark"


factor_covars_reduced <- c("atlas_block", "transition_state")

stable_covars_all <- c("sr_Diff", "pa_percent")

land_covars_reduced <- c("developed_total_base", 
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         
                         "wetlands_woody_base", "wetlands_herb_base",
                         
                         "developed_total_diff", "forest_evergreen_diff", 
                         "forest_total_diff", "wetlands_total_diff")

covars_numeric_reduced <- c(stable_covars_all, land_covars_reduced, climate_covars_all)


# State-specific Covariate Sets
### Separate data into bins where A2 detection is either 1 or 0; not necessarily
# in precise probabilities between colo, pers, abs, ext, but more what promotes
# 'new' v. 'continued' colonization, ie. what promotes det = 1, as opposed to 
# what promotes det = 0. 
### Scaling of covars w/in subsets for relevant normalized values

spp_name <- "Eastern Meadowlark"

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


corrplot::corrplot(M1, method = "color", tl.cex = 0.7, number.cex = 0.6) # visualize correlation plots
corrplot::corrplot(M2, method = "color", tl.cex = 0.7, number.cex = 0.6)
corrplot::corrplot(M3, method = "color", tl.cex = 0.7, number.cex = 0.6)

corrplot::corrplot(M4, method = "color", tl.cex = 0.7, number.cex = 0.6)
corrplot::corrplot(M5, method = "color", tl.cex = 0.7, number.cex = 0.6)
corrplot::corrplot(M6, method = "color", tl.cex = 0.7, number.cex = 0.6)


high_corr1 <- which(abs(M1) > 0.7 & abs(M1) < 1, arr.ind = TRUE) # list out high correlation pairs
apply(high_corr1, 1, function(i) cat(rownames(M1)[i[1]], "-", colnames(M1)[i[2]], "r =", M1[i[1],i[2]], "\n"))

high_corr2 <- which(abs(M2) > 0.7 & abs(M2) < 1, arr.ind = TRUE)
apply(high_corr2, 1, function(i) cat(rownames(M2)[i[1]], "-", colnames(M2)[i[2]], "r =", M2[i[1],i[2]], "\n"))

high_corr3 <- which(abs(M3) > 0.7 & abs(M3) < 1, arr.ind = TRUE)
apply(high_corr3, 1, function(i) cat(rownames(M3)[i[1]], "-", colnames(M3)[i[2]], "r =", M3[i[1],i[2]], "\n"))

high_corr4 <- which(abs(M4) > 0.7 & abs(M4) < 1, arr.ind = TRUE)
apply(high_corr4, 1, function(i) cat(rownames(M4)[i[1]], "-", colnames(M4)[i[2]], "r =", M4[i[1],i[2]], "\n"))

high_corr5 <- which(abs(M5) > 0.7 & abs(M5) < 1, arr.ind = TRUE)
apply(high_corr5, 1, function(i) cat(rownames(M5)[i[1]], "-", colnames(M5)[i[2]], "r =", M5[i[1],i[2]], "\n"))

high_corr6 <- which(abs(M6) > 0.7 & abs(M6) < 1, arr.ind = TRUE)
apply(high_corr6, 1, function(i) cat(rownames(M6)[i[1]], "-", colnames(M6)[i[2]], "r =", M6[i[1],i[2]], "\n"))


# Correlation Thinned Covariate Sets
### Removed one covar from highly correlated pairs when |r| > 0.7
factor_covars_reduced

land_covars_reduced <- c("developed_total_base", 
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_woody_base", "wetlands_herb_base",
                         
                         "forest_total_diff", "wetlands_total_diff")

climate_covars_reduced <- c("tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff")

stable_covars_reduced <- c("sr_Diff", "pa_percent")

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
                         
                         "forest_total_diff", "wetlands_total_diff")
land_covars_reduced_z <- paste0(land_covars_reduced, "_z")

climate_covars_reduced <- c("tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff")
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

# Model data (responses: col, abs, per)
mod_colabs_rll_z
mod_extper_rll_z
mod_perext_rll_z

mod_colabs_dnr_z
mod_extper_dnr_z
mod_perext_dnr_z

# Covariate IDs/vectors
land_covars_reduced_z
climate_covars_reduced_z

effort_covar <- "sr_Diff_z"
pa_covar <- "pa_percent_z"
required_covars <- c("sr_Diff_z", "pa_percent_z")


### --- STEP 1: COVARIATE PARTITIONED MODELS --- ###

### Full, additive candidate model sets for climate and land cover covariates
# w/ null included. AICc fitting, ranking; extract covars to pass on from "reference" 
# model (ie. delta = 0) to global model selection steps.


### FUNCTIONS ###

### Helper: generate all additive predictor combos w/in partitioned covariate sets
# required = structural effort control variable, PA% as primary effect of interest
BuildAdditiveRHS <- function(covariates, required) {
  
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
FitAdditiveModels <- function(response, 
                              data, 
                              covariates, 
                              required = required_covars,
                              family = binomial,
                              include_null = TRUE) {
  
  # Predictor structure
  covariates <- unique(covariates)
  
  rhs_terms <- BuildAdditiveRHS(
    covariates = covariates,
    required = required
  )
  
  formulas <- lapply(rhs_terms, function(rhs) {
    as.formula(paste(response, "~", rhs))
  })
  
  # Null structure
  if (include_null) {
    formulas <- c(list(as.formula(paste(response, "~ 1"))), formulas)
    modnames <- c("NULL", rhs_terms)
  } else {
    modnames <- rhs_terms
  }
  
  # Fit models
  models <- lapply(formulas, function(f) {
    glm(f, data = data, family = family)
  })
  
  # Identify models 
  names(models) <- modnames
  
  # AICc ranking
  aictab(
    cand.set = models,
    modnames = modnames,
    sort = TRUE
  )
}


### APPLICATION ###

### Automate function application to partitioned climate and land cover models 
# for each block/response subset simultaneously

### Create data directories for fx lookup
# Grid containing rows of all response x blocks x partition combos
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


# Apply: Index, fit partitions in single loop function
partitioned_add_models <- lapply(seq_len(nrow(partition_grid)), function(i) {
  
  # Extract ith row of partition_grid (unique response x block set x covariate partition combo) 
  row <- partition_grid[i, ]
  
  # Look-up key for named blocks x response combo
  key <- paste(row$blocks, row$response, sep = "_")
  
  # Fit, rank add models w/in each partition
  FitAdditiveModels(
    response     = row$response,
    data         = data_dir[[key]],
    covariates   = c(pa_covar, covar_dir[[row$partition]]),
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
ExtractTopModels <- function(aicc_table, delta = 2, drop_null = TRUE) {
  
  df <- as.data.frame(aicc_table)
  
  df <- df[df$Delta_AICc <= delta, ]
  
  if (drop_null) {
    df <- df[df$Modnames != "NULL",]
  }
  
  df
  
}

# Apply: threshold to all candidate models 
top_part_models <- lapply(partitioned_add_models, ExtractTopModels, delta = 2)


### Helper: Extract [REFERENCE] model [delta = 0]
ExtractReferenceModel <- function(aicc_table) {
  
  df <- as.data.frame(aicc_table)
  
  ref <- df[df$Delta_AICc == 0, ]
  
  if (nrow(ref) != 1) {
    warning("Expected [1] reference model, found ", nrow(ref))
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



### --- STEP 2: GLOBAL MODELS --- ###

### STEP 2A: Global Additive Models ###

### Combine climate, land covars lists from partitioned models into single covar list/pool for 
# new additive mod selection process; find new ref model to carry into interaction step (2B).

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
  
  # Assign names to blocks x response subset combos
  combos <- unique(names_dir[, c("blocks", "response")])
  
  # Merge covariate lists for blocks x response combos
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
merged_ref_covariates <- lapply(
  MergePartitionedCovariates(reference_part_covariates),
  function(covs) unique(c(required_covars, covs))
)


# Apply: add. model build, ref model and covar extraction workflow
global_add_models <- lapply(names(merged_ref_covariates), function(key) {
  
  FitAdditiveModels(
    response     = strsplit(key, "_")[[1]][2],
    data         = data_dir[[key]],
    covariates   = merged_ref_covariates[[key]],
    family       = binomial,
    include_null = TRUE
  )
})

names(global_add_models) <- names(merged_ref_covariates)


top_add_models <- lapply(global_add_models, ExtractTopModels, delta = 2)

reference_add_models <- lapply(top_add_models, ExtractReferenceModel)

reference_add_covariates <- lapply(reference_add_models, function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(character(0))
  ExtractReferenceCovariates(df$Modnames[1])
  
})





### --- STEP 2B: INTERACTIONS --- ###

### Custom list
pa_int_covs <- c(
  "tmax_38yr_z",
  "wetlands_woody_base_z"
)

# Helper: build custom interaction terms (w/ PA)
BuildPAInteractions <- function(covariates,
                                pa_var = "pa_percent_z",
                                pa_int_covs) {
  
  if (!pa_var %in% covariates) return(character(0))
  
  int_covs <- intersect(covariates, pa_int_covs)
  int_covs <- setdiff(int_covs, pa_var)
  
  paste0(pa_var, ":", int_covs)
}



### Helper: Construct global model
BuildInteractionRHS <- function(covariates,
                                required = required_covars,
                                pa_var = "pa_percent_z",
                                pa_int_covs) {
  
  covariates <- unique(c(covariates, required))
  
  # Failsafe: ensure PA present
  if (!pa_var %in% covariates) {
    stop("pa_percent_z missing from interaction covariates")
  }
  
  additive_rhs <- paste(covariates, collapse = " + ")
  
  int_terms <- BuildPAInteractions(
    covariates   = covariates,
    pa_var       = pa_var,
    pa_int_covs  = pa_int_covs
  )
  
  # additive-only model
  rhs_set <- additive_rhs
  
  # additive + each interaction individually
  if (length(int_terms) > 0) {
    rhs_set <- c(
      rhs_set,
      paste(additive_rhs, int_terms, sep = " + ")
    )
  }
  
  unique(rhs_set)
}


### Helper: Fit, Rank global and null models
FitInteractionModels <- function(response,
                                 data,
                                 covariates,
                                 required = required_covars,
                                 pa_int_covs,
                                 pa_var = "pa_percent_z",
                                 family = binomial) {
  
  rhs_terms <- BuildInteractionRHS(
    covariates  = covariates,
    required = required,
    pa_var      = pa_var,
    pa_int_covs = pa_int_covs
  )
  
  formulas <- c(
    list(as.formula(paste(response, "~ 1"))),
    lapply(rhs_terms, function(rhs) {
      as.formula(paste(response, "~", rhs))
    })
  )
  
  models <- lapply(formulas, function(f) {
    glm(f, data = data, family = family)
  })
  
  modnames <- c("NULL", rhs_terms)
  names(models) <- modnames
  
  aictab(
    cand.set = models,
    modnames = modnames,
    sort     = TRUE
  )
}




# Run selection 
global_int_models <- lapply(names(reference_add_covariates), function(key) {
  
  FitInteractionModels(
    response    = strsplit(key, "_")[[1]][2],
    data        = data_dir[[key]],
    covariates  = reference_add_covariates[[key]],
    pa_int_covs = pa_int_covs,
    family      = binomial
  )
})

names(global_int_models) <- names(reference_add_covariates)


top_int_models <- lapply(global_int_models, ExtractTopModels, delta = 2)

reference_int_models <- lapply(top_int_models, ExtractReferenceModel)

reference_int_covariates <- lapply(reference_int_models, function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(character(0))
  
  ExtractReferenceCovariates(df$Modnames[1])
})



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


final_glm_models <- lapply(names(reference_int_models), function(nm) {
  
  # nm format: "DNR_col", "RLL_ext", etc.
  response <- strsplit(nm, "_")[[1]][2]
  
  glm(
    formula = ExtractRefFormula(reference_int_models[[nm]], response),
    data    = data_dir[[nm]],
    family  = binomial
  )
})

names(final_glm_models) <- names(reference_int_models)

final_model_summaries <- lapply(final_glm_models, summary)





### --- VISUALIZE --- ###

### Formatted, exportable AICc tables for top model sets (deltaAICc < 2) for each 
# blocks x response subsets at each step.

FormatAICcTable <- function(df) {
  
  df %>%
    mutate(
      AICc       = round(AICc, 2),
      Delta_AICc = round(Delta_AICc, 2),
      AICcWt     = round(AICcWt, 3)
    ) %>%
    dplyr::select(
      Model = Modnames,
      K,
      AICc,
      Delta_AICc,
      AICcWt
    )
}


ParseModelName <- function(name) {
  
  parts <- strsplit(name, "_")[[1]]
  
  if (length(parts) == 3) {
    # partitioned: response_blocks_partition
    list(
      response  = parts[1],
      blocks    = parts[2],
      partition = parts[3]
    )
  } else if (length(parts) == 2) {
    # global models: blocks_response
    list(
      response  = parts[2],
      blocks    = parts[1],
      partition = NA
    )
  } else {
    stop("Unexpected model name format: ", name)
  }
}


PrepareAICcTables <- function(model_list, step_name) {
  
  lapply(names(model_list), function(nm) {
    
    meta <- ParseModelName(nm)
    
    FormatAICcTable(model_list[[nm]]) %>%
      mutate(
        Step      = step_name,
        Response  = meta$response,
        Blocks    = meta$blocks,
        Partition = meta$partition
      )
  }) |> setNames(names(model_list))
}


part_tables <- PrepareAICcTables(top_part_models, "Partitioned")
add_tables  <- PrepareAICcTables(top_add_models,  "Additive")
int_tables  <- PrepareAICcTables(top_int_models,  "Interaction")



SaveAICcPNG <- function(df, file, title, subtitle) {
  
  gt(df) %>%
    tab_header(
      title    = title,
      subtitle = subtitle
    ) %>%
    cols_align("center", everything()) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels(everything())
    ) %>%
    opt_all_caps() %>%
    opt_table_outline() %>%
    gtsave(file, expand = 15)
}


ExportAICcTables <- function(tables, out_dir) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (nm in names(tables)) {
    
    df <- tables[[nm]]
    
    title <- paste(
      unique(df$Step), "Models —",
      unique(df$Blocks), toupper(unique(df$Response))
    )
    
    subtitle <- if (!all(is.na(df$Partition))) {
      paste("Covariate set:", unique(df$Partition))
    } else {
      "ΔAICc ≤ 2"
    }
    
    SaveAICcPNG(
      df = df %>% dplyr::select(Model, K, AICc, Delta_AICc, AICcWt),
      file = file.path(out_dir, paste0("AICc_", nm, ".png")),
      title = title,
      subtitle = subtitle
    )
  }
}


ExportAICcTables(part_tables, "outputs/tables/aicc_partitioned")
ExportAICcTables(add_tables,  "outputs/tables/aicc_additive")
ExportAICcTables(int_tables,  "outputs/tables/aicc_interaction")