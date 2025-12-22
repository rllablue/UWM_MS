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


### --- FLEXIBLE SPECS --- ###

spp_name <- "Red-bellied Woodpecker"


### --- DATAFRAMES --- ###

# Carry-over DataFrames #
spp_zf_rll <- read.csv("data/summaries/spp_zf_rll.csv") 
covars_raw_rll <- read.csv("outputs/data/covars_raw_all.csv")

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
# RLL modeling df
mod_data_rll <- spp_zf_rll %>%
  filter(common_name == spp_name) %>%
  left_join(covars_raw_rll, by = "atlas_block")

write.csv(mod_data_rll, "outputs/data/mod_data_rll.csv", row.names = FALSE)

# DNR modeling df
mod_data_dnr <- spp_zf_rll %>%
  filter(atlas_block %in% blocks_dnr, 
         common_name == spp_name) %>%  
  left_join(covars_raw_rll, by = "atlas_block")

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
spp_name <- "Red-bellied Woodpecker"


factor_covars_reduced <- c("atlas_block", "transition_state")

stable_covars_all <- c("lon", "lat", "sr_Diff", "pa_percent")

land_covars_reduced <- c("shrub_scrub_base", "developed_total_base",
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

spp_name <- "Red-bellied Woodpecker"

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
  mutate(col_abs = ifelse(transition_state == "Colonization", 1, 0)) %>% # binomial response variable
  dplyr::select(all_of(factor_covars_reduced), # columns to keep
                ends_with("_z"),
                col_abs)

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
  mutate(ext_per = ifelse(transition_state == "Extinction", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                ext_per)

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
  mutate(per_ext = ifelse(transition_state == "Persistence", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                per_ext)


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
  mutate(col_abs = ifelse(transition_state == "Colonization", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                col_abs)

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
  mutate(ext_per = ifelse(transition_state == "Extinction", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                ext_per)

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
  mutate(per_ext = ifelse(transition_state == "Persistence", 1, 0)) %>%
  dplyr::select(all_of(factor_covars_reduced),
                ends_with("_z"),
                per_ext)


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

land_covars_reduced <- c("shrub_scrub_base", "developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_woody_base", "wetlands_herb_base",
                         
                         "developed_total_diff", "forest_total_diff", "wetlands_total_diff")

climate_covars_reduced <- c("tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff")

stable_covars_reduced <- c("sr_Diff", "pa_percent")

covars_numeric_reduced <- c(land_covars_reduced, climate_covars_reduced, stable_covars_reduced)
covars_numeric_reduced_z <- paste0(covars_numeric_reduced, "_z")



# VIF
vif_model1 <- glm(col_abs ~ ., data = mod_colabs_rll_z[, c("col_abs", covars_numeric_reduced_z)], family = binomial)
vif(vif_model1)
alias(vif_model1)

vif_model2 <- glm(ext_per ~ ., data = mod_extper_rll_z[, c("ext_per", covars_numeric_reduced_z)], family = binomial)
vif(vif_model2)
alias(vif_model2)

vif_model3 <- glm(per_ext ~ ., data = mod_perext_rll_z[, c("per_ext", covars_numeric_reduced_z)], family = binomial)
vif(vif_model3)
alias(vif_model3)

vif_model4 <- glm(col_abs ~ ., data = mod_colabs_dnr_z[, c("col_abs", covars_numeric_reduced_z)], family = binomial)
vif(vif_model4)
alias(vif_model4)

vif_model5 <- glm(ext_per ~ ., data = mod_extper_dnr_z[, c("ext_per", covars_numeric_reduced_z)], family = binomial)
vif(vif_model5)
alias(vif_model5)

vif_model6 <- glm(per_ext ~ ., data = mod_perext_dnr_z[, c("per_ext", covars_numeric_reduced_z)], family = binomial)
vif(vif_model6)
alias(vif_model6)


# VIF Thinned Covariate Sets
# Run 2+ x to check incoming and reduced/outgoing covariate set
### Keep relatively loose, keep when VIF < 10 to not thin data too much prior
# to first model ranking step

factor_covars_reduced <- c("atlas_block")

land_covars_reduced <- c("shrub_scrub_base", "developed_total_base",
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

# 2) Interaction Terms
# Generate plausible 2-way interactions w/ protected area

# 3) Global AICc Model Selection
# Models w/ best land, climate covariates + interaction terms w/ PA


### DATA ###
# Relevant carry-in data from above workflow

# Model data
mod_colabs_rll_z
mod_extper_rll_z
mod_perext_rll_z

mod_colabs_dnr_z
mod_extper_dnr_z
mod_perext_dnr_z

# Covariate IDs/vectors
land_covars_reduced
land_covars_reduced_z
climate_covars_reduced
climate_covars_reduced_z



### --- STEP 1: COVARIATE PARTITIONED MODELS --- ###

### Full, additive candidate model set for climate and land cover covariates
# w/ null model and no interaction terms included

### FUNCTIONS ###
# Helper: generate all additive combinations w/in partitioned covariate sets
GenerateAdditiveSets <- function(covariates) {
  unlist(
    lapply(seq_along(covariates), function(k) {
      combn(covariates, k, FUN = function(x) paste(x, collapse = "+" ))
    }),
    use.names = FALSE
  )
}


# Helper: build, fit, rank partitioned candidate models
FitPartitionedModels <- function(response, data, covariates, family = binomial) {
  
  # Model formula construction
  predictors <- GenerateAdditiveSets(covariates)
  
  formulas <- c(
    list(as.formula(paste(response, "~ 1"))),
    lapply(predictors, function(rhs) {
      as.formula(paste(response, "~", rhs))
    })
  )
  
  # Fit models
  models <- lapply(formulas, function(f) {
    glm(f, data = data, family = family)
  })
  
  # Identify models 
  modnames <- c("NULL", predictors)
  names(models) <- modnames
  
  # AICc ranking
  aictab(
    cand.set = models,
    modnames = modnames,
    sort = TRUE
  )
}


### APPLICATION ###

# Models by unique block x response subsets #
# Climate covariates partition
climate_col_rll <- FitPartitionedModels(
  response = "col_abs",
  data = mod_colabs_rll_z,
  covariates = climate_covars_reduced_z,
  family = binomial
)

climate_col_dnr <- FitPartitionedModels(
  response = "col_abs",
  data = mod_colabs_dnr_z,
  covariates = climate_covars_reduced_z,
  family = binomial
)


climate_ext_rll <- FitPartitionedModels(
  response = "ext_per",
  data = mod_extper_rll_z,
  covariates = climate_covars_reduced_z,
  family = binomial
)

climate_ext_dnr <- FitPartitionedModels(
  response = "ext_per",
  data = mod_extper_dnr_z,
  covariates = climate_covars_reduced_z,
  family = binomial
)


climate_per_rll <- FitPartitionedModels(
  response = "per_ext",
  data = mod_perext_rll_z,
  covariates = climate_covars_reduced_z,
  family = binomial
)

climate_per_dnr <- FitPartitionedModels(
  response = "per_ext",
  data = mod_perext_dnr_z,
  covariates = climate_covars_reduced_z,
  family = binomial
)


# Land cover covariates partition
land_col_rll <- FitPartitionedModels(
  response = "col_abs",
  data = mod_colabs_rll_z,
  covariates = land_covars_reduced_z,
  family = binomial
)

land_col_dnr <- FitPartitionedModels(
  response = "col_abs",
  data = mod_colabs_dnr_z,
  covariates = land_covars_reduced_z,
  family = binomial
)


land_ext_rll <- FitPartitionedModels(
  response = "ext_per",
  data = mod_extper_rll_z,
  covariates = land_covars_reduced_z,
  family = binomial
)

land_ext_dnr <- FitPartitionedModels(
  response = "ext_per",
  data = mod_extper_dnr_z,
  covariates = land_covars_reduced_z,
  family = binomial
)


land_per_rll <- FitPartitionedModels(
  response = "per_ext",
  data = mod_perext_rll_z,
  covariates = land_covars_reduced_z,
  family = binomial
)

land_per_dnr <- FitPartitionedModels(
  response = "per_ext",
  data = mod_perext_dnr_z,
  covariates = land_covars_reduced_z,
  family = binomial
)


# Store results #
partition_results <- list(
  col_rll_climate = climate_col_rll,
  col_rll_land    = land_col_rll,
  col_dnr_climate = climate_col_dnr,
  col_dnr_land    = land_col_dnr,
  
  ext_rll_climate = climate_ext_rll,
  ext_rll_land    = land_ext_rll,
  ext_dnr_climate = climate_ext_dnr,
  ext_dnr_land    = land_ext_dnr,
  
  per_rll_climate = climate_per_rll,
  per_rll_land    = land_per_rll,
  per_dnr_climate = climate_per_dnr,
  per_dnr_land    = land_per_dnr
)

# Check NULL model performance #
climate_col_rll[climate_col_rll$Modnames == "NULL", ]
land_col_rll[land_col_rll$Modnames == "NULL", ]

climate_col_dnr[climate_col_dnr$Modnames == "NULL", ]
land_col_dnr[land_col_dnr$Modnames == "NULL", ]


climate_ext_rll[climate_ext_rll$Modnames == "NULL", ]
land_ext_rll[land_ext_rll$Modnames == "NULL", ]

climate_ext_dnr[climate_ext_dnr$Modnames == "NULL", ]
land_ext_dnr[land_ext_dnr$Modnames == "NULL", ]


climate_per_rll[climate_per_rll$Modnames == "NULL", ]
land_per_rll[land_per_rll$Modnames == "NULL", ]

climate_per_dnr[climate_per_dnr$Modnames == "NULL", ]
land_per_dnr[land_per_dnr$Modnames == "NULL", ]

### RBWO: NULLs unsupported for all block x response subsets


# Evaluate top model partitions #
# Identify supported covariates from each partition (deltaAICc <= 2)

# Col-RLL
top_climate_col_rll <- climate_col_rll %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_climate_col_rll)

top_land_col_rll <- land_col_rll %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_land_col_rll)

# Col-DNR
top_climate_col_dnr <- climate_col_dnr %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_climate_col_dnr)

top_land_col_dnr <- land_col_dnr %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_land_col_dnr)


# Ext-RLL
top_climate_ext_rll <- climate_ext_rll %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_climate_ext_rll)

top_land_ext_rll <- land_ext_rll %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_land_ext_rll)

# Ext-DNR
top_climate_ext_dnr <- climate_ext_dnr %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_climate_ext_dnr)

top_land_ext_dnr <- land_ext_dnr %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_land_ext_dnr)


# Per-RLL
top_climate_per_rll <- climate_per_rll %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_climate_per_rll)

top_land_per_rll <- land_per_rll %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_land_per_rll)

# Per-DNR
top_climate_per_dnr <- climate_per_dnr %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_climate_per_dnr)

top_land_per_dnr <- land_per_dnr %>%
  as.data.frame() %>%
  filter(Delta_AICc <= 2, Modnames != "NULL")
View(top_land_per_dnr)



# Covariate sets per block x response subset #
### Chosen from [REFERENCE] among top models (deltaAICc <= 2)

RLL_col_abs_covs <- c(

)

RLL_ext_per_covs <- c(

)

RLL_per_ext_covs <- c(
  
)


DNR_col_abs_covs <- c(

)

DNR_ext_per_covs <- c(

)

DNR_per_ext_covs <- c(
  
)






### --- STEP 2: INTERACTIONS --- ###

# Custom list
pa_int_covs <- c(

  "sr_Diff_z"
)

# Helper: build interaction terms (w/ PA)
BuildInteractions <- function(response, covars, pa_int_list) {
  
  main_terms <- covars
  
  # only construct interactions for overlap of covars and chosen PA-int covars
  int_covs <- intersect(covars, pa_int_list)
  int_terms <- paste0("pa_percent_z:", int_covs)
  
  rhs <- paste(c(main_terms, int_terms), collapse = " + ")
  
  as.formula(paste(response, "~", rhs))
}








### --- STEP 3: GLOBAL MODELS --- ###

# Build global models
# Subset data
All_model_data <- list(
  RLL_col_abs  = mod_colabs_rll_z,
  RLL_ext_per  = mod_extper_rll_z,
  DNR_col_abs  = mod_colabs_dnr_z,
  DNR_ext_per  = mod_extper_dnr_z
)

All_covariates <- list(
  RLL_col_abs  = RLL_col_abs_covs,
  RLL_ext_per  = RLL_ext_per_covs,
  DNR_col_abs  = DNR_col_abs_covs,
  DNR_ext_per  = DNR_ext_per_covs
)

All_responses <- list(
  RLL_col_abs  = "col_abs",
  RLL_ext_per  = "ext_per",
  DNR_col_abs  = "col_abs",
  DNR_ext_per  = "ext_per"
)

# Construct global models
global_formulas <- lapply(names(All_covariates), function(n) {
  BuildInteractions(
    response = All_responses[[n]],
    covars   = All_covariates[[n]],
    pa_int_list = pa_int_covs
  )
})
names(global_formulas) <- names(All_covariates)



## Run, Fit ###
# Individual models for each block set, response

# Helper: automate model fit, rank
options(na.action = "na.fail") # required for dredge()

RunSelection <- function(formula, data) {
  fit <- glm(formula, data = data, family = binomial)
  dredge(fit, rank = "AICc")
}

Model_selection_tables <- lapply(names(global_formulas), function(n) {
  RunSelection(global_formulas[[n]], All_model_data[[n]])
})
names(Model_selection_tables) <- names(global_formulas)


## Extract Results ##
# Helper: Reformat dredge data to long
DredgeToLong <- function(dredge_df, model_name, response_lhs) {
  
  meta_cols <- c("df", "logLik", "AICc", "delta", "weight")
  pred_cols <- setdiff(colnames(dredge_df), meta_cols)
  
  # recover formulas based on variables present (non-NA)
  # rhs = right hand side, ie. covariates in model
  rhs_formulas <- apply(dredge_df[, pred_cols, drop = FALSE], 1, function(row) {
    included <- pred_cols[!is.na(row)]
    if (length(included) == 0) "1" else paste(included, collapse = " + ")
  })
  
  full_formulas <- paste(response_lhs, "~", rhs_formulas)
  
  data.frame(
    model_name = model_name,
    model_id   = paste0(model_name, "_m", seq_len(nrow(dredge_df))),
    formula    = full_formulas,
    logLik.    = as.numeric(dredge_df$logLik),
    K          = dredge_df$df,
    AICc       = dredge_df$AICc,
    delta      = dredge_df$delta,
    weight     = dredge_df$weight,
    stringsAsFactors = FALSE
  )
}

response_map <- All_responses

Long_models_list <- lapply(names(Model_selection_tables), function(n) {
  dredge_tbl <- as.data.frame(Model_selection_tables[[n]])
  DredgeToLong(dredge_tbl, n, response_map[[n]])
})

Long_models <- do.call(rbind, Long_models_list)

names(Long_models_list) <- names(Model_selection_tables)


# Full AICc results by block set, response
RLL_col_abs_models  <- Long_models_list$RLL_col_abs
RLL_ext_per_models  <- Long_models_list$RLL_ext_per
DNR_col_abs_models  <- Long_models_list$DNR_col_abs
DNR_ext_per_models  <- Long_models_list$DNR_ext_per

View(RLL_col_abs_models)
View(RLL_ext_per_models)
View(DNR_col_abs_models)
View(DNR_ext_per_models)


# Subset AICc results
# Helper: only top models delta < 2
FilterAICc <- function(df, threshold = 2, col = "delta") {
  if (!col %in% colnames(df)) {
    stop("Column '", col, "' not found in df. Available cols: ", paste(colnames(df), collapse = ", "))
  }
  df[!is.na(df[[col]]) & df[[col]] <= threshold, , drop = FALSE]
}

# Apply
RLL_col_abs_top  <- FilterAICc(Long_models_list$RLL_col_abs, threshold = 2)
RLL_ext_per_top  <- FilterAICc(Long_models_list$RLL_ext_per, threshold = 2)
DNR_col_abs_top  <- FilterAICc(Long_models_list$DNR_col_abs, threshold = 2)
DNR_ext_per_top  <- FilterAICc(Long_models_list$DNR_ext_per, threshold = 2)

View(RLL_col_abs_top)
View(RLL_ext_per_top)
View(DNR_col_abs_top)
View(DNR_ext_per_top)

# Put in a list for convenience
Top_models_list <- list(
  RLL_col_abs = RLL_col_abs_top,
  RLL_ext_per = RLL_ext_per_top,
  DNR_col_abs = DNR_col_abs_top,
  DNR_ext_per = DNR_ext_per_top
)

# Quick diagnostics: number of top models in each set and head()
lapply(Top_models_list, function(df) {
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