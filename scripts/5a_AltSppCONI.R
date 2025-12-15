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
library(glmnet)
library(lme4)
library(MuMIn)
library(pscl)
library(AICcmodavg)

# Visualization
library(ggplot2)
library(ggfortify)
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


# SPECIES RICHNESS / EFFORT PROXY (WIP, RE-ORDER)
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

land_covars_reduced <- c("water_open_base", "shrub_scrub_base", 
                         "grass_pasture_crop_base", "developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_total_base",
                         
                         "grass_pasture_crop_diff",
                         "developed_total_diff", "forest_total_diff", "wetlands_total_diff")

covars_numeric_reduced <- c(stable_covars_all, land_covars_reduced, climate_covars_all)


# State-specific Covariate Sets
### Separate data into bins where A2 detection is either 1 or 0; not necessarily
# in precise probabilities between colo, pers, abs, ext, but more what promotes
# 'new' v. 'continued' colonization, ie. what promotes det = 1, as opposed to 
# what promotes det = 0. 
### Scaling of covars w/in subsets for relevant normalized values

spp_name <- "Red-bellied Woodpecker"

# RLL
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

# DNR
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


# CORRELATIONS, COLLINEARITY # 

# Pairwise Correlations
covar_cols_colabs_rll <- grep("_z$", names(mod_colabs_rll_z), value = TRUE)
covar_cols_extper_rll <- grep("_z$", names(mod_extper_rll_z), value = TRUE)

covar_cols_colabs_dnr <- grep("_z$", names(mod_colabs_dnr_z), value = TRUE)
covar_cols_extper_dnr <- grep("_z$", names(mod_extper_dnr_z), value = TRUE)


M1 <- cor(mod_colabs_rll_z[, covar_cols_colabs_rll], use = "pairwise.complete.obs")
M2 <- cor(mod_extper_rll_z[, covar_cols_extper_rll], use = "pairwise.complete.obs")

M3 <- cor(mod_colabs_dnr_z[, covar_cols_colabs_dnr], use = "pairwise.complete.obs")
M4 <- cor(mod_extper_dnr_z[, covar_cols_extper_dnr], use = "pairwise.complete.obs")


corrplot(M1, method = "color", tl.cex = 0.7, number.cex = 0.6)
corrplot(M2, method = "color", tl.cex = 0.7, number.cex = 0.6)

corrplot(M3, method = "color", tl.cex = 0.7, number.cex = 0.6)
corrplot(M4, method = "color", tl.cex = 0.7, number.cex = 0.6)


high_corr1 <- which(abs(M1) > 0.7 & abs(M1) < 1, arr.ind = TRUE)
apply(high_corr1, 1, function(i) cat(rownames(M1)[i[1]], "-", colnames(M1)[i[2]], "r =", M1[i[1],i[2]], "\n"))

high_corr2 <- which(abs(M2) > 0.7 & abs(M2) < 1, arr.ind = TRUE)
apply(high_corr2, 1, function(i) cat(rownames(M2)[i[1]], "-", colnames(M2)[i[2]], "r =", M2[i[1],i[2]], "\n"))

high_corr3 <- which(abs(M3) > 0.7 & abs(M3) < 1, arr.ind = TRUE)
apply(high_corr3, 1, function(i) cat(rownames(M3)[i[1]], "-", colnames(M3)[i[2]], "r =", M3[i[1],i[2]], "\n"))

high_corr4 <- which(abs(M4) > 0.7 & abs(M4) < 1, arr.ind = TRUE)
apply(high_corr4, 1, function(i) cat(rownames(M4)[i[1]], "-", colnames(M4)[i[2]], "r =", M4[i[1],i[2]], "\n"))


# Correlation Thinned Covariate Sets
factor_covars_reduced

land_covars_reduced <- c("water_open_base", "shrub_scrub_base", 
                         "grass_pasture_crop_base", "developed_total_base",
                         "forest_deciduous_base", "forest_evergreen_base", "forest_mixed_base",
                         "wetlands_total_base",
                         
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

vif_model3 <- glm(col_abs ~ ., data = mod_colabs_dnr_z[, c("col_abs", covars_numeric_reduced_z)], family = binomial)
vif(vif_model3)
alias(vif_model3)

vif_model4 <- glm(ext_per ~ ., data = mod_extper_dnr_z[, c("ext_per", covars_numeric_reduced_z)], family = binomial)
vif(vif_model4)
alias(vif_model4)

# VIF Thinned Covariate Sets
# Run 2+ x to check incoming and reduced/outgoing covariate set
# Final set prior to AICc/Model Selection process

factor_covars_reduced <- c("atlas_block")

land_covars_reduced <- c("shrub_scrub_base", "grass_pasture_crop_base", 
                         "developed_total_base", "forest_mixed_base",
                         "wetlands_total_base",
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

### --- PROCESS --- ###

### AICc from pre-screened, correlation-controlled, biologically plausible 
# candidates with both additive and interaction terms.    
### Fixed Effect: atlas_block (cannot be random b/c each block only
# has one row of data/no replication, ie. effect is not estimatable).

# MuMIn::dredge() automates candidate set construction of all main effect covar combos
# but user must manually define, limit interaction terms. User should also further
# restrict candidate models prior to running dredge(), otherwise models will not 
# converge and/or overload processor. 

### Two Steps:
# 1) "Partitioned" Main-Effects Models (MuMIn::dredge())
# A) Land-only Model
# B) Climate-only Model
# C) Stable Numeric-only Model

# 2) "Manual" Interaction Models
# Generate plausible 2-way interactions
# Fit candidates, compare with AICc


### --- STEP 1: Partition Covariate Models --- ###

out_dir <- "outputs/model_selection"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Covariates #
factor_covars_reduced # atlas_block

land_covars_reduced_z # "shrub_scrub_base", "grass_pasture_crop_base", 
# "developed_total_base", "forest_mixed_base",
# "wetlands_total_base",
# "developed_total_diff", "forest_total_diff", "wetlands_total_diff"

climate_covars_reduced_z # "tmax_38yr", "prcp_38yr", "tmax_diff", "tmin_diff", "prcp_diff"

stable_covars_reduced_z # sr_Diff, pa_percent


# Datasets #
# Objects to iterate over
mod_datasets <- list(
  colabs_rll = list(df = mod_colabs_rll_z, resp = "col_abs"),
  extper_rll = list(df = mod_extper_rll_z, resp = "ext_per"),
  colabs_dnr = list(df = mod_colabs_dnr_z, resp = "col_abs"),
  extper_dnr = list(df = mod_extper_dnr_z, resp = "ext_per")
)


# Functions
# Helper: build model formula

BuildFormula <- function(response, covariates) {
  f <- as.formula(
    paste0(response, " ~ ", paste(covariates, collapse = " + "))
  )
  environment(f) <- environment()  # safe environment
  return(f)
}


# Helper: extract top covariates, importance metric
GetTopCovars <- function(dredge_obj, delta_cut = 2, imp_cut = 0.5) {
  
  # Top models within delta AICc
  top_models <- subset(dredge_obj, delta <= delta_cut)
  if (nrow(top_models) == 0) return(character(0))
  
  # Covariate columns (exclude metadata)
  covar_cols <- setdiff(
    colnames(top_models),
    c("(Int)", "df", "logLik", "AICc", "delta", "weight")
  )
  
  # Check if there are covariates at all
  if(length(covar_cols) == 0) return(character(0))
  
  # Try model-averaged importance safely
  imp <- tryCatch(importance(dredge_obj), error = function(e) NULL)
  
  # If importance failed, just return covariates present in top models
  if(is.null(imp)) {
    top_covs <- covar_cols[apply(top_models[, covar_cols, drop = FALSE], 2, function(x) any(!is.na(x)))]
    return(top_covs)
  }
  
  # Covariates appearing in at least one top model
  present_in_top <- covar_cols[apply(top_models[, covar_cols, drop = FALSE], 2, function(x) any(!is.na(x)))]
  
  # Filter by importance threshold
  top_covs <- present_in_top[present_in_top %in% names(imp[imp >= imp_cut])]
  
  return(top_covs)
}


# Apply
options(na.action = "na.fail")
results_partitioned <- list()
top_covars_partitioned <- list()

for(ds in names(mod_datasets)) {
  
  d <- mod_datasets[[ds]]$df
  y <- mod_datasets[[ds]]$resp
  message("=== Processing dataset: ", ds, " (response: ", y, ") ===")
  
  results_partitioned[[ds]] <- list()
  top_covars_partitioned[[ds]] <- list()
  
  for(partition_name in c("land", "climate", "stable")) {
    
    # Choose covariates for this partition
    covs <- switch(partition_name,
                   land    = land_covars_reduced_z,
                   climate = climate_covars_reduced_z,
                   stable  = stable_covars_reduced_z)
    
    # Build formula & fit global model
    f <- BuildFormula(y, covs)
    m <- glm(f, data = d, family = "binomial")
    
    # Dredge
    dredged <- dredge(m, rank = "AICc", trace = FALSE)
    
    # Save dredged results
    results_partitioned[[ds]][[partition_name]] <- dredged
    
    # Extract top covariates with importance
    top_covars_partitioned[[ds]][[partition_name]] <- GetTopCovars(dredged, delta_cut = 2, imp_cut = 0.5)
    
  }
}


# Assess Top Covariates
# Helper: summarize results across both responses, block sets

SummarizePartitionedResults <- function(top_covars_list, dredge_results_list, label = "") {
  
  cat("\n==============================\n")
  cat("     SUMMARY FOR:", label, "\n")
  cat("==============================\n\n")
  
  for (ds in names(top_covars_list)) {
    
    cat("\n--------------------------------------\n")
    cat(" Dataset:", ds, "\n")
    cat("--------------------------------------\n")
    
    for (partition in names(top_covars_list[[ds]])) {
      
      cat("\n  >> Partition:", toupper(partition), "\n")
      
      # Print top covariates
      cat("     Top covariates:\n")
      print(top_covars_list[[ds]][[partition]])
      
      # Print first rows of dredge table
      cat("\n     Top dredge models (head):\n")
      suppressMessages(print(head(dredge_results_list[[ds]][[partition]])))
      
      cat("\n")
    }
  }
}


# Apply
# DNR
SummarizePartitionedResults(
  top_covars_list = top_covars_partitioned[c("colabs_dnr", "extper_dnr")],
  dredge_results_list = results_partitioned[c("colabs_dnr", "extper_dnr")],
  label = "DNR BLOCKS"
)

# RLL
SummarizePartitionedResults(
  top_covars_list = top_covars_partitioned[c("colabs_rll", "extper_rll")],
  dredge_results_list = results_partitioned[c("colabs_rll", "extper_rll")],
  label = "RLL BLOCKS"
)


# Visualize
# Helper: build formatted importance summary tables

BuildTables <- function(top_covars_list, dredge_results_list) {
  
  tables_out <- list()
  
  for (ds in names(top_covars_list)) {
    
    dredge_parts <- dredge_results_list[[ds]]
    top_parts    <- top_covars_list[[ds]]
    
    block_set <- ifelse(grepl("dnr", ds), "DNR", "RLL")
    response  <- ifelse(grepl("colabs", ds), "col_abs", "ext_per")
    
    for (partition in names(top_parts)) {
      
      dredge_obj <- dredge_parts[[partition]]
      top_covs   <- top_parts[[partition]]
      
      # Identify covariate columns
      covar_cols <- setdiff(
        colnames(dredge_obj),
        c("(Int)", "df", "logLik", "AICc", "delta", "weight")
      )
      
      # ---- SAFE importance: sum of weights across models where covariate is present ----
      imp <- sapply(covar_cols, function(v) {
        present <- !is.na(dredge_obj[[v]])
        sum(dredge_obj$weight[present])
      })
      
      df <- data.frame(
        block_set = block_set,
        response  = response,
        partition = partition,
        covariate = names(imp),
        importance = as.numeric(imp),
        in_top_models = names(imp) %in% top_covs,
        rank = rank(-imp, ties.method = "first"),
        stringsAsFactors = FALSE
      )
      
      tbl_name <- paste(block_set, response, partition, sep = "_")
      tables_out[[tbl_name]] <- df
    }
  }
  
  return(tables_out)
}


# Apply
### Model-Averaged Variable Importance Score: sum of AIcc weights of all (dredged)
# models wherein a given covariate appears ("relative importance," Burnham & Anderson 2002).
# If all the highest-weight models include X, then importance ~= 1
# If some mid-weight models include X, then importance ~= 0.3-0.7
# If only low-weight models include X, then importance ~= 0.05-0.2
# If no top models include X, then importance ~= 0

formatted_tables <- BuildTables(
  top_covars_list = top_covars_partitioned,
  dredge_results_list = results_partitioned
)

# AICc, Importance-thinned Covariate Sets
### Burnham & Anderson 2002:
# Importance ≥ 0.80	--> Strong evidence the covariate belongs in the final model(s).
# Importance ~= 0.50–0.80	--> Moderate evidence; include if biologically meaningful.

# RLL, ColAbs
formatted_tables$RLL_col_abs_land
# >= 0.8: forest_mixed_base_z, developed_total_base_z, grass_pasture_crop_base_z, wetlands_total_base_z
# ~= 0.5-0.8: shrub_scrub_base_z, forest_total_diff_z, wetlands_total_diff_z
formatted_tables$RLL_col_abs_climate
# >= 0.8: tmax_38yr_z, prcp_38yr_z, prcp_diff_z
# ~= 0.5-0.8: none
formatted_tables$RLL_col_abs_stable
# >= 0.8: pa_percent_z, sr_Diff_z
# ~= 0.5-0.8: n/a

# RLL, ExtPer
formatted_tables$RLL_ext_per_land
# >= 0.8: none
# ~= 0.5-0.8: forest_total_diff_z, grass_pasture_crop_base_z
formatted_tables$RLL_ext_per_climate
# >= 0.8: none
# ~= 0.5-0.8: tmax_38yr_z, prcp_diff_z, tmin_diff_z
formatted_tables$RLL_ext_per_stable 
# >= 0.8: sr_Diff_z
# ~= 0.5-0.8: pa_percent_z 


# DNR, ColAbs
formatted_tables$DNR_col_abs_land
# >= 0.8: forest_mixed_base_z, developed_total_base_z, wetlands_total_base_z
# ~= 0.5-0.8: none
formatted_tables$DNR_col_abs_climate
# >= 0.8: tmax_38yr_z, prcp_38yr_z
# ~= 0.5-0.8: none
formatted_tables$DNR_col_abs_stable
# >= 0.8: pa_percent_z, sr_Diff_z
# ~= 0.5-0.8: n/a

# DNR, ExtPer
formatted_tables$DNR_ext_per_land
# >= 0.8: none
# ~= 0.5-0.8: forest_total_diff_z, developed_total_base_z
formatted_tables$DNR_ext_per_climate
# >= 0.8: tmax_38yr_z
# ~= 0.5-0.8: ptmax_diff_z
formatted_tables$DNR_ext_per_stable
# >= 0.8: pa_percent_z
# ~= 0.5-0.8: sr_Diff_z

################################ WIP FOR NOW #################################
### Don't use as another importance-based reduction step for now--think I'm oversimplifying model
# and creating too perfect of a fit (based on competative models in global ranking
# and ID of uninformative parameters)
# Ie. Just bring thru same covariates that went in


RLL_col_abs_covs <- c(
  "shrub_scrub_base_z",
  "grass_pasture_crop_base_z",
  "developed_total_base_z",
  "forest_mixed_base_z", 
  "wetlands_total_base_z", 
  "developed_total_base_z", 
  "forest_total_diff_z",
  "wetlands_total_diff_z",
  "tmax_38yr_z",
  "prcp_38yr_z",
  "tmax_diff_z",
  "tmin_diff_z",
  "prcp_diff_z",
  "pa_percent_z",
  "sr_Diff_z"
)

RLL_ext_per_covs <- c(
  "shrub_scrub_base_z",
  "grass_pasture_crop_base_z",
  "developed_total_base_z",
  "forest_mixed_base_z", 
  "wetlands_total_base_z", 
  "developed_total_base_z", 
  "forest_total_diff_z",
  "wetlands_total_diff_z",
  "tmax_38yr_z",
  "prcp_38yr_z",
  "tmax_diff_z",
  "tmin_diff_z",
  "prcp_diff_z",
  "pa_percent_z",
  "sr_Diff_z"
)


DNR_col_abs_covs <- c(
  "shrub_scrub_base_z",
  "grass_pasture_crop_base_z",
  "developed_total_base_z",
  "forest_mixed_base_z", 
  "wetlands_total_base_z", 
  "developed_total_base_z", 
  "forest_total_diff_z",
  "wetlands_total_diff_z",
  "tmax_38yr_z",
  "prcp_38yr_z",
  "tmax_diff_z",
  "tmin_diff_z",
  "prcp_diff_z",
  "pa_percent_z",
  "sr_Diff_z"
)

DNR_ext_per_covs <- c(
  "shrub_scrub_base_z",
  "grass_pasture_crop_base_z",
  "developed_total_base_z",
  "forest_mixed_base_z", 
  "wetlands_total_base_z", 
  "developed_total_base_z", 
  "forest_total_diff_z",
  "wetlands_total_diff_z",
  "tmax_38yr_z",
  "prcp_38yr_z",
  "tmax_diff_z",
  "tmin_diff_z",
  "prcp_diff_z",
  "pa_percent_z",
  "sr_Diff_z"
)


# Visualize

# Data formatting
all_imp <- dplyr::bind_rows(formatted_tables, .id = "table_name") # combine into single long table

all_imp <- all_imp %>% # cleaner id for x-axis
  dplyr::mutate(
    block_response_part = paste(block_set, response, partition, sep = "_"),
    covariate = factor(covariate)
  )

covs <- unique(all_imp$covariate) # covariate object

separator_df <- data.frame( # aesthetic: white space separating block sets
  block_response_part = "separator",
  covariate = covs,
  importance = NA,
  block_set = NA,
  response = NA,
  partition = NA
)

all_imp_sep <- dplyr::bind_rows(all_imp, separator_df) # aesthetic: put white space in table

all_imp_sep$block_response_part <- factor( # order factor levels
  all_imp_sep$block_response_part,
  levels = c(
    "RLL_col_abs_land",
    "RLL_col_abs_climate",
    "RLL_col_abs_stable",
    "RLL_ext_per_land",
    "RLL_ext_per_climate",
    "RLL_ext_per_stable",
    
    "separator",
    
    "DNR_col_abs_land",
    "DNR_col_abs_climate",
    "DNR_col_abs_stable",
    "DNR_ext_per_land",
    "DNR_ext_per_climate",
    "DNR_ext_per_stable"
  )
)

x_labels_pretty <- c( # aesthetic: asix labels
  
  "RLL_col_abs_land"    = "RLL: Colonization (Land)",
  "RLL_col_abs_climate" = "RLL: Colonization (Climate)",
  "RLL_col_abs_stable"  = "RLL: Colonization (Stable)",
  "RLL_ext_per_land"    = "RLL: Extinction (Land)",
  "RLL_ext_per_climate" = "RLL: Extinction (Climate)",
  "RLL_ext_per_stable"  = "RLL: Extinction (Stable)",
  
  "separator" = "",
  
  "DNR_col_abs_land"    = "DNR: Colonization (Land)",
  "DNR_col_abs_climate" = "DNR: Colonization (Climate)",
  "DNR_col_abs_stable"  = "DNR: Colonization (Stable)",
  "DNR_ext_per_land"    = "DNR: Extinction (Land)",
  "DNR_ext_per_climate" = "DNR: Extinction (Climate)",
  "DNR_ext_per_stable"  = "DNR: Extinction (Stable)"
)

y_labels <- c( # aesthetic: y-axis labels
  
  "(Intercept)"              = "(Intercept)",
  "developed_total_base_z"  = "Developed (%)",
  "forest_deciduous_base_z" = "Deciduous Forest (%)",
  "forest_evergreen_base_z" = "Evergreen Forest (%)",
  "forest_mixed_base_z"     = "Mixed Forest (%)",
  "forest_total_diff_z"     = "Forest Change (%)",
  "grassland_base_z"        = "Grassland (%)",
  "shrub_scrub_base_z"      = "Shrub/Scrub (%)",
  "wetlands_total_base_z"   = "Wetlands (%)",
  "prcp_38yr_z"             = "Avg Precip",
  "prcp_diff_z"             = "Δ Precip",
  "tmax_38yr_z"             = "Avg Tmax",
  "tmax_diff_z"             = "Δ Tmax",
  "tmin_diff_z"             = "Δ Tmin",
  "pa_percent_z"            = "Protected Area (%)",
  "sr_Diff_z"               = "Δ Species Richness"
)

# ggplot() Vis
# Heatmap of importance score by covariate, block set, response
ggplot(all_imp_sep, aes(x = block_response_part,
                        y = covariate,
                        fill = importance)) +
  geom_tile(color = NA) +   # <--- No borders anywhere
  scale_fill_viridis_c(option = "viridis", name = "Importance",
                       na.value = "white") +
  scale_x_discrete(labels = x_labels_pretty) +
  scale_y_discrete(labels = y_labels) + 
  labs(
    x = "Model Set (Response, Block Set)",
    y = "Covariate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )







### --- STEP 2: "Manual" Interaction Models --- ###

### Carry-ins ###
# Model data
mod_colabs_rll_z
mod_extper_rll_z
mod_colabs_dnr_z
mod_extper_dnr_z

# Thinned covar sets by block set, response
RLL_col_abs_covs
RLL_ext_per_covs
DNR_col_abs_covs
DNR_ext_per_covs


### Model Construction ###

# Interactions #
# Custom list
pa_int_covs <- c(
  "developed_total_base_z",
  "forest_total_diff_z",
  "wetlands_total_base_z",
  "forest_mixed_base_z",
  "tmax_38yr_z",
  "tmin_diff_z"
  #"sr_Diff_z"
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

# Construct model
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
  
  # recover RHS formulas based on variables present (non-NA)
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





# Visualize
### Tidier AICc tables
PrettyModels <- function(df) {
  if (nrow(df) == 0) return(df)
  
  df |>
    dplyr::mutate(
      AICc  = round(AICc, 3),
      delta = round(delta, 3),
      weight = round(weight, 4)
    ) |>
    dplyr::rename(
      Model      = model_id,
      Set        = model_name,
      Formula    = formula,
      ΔAICc      = delta,
      Weight     = weight
    ) |>
    dplyr::select(Set, Model, Formula, AICc, ΔAICc, Weight)
}

Pretty_list <- lapply(Top_models_list, PrettyModels)
Pretty_list


