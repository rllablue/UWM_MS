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




#####################################
### MODEL CONSTRUCTION, SELECTION ### 
#####################################

### --- PROCESS, SET UP --- ###

### Steps ###

# 1) Partitioned Additive ENV Model Selection
# Construct, fit, rank separate models partitioning land and climate covariates;
# carry over best covariates for each set into full models.

# 2) Full Model Selection
# Construct, fit, rank top models AICc < 2 w/ additive land, climate, effort, PA terms; interaction terms w PA;
# extract reference model for each response subset.

# 3) Full Model Diagnostics
# Assess top candidate models for each response for correlations, dispersion, fit, uninformative parameters, etc.




### --- STEP 1: ENV PARTITIONED MODELS --- ###

### Full, additive candidate model sets for climate and land cover covariates
# w/ null included. AICc fitting, ranking; extract covars to pass on from "reference" 
# model (ie. delta = 0) to global model selection steps.


effort_covs <- c("sr_diff", "priority_block")
sac_cov <- "sac_kernel"
pa_cov <- "pa_prop"


# lookup grid for model response-specific cov subsets
partition_grid <- expand.grid(
  response = c("col", "ext"),
  blocks = "RLL",
  partition = c("climate", "land"),
  stringsAsFactors = FALSE
)


### Helper: generate all additive predictor combos w/in partitioned covariate sets
# required = structural effort- and spatial-control variable
BuildPartitionedRHS <- function(covariates,
                                effort_covs,
                                sac_cov) {
  
  covariates <- unique(covariates)
  
  forced <- unique(c(effort_covs, sac_cov))
  
  optional <- setdiff(covariates, forced)
  
  rhs <- paste(forced, collapse = " + ")
  
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


### Helper: Fit, rank partitioned candidate models
FitPartitionedModels <- function(response,
                                 data,
                                 covariates,
                                 effort_covs,
                                 sac_cov,
                                 family = binomial,
                                 include_null = TRUE) {
  
  rhs_terms <- BuildPartitionedRHS(
    covariates  = covariates,
    effort_covs = effort_covs,
    sac_cov     = sac_cov
  )
  
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
    effort_covs  = effort_covs,
    sac_cov      = sac_cov,
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








### Helper: tidy table of top models (delta <= 2)
ExtractTopModelTable <- function(model_sel_obj, delta = 2) {
  
  top_mods <- subset(model_sel_obj, delta <= delta)
  
  out <- as.data.frame(top_mods)
  
  out$Model <- rownames(out)
  
  rownames(out) <- NULL
  
  out |>
    dplyr::transmute(
      Rank      = dplyr::row_number(),
      Model,
      Covariates = Model,
      K         = df,
      LogLik    = logLik,
      AICc,
      Delta_AICc = delta,
      Weight    = weight
    )
}

top_part_tables <- lapply(
  top_part_models,
  ExtractTopModelTable,
  delta = 2
)



### Helper: model covs 1/0 matrix
ExtractModelMatrix <- function(top_tbl,
                               forced_covs = c("sr_diff",
                                               "priority_block",
                                               "sac_kernel")) {
  
  cov_lists <- strsplit(top_tbl$Model, " \\+ ")
  
  all_covs <- sort(unique(unlist(cov_lists)))
  all_covs <- setdiff(all_covs, forced_covs)
  
  mat <- sapply(all_covs, function(v)
    sapply(cov_lists, function(x) as.integer(v %in% x))
  )
  
  cbind(
    top_tbl[, c(
      "Rank",
      "K",
      "LogLik",
      "AICc",
      "Delta_AICc",
      "Weight"
    )],
    as.data.frame(mat)
  )
}

top_part_matrices <- lapply(
  top_part_tables,
  ExtractModelMatrix
)

top_part_matrices



### Helper: variable frequency table
CovariateSupport <- function(top_tbl,
                             forced_covs = c("sr_diff",
                                             "priority_block",
                                             "sac_kernel")) {
  
  cov_lists <- strsplit(top_tbl$Model, " \\+ ")
  
  all_covs <- sort(unique(unlist(cov_lists)))
  
  # remove forced terms if desired
  all_covs <- setdiff(all_covs, forced_covs)
  
  support_tbl <- lapply(all_covs, function(v) {
    
    present <- sapply(cov_lists, function(x) v %in% x)
    
    data.frame(
      Covariate       = v,
      Model_Frequency = sum(present),
      Weight_Sum      = sum(top_tbl$Weight[present]),
      stringsAsFactors = FALSE
    )
    
  }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(
      dplyr::desc(Weight_Sum),
      dplyr::desc(Model_Frequency)
    )
  
  support_tbl
}

part_cov_support <- lapply(
  top_part_tables,
  CovariateSupport
)

part_cov_support








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


### Helper: Extract reference model covariates
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
merged_ref_covariates







### --- STEP 2: ENV MODEL SELECTION --- #################################################################################### WIP

### Combine climate, land covars lists from partitioned models into single cov pool 
# per response; carry into bext step to add PA, PA:interactions

### Helper: 
BuildGlobalRHS <- function(covariates,
                           effort_covs,
                           sac_cov) {
  
  covariates <- unique(covariates)
  
  forced <- unique(c(effort_covs, sac_cov))
  
  optional <- setdiff(covariates, forced)
  
  rhs <- paste(forced, collapse = " + ")
  
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


### Helper: 
FitGlobalModels <- function(response,
                            data,
                            covariates,
                            effort_covs,
                            sac_cov,
                            family = binomial,
                            include_null = TRUE) {
  
  rhs_terms <- BuildGlobalRHS(
    covariates  = covariates,
    effort_covs = effort_covs,
    sac_cov     = sac_cov
  )
  
  models <- list()
  
  if (include_null) {
    models[["NULL"]] <- glm(
      as.formula(paste(response, "~ 1")),
      data = data,
      family = family
    )
  }
  
  for (rhs in rhs_terms) {
    
    models[[rhs]] <- glm(
      as.formula(paste(response, "~", rhs)),
      data = data,
      family = family
    )
    
  }
  
  MuMIn::model.sel(models, rank = "AICc")
}



# Filter out response x block models w 0 covars
valid_keys <- names(merged_ref_covariates)[
  lengths(merged_ref_covariates) > 0
]

### Apply
global_models <- lapply(valid_keys, function(key) { 
  
  FitGlobalModels(
    response    = strsplit(key, "_")[[1]][2],
    data        = data_dir[[key]],
    covariates  = merged_ref_covariates[[key]],
    effort_covs = effort_covs,
    sac_cov     = sac_cov,
    family      = binomial
  )
  
})

names(global_models) <- valid_keys

top_global_models <- lapply(global_models, ExtractTopModels, delta = 2)
reference_global_models <- lapply(top_global_models, ExtractReferenceModel) 








### Helper: tidy table of top models (delta <= 2)
ExtractTopModelTable <- function(model_sel_obj, delta = 2) {
  
  top_mods <- subset(model_sel_obj, delta <= delta)
  
  out <- as.data.frame(top_mods)
  
  out$Model <- rownames(out)
  
  rownames(out) <- NULL
  
  out |>
    dplyr::transmute(
      Rank      = dplyr::row_number(),
      Model,
      Covariates = Model,
      K         = df,
      LogLik    = logLik,
      AICc,
      Delta_AICc = delta,
      Weight    = weight
    )
}

top_global_tables <- lapply(
  top_global_models,
  ExtractTopModelTable,
  delta = 2
)



### Helper: model covs 1/0 matrix
ExtractModelMatrix <- function(top_tbl,
                               forced_covs = c("sr_diff",
                                               "priority_block",
                                               "sac_kernel")) {
  
  cov_lists <- strsplit(top_tbl$Model, " \\+ ")
  
  all_covs <- sort(unique(unlist(cov_lists)))
  all_covs <- setdiff(all_covs, forced_covs)
  
  mat <- sapply(all_covs, function(v)
    sapply(cov_lists, function(x) as.integer(v %in% x))
  )
  
  cbind(
    top_tbl[, c(
      "Rank",
      "K",
      "LogLik",
      "AICc",
      "Delta_AICc",
      "Weight"
    )],
    as.data.frame(mat)
  )
}

top_global_matrices <- lapply(
  top_global_tables,
  ExtractModelMatrix
)

top_global_matrices



### Helper: variable frequency table
CovariateSupport <- function(top_tbl,
                             forced_covs = c("sr_diff",
                                             "priority_block",
                                             "sac_kernel")) {
  
  cov_lists <- strsplit(top_tbl$Model, " \\+ ")
  
  all_covs <- sort(unique(unlist(cov_lists)))
  
  # remove forced terms if desired
  all_covs <- setdiff(all_covs, forced_covs)
  
  support_tbl <- lapply(all_covs, function(v) {
    
    present <- sapply(cov_lists, function(x) v %in% x)
    
    data.frame(
      Covariate       = v,
      Model_Frequency = sum(present),
      Weight_Sum      = sum(top_tbl$Weight[present]),
      stringsAsFactors = FALSE
    )
    
  }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(
      dplyr::desc(Weight_Sum),
      dplyr::desc(Model_Frequency)
    )
  
  support_tbl
}

global_cov_support <- lapply(
  top_global_tables,
  CovariateSupport
)

global_cov_support








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
