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

# 1) Partitioned AICc Model Selection
# Construct, fit, rank separate models partitioning land and climate covariates;
# carry over best covariates for each set into full ENV models.

# 2) Additive ENV Model Selection
# Construct, fit, rank top models AICc < 2 w/ additive land, climate, effort terms;
# extract reference model for each response subset.

# 3) Model Diagnostics
# Assess top candidate models for each response for correlations, dispersion, fit, uninformative parameters, etc.








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
required_covs <- c("sr_diff", "priority_block", "sac_kernel")


FitPartitionedModels <- function(response, 
                                 data, 
                                 covariates, 
                                 required = required_covs,
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
merged_ref_covariates




### --- STEP 2: ENV MODEL SELECTION --- ###

### Combine climate, land covars lists from partitioned models into single covar list/pool for 
# new additive mod selection process; find new ref model to carry into interaction step (2B).

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


# Predicted Probabilities
# histograms
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

for (nm in names(models)) {
  hist(
    predict(models[[nm]], type = "response"),
    main = paste("Predicted probs:", nm),
    xlab = "Probability",
    breaks = 30
  )
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


