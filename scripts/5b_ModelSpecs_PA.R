### WIP NOTES
# 6.8: Added model fit/diagnostics section; need to store important model coefficients
# (for ease of access during results writing); need to re-build out PA sub-type covariate 
# terms, e.g. GAP status, PA size, PA connectivity




### Script: 5a_ModelSpecs_PA.R
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

# Model Objects: global_glm_models

# Data Objects: data_dir <- list(
#RLL_col = modcol_rll_z_df,
#RLL_ext = modext_rll_z_df
#)





######################
### BASE PA MODELS ###
######################

pa_covar <- c("pa_prop")


### --- A: Construct PA Models
### Add PA term/s to top candidates from ENV model script, i.e. global_glm_models

## Helper: parse ENV model string, add PA covariate
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




########################
### MODEL ASSESSMENT ###
########################

models <- pa_glm_models


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




##########################
### PA CHARACTERISTICS ###
##########################


### GAP STATUS #################################################################



### --- PROTECTION STRINGENCY --- ###

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



### --- PA MANAGER --- ###

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











########################
### MODEL ASSESSMENT ###
########################

### --- 1: STATS CHECKS --- ###

### --- 1A: Correlated Pairs
GetCorrelatedPairs <- function(model,
                               cutoff = 0.7,
                               digits = 2) {
  
  mm <- model.matrix(model)
  
  mm <- mm[
    ,
    colnames(mm) != "(Intercept)",
    drop = FALSE
  ]
  
  cor_mat <- cor(
    mm,
    use = "complete.obs"
  )
  
  cor_df <- as.data.frame(
    as.table(cor_mat)
  )
  
  names(cor_df) <- c(
    "var1",
    "var2",
    "correlation"
  )
  
  cor_df <- cor_df %>%
    filter(var1 != var2) %>%
    rowwise() %>%
    mutate(
      pair = paste(
        sort(c(var1, var2)),
        collapse = "__"
      )
    ) %>%
    ungroup() %>%
    distinct(
      pair,
      .keep_all = TRUE
    ) %>%
    dplyr::select(-pair)
  
  cor_df %>%
    filter(
      abs(correlation) >= cutoff
    ) %>%
    mutate(
      correlation = round(
        correlation,
        digits
      )
    ) %>%
    arrange(
      desc(abs(correlation))
    )
  
}

global_cor_pairs <- lapply(
  global_glm_models,
  GetCorrelatedPairs,
  cutoff = 0.7
)

global_cor_pairs



### --- 1B: VIF
vif_global_models <- lapply(
  global_glm_models,
  car::vif
)

vif_global_models



### --- 2: MODEL FIT --- ###

### --- 2A: AUC

global_auc <- lapply(
  names(global_glm_models),
  function(nm) {
    
    model <- global_glm_models[[nm]]
    
    observed <- model$model[[1]]
    
    predicted <- predict(
      model,
      type = "response"
    )
    
    auc_value <- pROC::auc(
      observed,
      predicted
    )
    
    data.frame(
      model = nm,
      AUC = as.numeric(
        auc_value
      )
    )
  }
)

global_auc <- bind_rows(
  global_auc
)

global_auc



### --- 2B: Pseudo-R2
global_r2 <- lapply(
  names(global_glm_models),
  function(nm) {
    
    model <- global_glm_models[[nm]]
    
    r2_vals <- pscl::pR2(model)
    
    data.frame(
      model = nm,
      McFadden_R2 =
        r2_vals["McFadden"]
    )
    
  }
)

global_r2 <- bind_rows(
  global_r2
)

global_r2



### --- 2C: Over-dispersion
### disp. ratio approx. = 1 (none); disp. ration > 1.2 (possible); p < 0.05 (sig. over-dispersion)

overdisp_results <- bind_rows(
  
  lapply(
    names(global_glm_models),
    function(nm) {
      
      chk <- performance::check_overdispersion(
        global_glm_models[[nm]]
      )
      
      tibble(
        model = nm,
        dispersion_ratio = chk$dispersion_ratio,
        p_value = chk$p_value
      )
      
    }
  )
)

overdisp_results



### --- RESIDUALS --- ###

### --- 1A: DHARMa Residuals
set.seed(123)

global_dharma <- lapply(
  global_glm_models,
  DHARMa::simulateResiduals,
  plot = FALSE
)


### --- 1B: Residual Uniformity
### Checks whether residuals follow expected uniform distribution;
# p > 0.05 (acceptable distribution), p < 0.05 (possible model misspecification) 

uniformity_results <- lapply(
  names(global_dharma),
  function(nm) {
    
    tst <- DHARMa::testUniformity(
      global_dharma[[nm]]
    )
    
    data.frame(
      model = nm,
      statistic = unname(tst$statistic),
      p_value = tst$p.value
    )
    
  }
)

uniformity_results <- bind_rows(
  uniformity_results
)

uniformity_results



### --- 1C: Residual Dispersion
### p > 0.05 (dispersion acceptable), p < 0.005 (possible over- or under-dispersion)

dispersion_results <- lapply(
  names(global_dharma),
  function(nm) {
    
    tst <- DHARMa::testDispersion(
      global_dharma[[nm]]
    )
    
    data.frame(
      model = nm,
      statistic = unname(tst$statistic),
      p_value = tst$p.value
    )
    
  }
)

dispersion_results <- bind_rows(
  dispersion_results
)

dispersion_results


### --- 1D: Residual Outlier Test
outlier_results <- bind_rows(
  
  lapply(
    names(global_dharma),
    function(nm) {
      
      tst <- DHARMa::testOutliers(
        global_dharma[[nm]]
      )
      
      tibble(
        model = nm,
        
        outlier_frequency =
          if (!is.null(tst$estimate)) {
            as.numeric(tst$estimate)
          } else if (!is.null(tst$outliers)) {
            mean(tst$outliers, na.rm = TRUE)
          } else {
            NA_real_
          },
        
        p_value =
          if (!is.null(tst$p.value)) {
            tst$p.value
          } else {
            NA_real_
          },
        
        method =
          if (!is.null(tst$method)) {
            tst$method
          } else {
            NA_character_
          }
      )
      
    }
  )
)

outlier_results




### --- 5E: Visualizations

### Diagnostic Plots
par(mfrow = c(1, 2))

plot(
  global_dharma$RLL_col,
  main = "Colonization"
)

plot(
  global_dharma$RLL_ext,
  main = "Extirpation"
)

par(mfrow = c(1, 1))



### ROC Curves
global_predictions <- lapply(
  global_glm_models,
  predict,
  type = "response"
)

roc_col <- pROC::roc(
  response =
    global_glm_models$RLL_col$model[[1]],
  predictor =
    global_predictions$RLL_col
)

roc_ext <- pROC::roc(
  response =
    global_glm_models$RLL_ext$model[[1]],
  predictor =
    global_predictions$RLL_ext
)


# dataframe for plotting
roc_col_df <- data.frame(
  FPR = 1 - roc_col$specificities,
  TPR = roc_col$sensitivities,
  Model = "Colonization"
)

roc_ext_df <- data.frame(
  FPR = 1 - roc_ext$specificities,
  TPR = roc_ext$sensitivities,
  Model = "Extirpation"
)

roc_df <- bind_rows(
  roc_col_df,
  roc_ext_df
)

auc_col <- round(
  as.numeric(
    pROC::auc(roc_col)
  ),
  3
)

auc_ext <- round(
  as.numeric(
    pROC::auc(roc_ext)
  ),
  3
)


# plot
ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  scale_color_manual(
    values = c(
      "Colonization" = "#2C7BB6",
      "Extinction"   = "#D7191C"
    ),
    labels = c(
      paste0("Colonization (AUC = ", auc_col, ")"),
      paste0("Extirpation (AUC = ", auc_ext, ")")
    )
  ) +
  
  labs(
    title = "ROC Curves for ENV Logistic GLMs",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL
  ) +
  
  coord_equal() +
  
  theme_minimal(base_size = 13) +
  
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

