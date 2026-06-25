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
  dplyr, tidyverse, magrittr, purrr, stringr, readxl,
  
  # Diagnostics
  nnet, stats, broom, corrplot, car, performance, DescTools, Metrics, pROC, rsample,
  lme4, pscl, AICcmodavg, MuMIn, arm, ncf, DHARMa, psych, usdm, glmmTMB,
  
  # Visualization
  ggplot2, viridis, gt, sf, gpkg, webshot2
  
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

### --- A: Construct PA Models
### Add PA term/s to top candidates from ENV model script, i.e. global_glm_models

pa_glm_models <- lapply(names(global_glm_models), function(nm) {
  
  glm(
    formula = update(
      formula(global_glm_models[[nm]]),
      . ~ . + pa_prop
    ),
    data = data_dir[[nm]],
    family = binomial
  )
  
})

names(pa_glm_models) <- names(global_glm_models)

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
    title = "ROC Curves for PA Logistic GLMs",
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

### GAP STATUS ############################################################################################################################ WIP



### --- PROTECTION STRINGENCY --- ###

gap_covars <- c("gap1_prop", "gap2_prop", "gap3_prop", "gap12combo_prop")


pasub_glm_models <- lapply(names(pa_glm_models), function(nm) {
  
  glm(
    formula = update(
      formula(pa_glm_models[[nm]]),
      . ~ . + gap12combo_prop + gap3_prop
    ),
    data = data_dir[[nm]],
    family = binomial
  )
  
})

names(pasub_glm_models) <- names(pa_glm_models)

pasub_glm_summaries <- lapply(pasub_glm_models, summary)
pasub_glm_summaries






########################
### MODEL ASSESSMENT ###
########################

models <- pasub_glm_models


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
    title = "ROC Curves for PA Logistic GLMs",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))


















### --- PA MANAGER --- #########################################################

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