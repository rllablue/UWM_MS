## WIP Notes: 
# 6.8: Separated/Created script; need to add int terms, e.g. PA:ENV, PA:SAC

### --- Simultaneously fitted models
int_glm_models <- lapply(names(sac_glm_models), function(nm) {
  
  glm(
    formula = update(
      formula(sac_glm_models[[nm]]),
      . ~ . + pa_prop*sac_kernel + pa_prop*tmax_38yr + pa_prop*forest_deciduous_base + pa_prop*wetlands_woody_base
    ),
    data = data_dir[[nm]],
    family = binomial
  )
  
})

names(int_glm_models) <- names(sac_glm_models)

int_glm_summaries <- lapply(int_glm_models, summary)
int_glm_summaries


# w pa subsets
int2_glm_models <- lapply(names(sac2_glm_models), function(nm) { # or pasub_glm_models
  
  glm(
    formula = update(
      formula(sac2_glm_models[[nm]]),
      . ~ . + pa_prop*sac_kernel + pa_prop*tmax_38yr + pa_prop*forest_deciduous_base + pa_prop*wetlands_woody_base
    ),
    data = data_dir[[nm]],
    family = binomial
  )
  
})

names(int2_glm_models) <- names(sac2_glm_models)

int2_glm_summaries <- lapply(int2_glm_models, summary)
int2_glm_summaries








### --- Individually fitted models
sac_glm_models_col <- sac_glm_models[grepl("_col$", names(sac_glm_models))]
sac_glm_models_ext <- sac_glm_models[grepl("_ext$", names(sac_glm_models))]

data_dir_col <- data_dir[grepl("_col$", names(data_dir))]
data_dir_ext <- data_dir[grepl("_ext$", names(data_dir))]


# Colonization
int_glm_models_col <- lapply(names(sac_glm_models_col), function(nm) {
  
  glm(
    formula = update(
      formula(sac_glm_models_col[[nm]]),
      . ~ . + pa_prop*sac_kernel
    ),
    data = data_dir_col[[nm]],
    family = binomial
  )
  
})

names(int_glm_models_col) <- names(sac_glm_models_col)
int_glm_summaries_col <- lapply(int_glm_models_col, summary)
int_glm_summaries_col


# Extirpation
int_glm_models_ext <- lapply(names(sac_glm_models_ext), function(nm) {
  
  glm(
    formula = update(
      formula(sac_glm_models_ext[[nm]]),
      . ~ . + pa_prop*sac_kernel
    ),
    data = data_dir_ext[[nm]],
    family = binomial
  )
  
})

names(int_glm_models_ext) <- names(sac_glm_models_ext)
int_glm_summaries_ext <- lapply(int_glm_models_ext, summary)
int_glm_summaries_ext






########################
### MODEL ASSESSMENT ###
########################

models <- int_glm_models


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
    title = "ROC Curves for SAC Logistic GLMs",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))




