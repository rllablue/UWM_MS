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

# Model and covariate fitting, evaluation
library(nnet)
library(stats)
library(broom)
library(car)
library(performance)
library(DescTools)
library(Metrics)
library(pROC)
library(rsample)
library(glmnet)
library(forcats)

# Spatial
library(sf)
library(units)

# Visualization
library(ggplot2)
library(viridis)
library(RColorBrewer)



########################
### MODEL COMPARISON ###
########################

### --- EXTRACT EVALUATION METRICS --- ###

Compare_models <- function(species_name, zf_summary, covariates) {
  
  # --- Fit multinomial model --- #
  multi_results <- Model_multinom(species_name, zf_summary, covariates)
  
  # Extract fit metrics
  multi_aic <- AIC(multi_results$model)
  multi_r2  <- DescTools::PseudoR2(multi_results$model, which = "McFadden")
  
  # --- Fit three separate logistic models --- #
  log_results <- Model_three_logistic(species_name, zf_summary, covariates)
  
  # Extract AIC and pseudo-R2 per logistic model
  log_metrics <- map(log_results$models, ~list(
    AIC = AIC(.x),
    pseudoR2 = DescTools::PseudoR2(.x, which = "McFadden")
  ))
  
  # --- Combine predicted probabilities --- #
  pred_comparison <- multi_results$data %>%
    select(atlas_block, transition_state) %>%
    mutate(row_id = row_number()) %>%
    left_join(
      multi_results$pred_probs %>% as.data.frame() %>% mutate(row_id = row_number()) %>%
        rename_with(~paste0("multi_", .), -row_id),
      by = "row_id"
    ) %>%
    left_join(log_results$pred_probs, by = "atlas_block") %>%
    select(-row_id)
  
  # --- Return results --- #
  return(list(
    multinomial = list(model = multi_results$model, pred_probs = multi_results$pred_probs, AIC = multi_aic, pseudoR2 = multi_r2),
    logistic = list(models = log_results$models, pred_probs = log_results$pred_probs, metrics = log_metrics),
    comparison_df = pred_comparison
  ))
}


### --- APPLY --- ###

# Example: Red-bellied Woodpecker
rbwo_comparison <- Compare_models("Red-bellied Woodpecker", breeders_zf_summary, wibba_modeling_comp)

# --- Inspect metrics --- #
rbwo_comparison$multinomial$AIC
rbwo_comparison$multinomial$pseudoR2

rbwo_comparison$logistic$metrics$colonization
rbwo_comparison$logistic$metrics$extinction
rbwo_comparison$logistic$metrics$persistence

# --- Inspect predicted probabilities --- #
head(rbwo_comparison$comparison_df)


