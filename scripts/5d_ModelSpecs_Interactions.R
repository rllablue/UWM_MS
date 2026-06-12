## WIP Notes: 
# 6.8: Separated/Created script; need to add int terms, e.g. PA:ENV, PA:SAC



################################ PA: Interactions
mod_col_pa <- pa_glm_models[["RLL_col"]]
mod_ext_pa <- pa_glm_models[["RLL_ext"]]

form_col_base <- formula(mod_col_pa)
form_ext_base <- formula(mod_ext_pa)

form_col_int <- update(form_col_base, . ~ . + pa_prop:forest_total_base) # :tmax_38yr, :forest_total_base, :grassland_base, :wetlands_total_base
form_ext_int <- update(form_ext_base, . ~ . + pa_prop:forest_total_base)


mod_col_pa_int <- glm(
  formula = form_col_int,
  data    = mod_col_rll,
  family  = binomial
)
summary(mod_col_pa_int)
car::vif(mod_col_pa_int, type = "terms")


mod_ext_pa_int <- glm(
  formula = form_ext_int,
  data    = mod_ext_rll,
  family  = binomial
)
summary(mod_ext_pa_int)
car::vif(mod_ext_pa_int, type = "terms")






