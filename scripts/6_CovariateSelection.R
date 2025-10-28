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



###########################
### COVARIATE SELECTION ###
###########################

### Use penalized regression techniques and other threshold structures for 
# determining which covariate sets to use for each species in modeling protocol







