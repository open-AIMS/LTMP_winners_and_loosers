## ---- load_libraries
library(tidyverse)
library(INLA)
library(readxl)
library(glmmTMB)
library(brms)
library(emmeans)
## ----end

## ---- load_functions
##setwd("R")
source("../R/helper_functions.R")
## ----end

## ---- prepare_paths_module
source("../R/02_prepare_paths.R")
## ----end

## ---- process_data
source("10_process_data.R")
## ----end

## ---- fit_models
source("20_fit_models.R")
## ----end
