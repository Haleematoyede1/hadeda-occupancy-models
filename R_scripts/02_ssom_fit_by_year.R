# ---------------------------------------------------------
# Script: 02_ssom_fit_by_year.R
# Purpose: Fit single-season occupancy models for ONE year (2008-2023)
# Project: MSc
# Author: Oyede Haleemat Abosede
# ---------------------------------------------------------

# -----------------------------
# 1. Load packages
# -----------------------------
library(unmarked)
library(AICcmodavg)
library(dplyr)
library(readr)
library(AICcmodavg)

# -----------------------------
# 2. User-defined year
# -----------------------------
YEAR <- 2008   # <<< CHANGE THIS ONLY

# -----------------------------
# 3. File paths
# -----------------------------
input_dir  <- "outputData/SSOM/models"
output_tbl <- "outputData/SSOM/tables"
output_best <- "outputData/SSOM/best_models"

if (!dir.exists(output_tbl)) dir.create(output_tbl, recursive = TRUE)
if (!dir.exists(output_best)) dir.create(output_best, recursive = TRUE)

# -----------------------------
# 4. Load prepared SSOM inputs
# -----------------------------
ssom_inputs <- readRDS(
  file.path(input_dir, paste0("ssom_inputs_", YEAR, ".rds"))
)

umf <- ssom_inputs$umf

# -----------------------------
# 5. Fit candidate SSOMs
# -----------------------------
# Detection structure is held constant across all models

M0 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~ 1, umf)

M1 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             Tmax_spline1 + Tmax_spline2 +
             Rf_spline1   + Rf_spline2,
           umf)

M2 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             ndvi_spline1 + ndvi_spline2,
           umf)

M3 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             AG_z,
           umf)

M4 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             UB_z,
           umf)

M5 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             WB_z + NV_z,
           umf)

M6 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             Tmax_spline1 + Tmax_spline2 +
             Rf_spline1   + Rf_spline2 +
             ndvi_spline1 + ndvi_spline2,
           umf)

M7 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             ndvi_spline1 + ndvi_spline2 +
             UB_z,
           umf)

M8 <- occu(~ TotalHours_z + TotalSpp_z + Road_z ~
             Tmax_spline1 + Tmax_spline2 +
             AG_z,
           umf)

# -----------------------------
# 6. Model selection (AIC)
# -----------------------------
model_set <- fitList(
  baseline               = M0,
  climate                = M1,
  productivity           = M2,
  agriculture            = M3,
  urban                  = M4,
  wetland_natural        = M5,
  climate_productivity   = M6,
  urban_productivity     = M7,
  climate_agriculture    = M8
)

aic_results <- modSel(model_set)

aic_table <- aic_results@Full %>%
  select(model, nPars, AIC, delta, AICwt, cumltvWt) %>%
  mutate(Year = YEAR) %>%
  relocate(Year)

# -----------------------------
# 7. Save outputs

# -----------------------------
write_csv(
  aic_table,
  file.path(output_tbl, paste0("SSOM_AIC_", YEAR, ".csv"))
)

best_model_name <- aic_results@Full$model[1]
best_model <- model_set@fits[[best_model_name]]

saveRDS(
  best_model,
  file.path(output_best, paste0("SSOM_best_", YEAR, ".rds"))
)

message("SSOM models fitted and saved for year ", YEAR)
message("Best model: ", best_model_name)
