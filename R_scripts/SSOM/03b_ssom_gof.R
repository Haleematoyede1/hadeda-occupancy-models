# ===============================================================
# Script: 03b_ssom_gof.R
# Purpose: Goodness-of-fit for best SSOM using parametric bootstrap
# Author: Oyede Haleemat
# ===============================================================

# -----------------------------
# 1. Load packages
# -----------------------------
library(unmarked)
library(dplyr)
library(readr)

# -----------------------------
# 2. User-defined settings
# -----------------------------
YEAR  <- 2008        # <<< CHANGE OR LOOP LATER
NSIM  <- 999         # Number of bootstrap simulations
SEED  <- 123         # Reproducibility

set.seed(SEED)

# -----------------------------
# 3. File paths (HPC-safe)
# -----------------------------
base_dir   <- getwd()

input_mod  <- file.path(base_dir, "outputData/best_models")
input_umf  <- file.path(base_dir, "outputData/SSOM/models")
output_gof <- file.path(base_dir, "outputData/SSOM_GOF")

if (!dir.exists(output_gof)) dir.create(output_gof, recursive = TRUE)

# -----------------------------
# 4. Load best model + UMF
# -----------------------------
best_model <- readRDS(
  file.path(input_mod, paste0("SSOM_best_", YEAR, ".rds"))
)

ssom_inputs <- readRDS(
  file.path(input_umf, paste0("ssom_inputs_", YEAR, ".rds"))
)

umf <- ssom_inputs$umf

# -----------------------------
# 5. MacKenzie–Bailey GOF
# -----------------------------
mb_gof <- mb.gof.test(
  object    = best_model,
  nsim      = NSIM,
  plot.hist = FALSE
)

# -----------------------------
# 6. Extract GOF statistics
# -----------------------------
gof_summary <- tibble(
  Year           = YEAR,
  ChiSquare_obs  = mb_gof$Chi.square,
  ChiSquare_mean = mean(mb_gof$Chi.square.sim),
  p_value        = mb_gof$p.value,
  chat           = mb_gof$Chi.square / mean(mb_gof$Chi.square.sim),
  nsim           = NSIM
)

# -----------------------------
# 7. Save results
# -----------------------------
write_csv(
  gof_summary,
  file.path(output_gof, paste0("SSOM_GOF_", YEAR, ".csv"))
)

saveRDS(
  mb_gof,
  file.path(output_gof, paste0("SSOM_GOF_object_", YEAR, ".rds"))
)

message("GOF completed successfully for year ", YEAR)
print(gof_summary)
