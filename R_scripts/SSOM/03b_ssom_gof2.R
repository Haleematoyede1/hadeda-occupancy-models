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
library(AICcmodavg)

# -----------------------------
# 2. User-defined settings
# -----------------------------
YEAR  <- 2008        # <<< CHANGE OR LOOP LATER
NSIM  <- 100         # Number of bootstrap simulations
SEED  <- 123         # Reproducibility

set.seed(SEED)

# -----------------------------
# 3. File paths (HPC-safe)
# -----------------------------
base_dir   <- getwd()

input_mod  <- file.path(base_dir, "outputData/SSOM/best_models")
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
  mod       = best_model,
  nsim      = NSIM,
  plot.hist = FALSE
)

summary(best_model)

# -----------------------------
# 6. Extract GOF statistics
# -----------------------------
gof_summary <- tibble::tibble(
  Model_Type      = mb_gof$model.type,
  Test            = "MacKenzie–Bailey GOF",
  ChiSq_Observed  = mb_gof$chi.square,
  ChiSq_Mean_Sim  = mean(mb_gof$t.star),
  ChiSq_SD_Sim    = sd(mb_gof$t.star),
  p_value         = mb_gof$p.value,
  c_hat           = mb_gof$c.hat.est,
  nsim            = mb_gof$nsim
)



gof_summary <- gof_summary %>%
  mutate(
    GOF_Informative = ifelse(
      p_value > 0 & is.finite(c_hat) & c_hat < 10,
      "Yes",
      "No (sparse detection histories)"
    )
  )

# -----------------------------
# 7. Save results
# -----------------------------
write.csv(
  gof_summary,
  file = file.path(
    output_gof,
    paste0("SSOM_GOF_summary_", YEAR, ".csv")
  ),
  row.names = FALSE
)

saveRDS(
  mb_gof,
  file.path(output_gof, paste0("SSOM_GOF_object_", YEAR, ".rds"))
)

message("GOF completed successfully for year ", YEAR)
print(gof_summary)
