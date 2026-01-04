# ===============================================================
# 04_ssom_covariate_effects.R
# Purpose: Marginal occupancy effects from best SSOM
# ===============================================================

library(unmarked)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)

YEAR <- 2008

# -----------------------------
# Paths
# -----------------------------
best_dir   <- "outputData/SSOM/best_models"
input_data <- "inputData/processedData"
output_dir <- "outputData/SSOM/figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load best model
# -----------------------------
best_model <- readRDS(
  file.path(best_dir, paste0("SSOM_best_", YEAR, ".rds"))
)

best_model@formula
# -----------------------------
# Load site covariates (training year)
# -----------------------------
climate <- read_csv(
  file.path(input_data, "Hadeda_ClimateCovariates_Data_Annual.csv"),
  show_col_types = FALSE
) %>%
  filter(Year == YEAR)

road <- read_csv(
  file.path(input_data, "Hadeda_Accessibility_Data.csv"),
  show_col_types = FALSE
) %>%
  rename(Road = RoadDistance)

site_data <- climate %>%
  left_join(road, by = "Pentad") %>%
  drop_na()

# ----------------------------
# Required functions
# ----------------------------

ssom_inputs <- readRDS(
  file.path("outputData/SSOM/models",
            paste0("ssom_inputs_", YEAR, ".rds"))
)

# Get the mean for all variables in the best model
site_covariate <- ssom_inputs[["site_covs"]]####
mean_covs <- site_covariate %>%
  summarise(
    Tmax_spline1 = mean(Tmax_spline1, na.rm = TRUE),
    Tmax_spline2 = mean(Tmax_spline2, na.rm = TRUE),
    Rf_spline1   = mean(Rf_spline1, na.rm = TRUE),
    Rf_spline2   = mean(Rf_spline2, na.rm = TRUE),
    ndvi_spline1 = mean(ndvi_spline1, na.rm = TRUE),
    ndvi_spline2 = mean(ndvi_spline2, na.rm = TRUE)
  )


#scale transformation
transform_info <- ssom_inputs[["transform_info"]]
# Log + z transform using stored training info
transform_var <- function(x, info) {
  if (info$log) x <- log(x + 1)
  (x - info$mean) / info$sd
}

# spline function
sortBase <- function(vec, n.knots = 2) {
  knots <- quantile(
    unique(vec),
    probs = seq(0, 1, length.out = n.knots + 2)[-c(1, n.knots + 2)],
    na.rm = TRUE
  )
  ZK <- abs(outer(vec, knots, "-"))^3
  OM <- abs(outer(knots, knots, "-"))^3
  sv <- svd(OM)
  sqrtOM <- t(sv$v %*% (t(sv$u) * sqrt(sv$d)))
  t(solve(sqrtOM, t(ZK)))
}

# -----------------------------
# NDVI Marginal covariate effect plot
# -----------------------------
# -----------------------------
# NDVI sequence for plot ranges
# -----------------------------
ndvi_seq <- seq(
  min(site_data$ndvi, na.rm = TRUE),
  max(site_data$ndvi, na.rm = TRUE),
  length.out = 100
)

# -----------------------------
# Apply SAME transformations as training
# -----------------------------
ndvi_z <- transform_var(ndvi_seq, transform_info[["ndvi"]])

# -----------------------------
# Apply SAME Spline transformations as training
# -----------------------------
ndvi_base <- scale(sortBase(ndvi_z))

# -----------------------------
# Hold other covariates constant
# -----------------------------
newdata <-  data.frame(
  ndvi_spline1 = ndvi_base[,1],
  ndvi_spline2 = ndvi_base[,2],
  Tmax_spline1 = mean_covs$Tmax_spline1,
  Tmax_spline2 = mean_covs$Tmax_spline2,
  Rf_spline1   = mean_covs$Rf_spline1,
  Rf_spline2   = mean_covs$Rf_spline2
)

# -----------------------------
# Predict occupancy
# -----------------------------
pred <- predict(
  best_model,
  newdata = newdata,
  type = "state"
)


plot_df <- data.frame(
  ndvi = ndvi_seq,
  fit  = pred$Predicted,
  lo   = pred$lower,
  hi   = pred$upper
)

write_csv(
  plot_df,
  file.path(output_dir,
            paste0("SSOM_effect_NDVI_", YEAR, ".csv"))
)

# -----------------------------
# Plot
# -----------------------------
p1 = ggplot(plot_df, aes(x = ndvi, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(size = 1) +
  labs(
    x = "NDVI",
    y = "Occupancy probability",
    title = paste("Effect of NDVI on Hadeda occupancy –", YEAR)
  ) +
  theme_minimal()

# -----------------------------
# Maximum Temperature Marginal covariate effect plot
# -----------------------------
# -----------------------------
# Tmax sequence for plot ranges
# -----------------------------
Tmax_seq <- seq(
  min(site_data$Tmax, na.rm = TRUE),
  max(site_data$Tmax, na.rm = TRUE),
  length.out = 100
)

# -----------------------------
# Apply SAME transformations as training
# -----------------------------
Tmax_z <- transform_var(Tmax_seq, transform_info[["Tmax"]])

# -----------------------------
# Apply SAME Spline transformations as training
# -----------------------------
Tmax_base <- scale(sortBase(Tmax_z))

# -----------------------------
# Hold other covariates constant
# -----------------------------
newdata <-  data.frame(
  Tmax_spline1 = Tmax_base[,1],
  Tmax_spline2 = Tmax_base[,2],
  ndvi_spline1 = mean_covs$ndvi_spline1,
  ndvi_spline2 = mean_covs$ndvi_spline2,
  Rf_spline1   = mean_covs$Rf_spline1,
  Rf_spline2   = mean_covs$Rf_spline2
)

# -----------------------------
# Predict occupancy
# -----------------------------
pred <- predict(
  best_model,
  newdata = newdata,
  type = "state"
)


plot_df <- data.frame(
  Tmax = Tmax_seq,
  fit  = pred$Predicted,
  lo   = pred$lower,
  hi   = pred$upper
)

write_csv(
  plot_df,
  file.path(output_dir,
            paste0("SSOM_effect_Tmax_", YEAR, ".csv"))
)

# -----------------------------
# Plot
# -----------------------------
p2 = ggplot(plot_df, aes(x = Tmax, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(size = 1) +
  labs(
    x = "Tmax",
    y = "Occupancy probability",
    title = paste("Effect of Temperature on Hadeda occupancy –", YEAR)
  ) +
  theme_minimal()

# -----------------------------
# Rainfall Marginal covariate effect plot
# -----------------------------
# -----------------------------
# Rf sequence for plot ranges
# -----------------------------
Rf_seq <- seq(
  min(site_data$Rf, na.rm = TRUE),
  max(site_data$Rf, na.rm = TRUE),
  length.out = 100
)

# -----------------------------
# Apply SAME transformations as training
# -----------------------------
Rf_z <- transform_var(Rf_seq, transform_info[["Rf"]])

# -----------------------------
# Apply SAME Spline transformations as training
# -----------------------------
Rf_base <- scale(sortBase(Rf_z))

# -----------------------------
# Hold other covariates constant
# -----------------------------
newdata <-  data.frame(
  Rf_spline1 = Rf_base[,1],
  Rf_spline2 = Rf_base[,2],
  ndvi_spline1 = mean_covs$ndvi_spline1,
  ndvi_spline2 = mean_covs$ndvi_spline2,
  Tmax_spline1   = mean_covs$Tmax_spline1,
  Tmax_spline2   = mean_covs$Tmax_spline2
)

# -----------------------------
# Predict occupancy
# -----------------------------
pred <- predict(
  best_model,
  newdata = newdata,
  type = "state"
)


plot_df <- data.frame(
  Rf = Rf_seq,
  fit  = pred$Predicted,
  lo   = pred$lower,
  hi   = pred$upper
)

write_csv(
  plot_df,
  file.path(output_dir,
            paste0("SSOM_effect_Rf_", YEAR, ".csv"))
)

# -----------------------------
# Plot
# -----------------------------
p3 = ggplot(plot_df, aes(x = Rf, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
  geom_line(size = 1) +
  labs(
    x = "Rf",
    y = "Occupancy probability",
    title = paste("Effect of Rainfall on Hadeda occupancy –", YEAR)
  ) +
  theme_minimal()

##########################################################################

