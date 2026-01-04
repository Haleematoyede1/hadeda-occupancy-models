# ===============================================================
# 04_ssom_covariate_effects.R
# Purpose: Marginal occupancy effects from best SSOM (year-by-year)
# ===============================================================

library(unmarked)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyverse)

YEAR <- 2023 # THIS WOULDNT WORK FOR SEASONS/YEARS 2010 AND 2014 BEST MODELS WAS DIFFERENT CHECK BELOW FOR CONTINUATION

# -----------------------------
# Paths
# -----------------------------
best_dir    <- "outputData/SSOM/best_models"
input_data  <- "inputData/processedData"
inputs_dir  <- "outputData/SSOM/models"
output_dir  <- "outputData/SSOM/marginal_Effect"
output_dir_plot  <- "outputData/SSOM/marginal_Effect_plot"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load best model + ssom inputs
# -----------------------------
best_model <- readRDS(file.path(best_dir, paste0("SSOM_best_", YEAR, ".rds")))

ssom_inputs <- readRDS(file.path(inputs_dir, paste0("ssom_inputs_", YEAR, ".rds")))
site_covs_train  <- ssom_inputs$site_covs
transform_info   <- ssom_inputs$transform_info

# What occupancy-side covariates does the best model use?
form <- best_model@formula
occ_terms <- all.vars(form[[3]])
message("Occupancy terms in best model: ", paste(occ_terms, collapse = ", "))

# -----------------------------
# Load raw covariates for ranges (untransformed)
# -----------------------------
climate <- read_csv(file.path(input_data, "Hadeda_ClimateCovariates_Data_Annual.csv"),
                    show_col_types = FALSE) %>%
  filter(Year == YEAR)

road <- read_csv(file.path(input_data, "Hadeda_Accessibility_Data.csv"),
                 show_col_types = FALSE) %>%
  rename(Road = RoadDistance)

site_data <- climate %>%
  left_join(road, by = "Pentad") %>%
  drop_na()

# -----------------------------
# Helpers: same transforms as training
# -----------------------------
transform_var <- function(x, info) {
  if (is.null(info)) stop("Missing transform_info for a variable used in plotting.")
  if (isTRUE(info$log)) x <- log(x + 1)
  (x - info$mean) / info$sd
}

sortBase <- function(vec, n.knots = 2) {
  knots <- quantile(unique(vec),
                    probs = seq(0, 1, length.out = n.knots + 2)[-c(1, n.knots + 2)],
                    na.rm = TRUE)
  ZK <- abs(outer(vec, knots, "-"))^3
  OM <- abs(outer(knots, knots, "-"))^3
  sv <- svd(OM)
  sqrtOM <- t(sv$v %*% (t(sv$u) * sqrt(sv$d)))
  t(solve(sqrtOM, t(ZK)))
}

# Build baseline newdata: all occupancy terms held at their means
mean_covs <- site_covs_train %>%
  summarise(across(all_of(occ_terms), ~ mean(.x, na.rm = TRUE)))

make_newdata_base <- function(n) {
  # replicate mean row n times
  as.data.frame(mean_covs[rep(1, n), , drop = FALSE])
}

# Core function: plot effect of a spline variable (2-term spline)
plot_spline_effect <- function(var_raw, spline1, spline2, xlab, title_prefix) {
  
  # Only run if the model uses these spline terms
  if (!all(c(spline1, spline2) %in% occ_terms)) {
    message("Skipping ", var_raw, " effect (", spline1, "/", spline2, " not in model).")
    return(NULL)
  }
  
  x_seq <- seq(min(site_data[[var_raw]], na.rm = TRUE),
               max(site_data[[var_raw]], na.rm = TRUE),
               length.out = 100)
  
  x_z <- transform_var(x_seq, transform_info[[var_raw]])
  base <- scale(sortBase(x_z))
  
  newdata <- make_newdata_base(n = length(x_seq))
  newdata[[spline1]] <- base[, 1]
  newdata[[spline2]] <- base[, 2]
  
  pred <- predict(best_model, newdata = newdata, type = "state")
  
  df <- data.frame(
    x   = x_seq,
    fit = pred$Predicted,
    lo  = pred$lower,
    hi  = pred$upper
  )
  
  # save csv
  write_csv(df, file.path(output_dir, paste0("SSOM_effect_", var_raw, "_", YEAR, ".csv")))
  
  # plot
  p <- ggplot(df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
    geom_line(linewidth = 1) +
    labs(
      x = xlab,
      y = "Occupancy probability",
      title = paste0(title_prefix, " – ", YEAR)
    ) +
    theme_classic()
  
  # save png
  ggsave(file.path(output_dir_plot, paste0("SSOM_effect_", var_raw, "_", YEAR, ".png")),
         plot = p, width = 7, height = 5, dpi = 300)
  
  return(p)
}

# -----------------------------
# Generate effects (only if in model)
# -----------------------------
p_ndvi <- plot_spline_effect(
  var_raw = "ndvi",
  spline1 = "ndvi_spline1",
  spline2 = "ndvi_spline2",
  xlab = "NDVI",
  title_prefix = "Marginal effect of NDVI on occupancy"
)

p_tmax <- plot_spline_effect(
  var_raw = "Tmax",
  spline1 = "Tmax_spline1",
  spline2 = "Tmax_spline2",
  xlab = "Maximum temperature (°C)",
  title_prefix = "Marginal effect of Tmax on occupancy"
)

p_rf <- plot_spline_effect(
  var_raw = "Rf",
  spline1 = "Rf_spline1",
  spline2 = "Rf_spline2",
  xlab = "Rainfall",
  title_prefix = "Marginal effect of rainfall on occupancy"
)

message("Done: covariate effects for year ", YEAR)


#######################################################
#######################################################
#######################################################
# MARGINAL EFFECT PLOTS FOR 2010 AND 2014
#######################################################


YEAR <- 2010 # FOR SEASONS/YEARS 2010 AND 2014 

# -----------------------------
# Paths
# -----------------------------
best_dir    <- "outputData/SSOM/best_models"
input_data  <- "inputData/processedData"
inputs_dir  <- "outputData/SSOM/models"
output_dir  <- "outputData/SSOM/marginal_Effect"
output_dir_plot  <- "outputData/SSOM/marginal_Effect_plot"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load best model + ssom inputs
# -----------------------------
best_model <- readRDS(file.path(best_dir, paste0("SSOM_best_", YEAR, ".rds")))

ssom_inputs <- readRDS(file.path(inputs_dir, paste0("ssom_inputs_", YEAR, ".rds")))
site_covs_train  <- ssom_inputs$site_covs
transform_info   <- ssom_inputs$transform_info

# What occupancy-side covariates does the best model use?
form <- best_model@formula
occ_terms <- all.vars(form[[3]])
message("Occupancy terms in best model: ", paste(occ_terms, collapse = ", "))

# -----------------------------
# Load raw covariates for ranges (untransformed)
# -----------------------------
climate <- read_csv(file.path(input_data, "Hadeda_ClimateCovariates_Data_Annual.csv"),
                    show_col_types = FALSE) %>%
  filter(Year == YEAR)

road <- read_csv(file.path(input_data, "Hadeda_Accessibility_Data.csv"),
                 show_col_types = FALSE) %>%
  rename(Road = RoadDistance)

site_data <- climate %>%
  left_join(road, by = "Pentad") %>%
  drop_na()

# -----------------------------
# Helpers: same transforms as training
# -----------------------------
transform_var <- function(x, info) {
  if (is.null(info)) stop("Missing transform_info for a variable used in plotting.")
  if (isTRUE(info$log)) x <- log(x + 1)
  (x - info$mean) / info$sd
}

sortBase <- function(vec, n.knots = 2) {
  knots <- quantile(unique(vec),
                    probs = seq(0, 1, length.out = n.knots + 2)[-c(1, n.knots + 2)],
                    na.rm = TRUE)
  ZK <- abs(outer(vec, knots, "-"))^3
  OM <- abs(outer(knots, knots, "-"))^3
  sv <- svd(OM)
  sqrtOM <- t(sv$v %*% (t(sv$u) * sqrt(sv$d)))
  t(solve(sqrtOM, t(ZK)))
}

# Build baseline newdata: all occupancy terms held at their means
mean_covs <- site_covs_train %>%
  summarise(across(all_of(occ_terms), ~ mean(.x, na.rm = TRUE)))

make_newdata_base <- function(n) {
  as.data.frame(mean_covs[rep(1, n), , drop = FALSE])
}

# -----------------------------
# Core function: SPLINE effect
# -----------------------------
plot_spline_effect <- function(var_raw, spline1, spline2, xlab, title_prefix) {
  
  if (!all(c(spline1, spline2) %in% occ_terms)) {
    message("Skipping ", var_raw, " effect (", spline1, "/", spline2, " not in model).")
    return(NULL)
  }
  
  x_seq <- seq(min(site_data[[var_raw]], na.rm = TRUE),
               max(site_data[[var_raw]], na.rm = TRUE),
               length.out = 100)
  
  x_z <- transform_var(x_seq, transform_info[[var_raw]])
  base <- scale(sortBase(x_z))
  
  newdata <- make_newdata_base(n = length(x_seq))
  newdata[[spline1]] <- base[, 1]
  newdata[[spline2]] <- base[, 2]
  
  pred <- predict(best_model, newdata = newdata, type = "state")
  
  df <- data.frame(
    x   = x_seq,
    fit = pred$Predicted,
    lo  = pred$lower,
    hi  = pred$upper
  )
  
  write_csv(df, file.path(output_dir, paste0("SSOM_effect_", var_raw, "_", YEAR, ".csv")))
  
  p <- ggplot(df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
    geom_line(linewidth = 1) +
    labs(x = xlab, y = "Occupancy probability",
         title = paste0(title_prefix, " – ", YEAR)) +
    theme_classic()
  
  dir.create(output_dir_plot, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output_dir_plot, paste0("SSOM_effect_", var_raw, "_", YEAR, ".png")),
         plot = p, width = 7, height = 5, dpi = 300)
  
  return(p)
}

# -----------------------------
# Core function: LINEAR effect (e.g., AG_z, UB_z)
# -----------------------------
plot_linear_effect <- function(var_raw, var_z, xlab, title_prefix) {
  
  if (!(var_z %in% occ_terms)) {
    message("Skipping ", var_raw, " effect (", var_z, " not in model).")
    return(NULL)
  }
  
  x_seq <- seq(min(site_data[[var_raw]], na.rm = TRUE),
               max(site_data[[var_raw]], na.rm = TRUE),
               length.out = 100)
  
  # transform raw -> z using training transform_info
  x_z <- transform_var(x_seq, transform_info[[var_raw]])
  
  newdata <- make_newdata_base(n = length(x_seq))
  newdata[[var_z]] <- x_z
  
  pred <- predict(best_model, newdata = newdata, type = "state")
  
  df <- data.frame(
    x   = x_seq,
    fit = pred$Predicted,
    lo  = pred$lower,
    hi  = pred$upper
  )
  
  write_csv(df, file.path(output_dir, paste0("SSOM_effect_", var_raw, "_", YEAR, ".csv")))
  
  p <- ggplot(df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3) +
    geom_line(linewidth = 1) +
    labs(x = xlab, y = "Occupancy probability",
         title = paste0(title_prefix, " – ", YEAR)) +
    theme_classic()
  
  dir.create(output_dir_plot, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output_dir_plot, paste0("SSOM_effect_", var_raw, "_", YEAR, ".png")),
         plot = p, width = 7, height = 5, dpi = 300)
  
  return(p)
}


p_tmax <- plot_spline_effect(
  var_raw = "Tmax",
  spline1 = "Tmax_spline1",
  spline2 = "Tmax_spline2",
  xlab = "Maximum temperature (°C)",
  title_prefix = "Marginal effect of Tmax on occupancy"
)

p_ag <- plot_linear_effect(
  var_raw = "AG",
  var_z   = "AG_z",
  xlab = "Agriculture (%)",
  title_prefix = "Marginal effect of agriculture on occupancy"
)

p_tmax <- plot_spline_effect("Tmax", "Tmax_spline1", "Tmax_spline2",
                             "Maximum temperature (°C)",
                             "Marginal effect of Tmax on occupancy")

p_ag <- plot_linear_effect("AG", "AG_z",
                           "Agriculture (%)",
                           "Marginal effect of agriculture on occupancy")
