# -----------------------------
# 1. Load packages
# -----------------------------
library(unmarked)
library(tidyverse)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

# -----------------------------
# 2. User-defined year
# -----------------------------
YEAR <- 2008   # <<< CHANGE ONLY THIS

# -----------------------------
# 3. File paths
# -----------------------------
input_models <- "outputData/SSOM/models"
input_data   <- "inputData/processedData"
output_pred  <- "outputData/SSOM/tables"
output_fig   <- "outputData/SSOM/figures"
best_dir   <- "outputData/SSOM/best_models"

dir.create(output_pred, showWarnings = FALSE, recursive = TRUE)
dir.create(output_fig,  showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 4. Load best SSOM for this year
# -----------------------------
best_model <- readRDS(
  file.path(best_dir, paste0("SSOM_best_", YEAR, ".rds"))
)

# -----------------------------
# 5. Load transformation info from training year
# -----------------------------
ssom_inputs <- readRDS(
  file.path(input_models, paste0("ssom_inputs_", YEAR, ".rds"))
)

transform_info <- ssom_inputs$transform_info

# -----------------------------
# 6. Load FULL site-level covariates (ALL pentads)
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


# -----------------------------
# 7. Apply SAME transformations as training
# -----------------------------
transform_var <- function(x, info) {
  if (info$log) x <- log(x + 1)
  (x - info$mean) / info$sd
}

vars_to_transform <- intersect(
  names(transform_info),
  names(site_data)
)

for (v in vars_to_transform) {
  site_data[[paste0(v, "_z")]] <-
    transform_var(site_data[[v]], transform_info[[v]])
}

# -----------------------------
# 8. Rebuild spline bases 
# -----------------------------
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

spline_vars <- c("Tmax", "Tmin", "Rf", "Soil", "ndvi")

for (v in spline_vars) {
  
  z <- site_data[[paste0(v, "_z")]]
  ok <- is.finite(z)
  
  # Guard against poor variation
  if (sum(ok) < 10 || length(unique(z[ok])) < 5) {
    site_data[[paste0(v, "_spline1")]] <- NA
    site_data[[paste0(v, "_spline2")]] <- NA
    next
  }
  
  base <- scale(sortBase(z[ok]))
  
  site_data[[paste0(v, "_spline1")]] <- NA
  site_data[[paste0(v, "_spline2")]] <- NA
  
  site_data[[paste0(v, "_spline1")]][ok] <- base[, 1]
  site_data[[paste0(v, "_spline2")]][ok] <- base[, 2]
}

# -----------------------------
# 9. Build prediction covariate matrix
# -----------------------------
site_covs <- site_data %>%
  select(
    Pentad,
    Tmax_spline1, Tmax_spline2,
    Tmin_spline1, Tmin_spline2,
    Rf_spline1,   Rf_spline2,
    Soil_spline1, Soil_spline2,
    ndvi_spline1, ndvi_spline2,
    NV_z, AG_z, UB_z, WB_z, OTH_z,
    Road_z
  ) %>%
  distinct(Pentad, .keep_all = TRUE) %>%
  column_to_rownames("Pentad")


# -----------------------------
# 10. Predict occupancy for ALL pentads
# -----------------------------
psi_hat <- predict(
  best_model,
  newdata = site_covs,
  type = "state"
)

occ_results <- data.frame(
  Pentad         = rownames(site_covs),
  occupancy_prob = psi_hat$Predicted,
  lower_ci       = psi_hat$lower,
  upper_ci       = psi_hat$upper,
  Year           = YEAR
)

# -----------------------------
# 11. Create coordinates from pentad ID
# -----------------------------
occ_results <- occ_results %>%
  mutate(
    Latitude  = -(as.numeric(substr(Pentad, 1, 2)) +
                    (as.numeric(substr(Pentad, 3, 4)) + 2.5) / 60),
    Longitude =  (as.numeric(substr(Pentad, 6, 7)) +
                    (as.numeric(substr(Pentad, 8, 9)) + 2.5) / 60)
  )

# -----------------------------
# 12. Save prediction table
# -----------------------------
write_csv(
  occ_results,
  file.path(output_pred,
            paste0("SSOM_Occupancy_Predictions_", YEAR, ".csv"))
)

# -----------------------------
# 13. Plot predicted occupancy map
# -----------------------------
ggplot(occ_results,
       aes(x = Longitude, y = Latitude, fill = occupancy_prob)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = brewer.pal(7, "Spectral"),
    name = "Occupancy\nProbability"
  ) +
  coord_equal() +
  labs(
    title = paste("Predicted Occupancy Probability –", YEAR),
    x = "Longitude", y = "Latitude"
  ) +
  theme_void()

ggsave(
  file.path(output_fig,
            paste0("SSOM_Occupancy_Map_", YEAR, ".png")),
  width = 8, height = 6, dpi = 300
)

message("SSOM occupancy predictions completed for year ", YEAR)
